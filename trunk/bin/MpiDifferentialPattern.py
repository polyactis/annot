#!/usr/bin/env mpipython
"""
Usage: MpiDifferentialPattern -k SCHEMA -g -i -o -s  [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-g ...,	gim(gene incidence matrix) inputfile
	-i ...,	the input file, dataset signature output by fim
	-l ...,	lm_bit of setting1, ('00001', default)
	-a ...,	acc_cutoff of setting1, (0.6 default)
	-o ...,	pic_output_dir.(output file and graph pictures)
	-p ...,	p value cutoff(0.01, default)
	-s ...,	dataset_signature, like '1,3,30'
	-v ...,	message size(10,000,000, default)
	-b,	debug version.
	-r,	enable report flag, WATCH: enable graph drawing
	-h, --help	Display the usage infomation.
	
Examples:
	mpirun -np 20 -machinefile ~/hostfile /usr/bin/mpipython ~/script/annot/bin/MpiDifferentialPattern.py
	-k -g -i -o -s
	
Description:
	Find patterns whose recurrence shows difference between the given
	  dataset_signature and the rest.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import sys, os, getopt, csv, math, cPickle
from Scientific import MPI
from codense.common import mpi_synchronize, db_connect, output_node, draw_pattern,form_schema_tables,\
	computing_node, input_node, get_gene_id2gene_no, get_gene_no2gene_id, get_gene_no2go_no
from sets import Set
from rpy import r
from MpiFromDatasetSignatureToPattern import encodeOccurrenceBv, encodeOccurrence
from cluster_info import cluster_info


class MpiDifferentialPattern:
	def __init__(self,hostname='zhoudb', dbname='graphdb', schema=None, gim_inputfname=None, fname=None,\
		lm_bit='00001', acc_cutoff=0.6, p_value_cut_off=0.01, pic_output_dir=None, dataset_signature=[],\
		message_size=10000000, debug=0, report=0):
		"""
		10-22-05
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.gim_inputfname = gim_inputfname
		self.fname = fname
		self.lm_bit = lm_bit
		self.acc_cutoff = float(acc_cutoff)
		self.p_value_cut_off = float(p_value_cut_off)
		self.pic_output_dir = pic_output_dir
		self.dataset_signature_set = Set(dataset_signature)
		self.message_size = int(message_size)
		self.debug = int(debug)
		self.report = int(report)
	
	def get_gene2enc_array(self, gim_inputfname, gene_id2no):
		sys.stderr.write("Getting gene2enc_array...\n")
		reader = csv.reader(open(gim_inputfname), delimiter='\t')
		gene2enc_array = {}
		for row in reader:
			no_of_occurrences, occ_array, gene_id = row[0], row[1:-1], row[-1]
			if gene_id in gene_id2no:
				gene_no = gene_id2no[gene_id]
				occ_array = map(int, occ_array)
				gene2enc_array[gene_no] = encodeOccurrenceBv(occ_array)
		sys.stderr.write("End getting gene2enc_array.\n")
		return gene2enc_array

	def get_effective_vertex_set(self, vertex_list, gene2enc_array, which_dataset):
		effective_vertex_set = Set()
		enc_dataset = encodeOccurrence([which_dataset])
		for vertex in vertex_list:
			if gene2enc_array[vertex]&enc_dataset == enc_dataset:
				effective_vertex_set.add(vertex)
		return effective_vertex_set
	
	def cal_effective_rec_array(self, edge_set, original_rec_array, vertex_list, gene2enc_array, debug=0):
		"""
		11-03-05
			based on the present vertices in each dataset, calculate no_of_effective_possible_edges,
			adjusted occurrence = no_of_real_edges_int/no_of_effective_possible_edges
		"""
		effective_rec_array = []
		no_of_edges = len(edge_set)
		for i in range(len(original_rec_array)):
			no_of_real_edges_float = original_rec_array[i]*no_of_edges
			no_of_real_edges_int = int(no_of_real_edges_float)	#need to compensate the loss because of integer cast
			if (no_of_real_edges_float-no_of_real_edges_int)>0.5:
				no_of_real_edges_int += 1
			effective_vertex_set = self.get_effective_vertex_set(vertex_list, gene2enc_array, i+1)	#WATCH i+1
			no_of_effective_possible_edges = 0
			for edge in edge_set:
				if edge[0] in effective_vertex_set and edge[1] in effective_vertex_set:
					no_of_effective_possible_edges += 1
			"""
			if debug:
				print "which_dataset",i+1
				print "effective_vertex_set",effective_vertex_set
				print "no_of_real_edges_int",no_of_real_edges_int
				print "no_of_effective_possible_edges",no_of_effective_possible_edges
			"""
			if no_of_effective_possible_edges:
				effective_rec_array.append(no_of_real_edges_int/float(no_of_effective_possible_edges))
			else:
				effective_rec_array.append(1)
		return effective_rec_array
	
	def get_t_test_p_value(self, effective_rec_array, interesting_dataset_set, debug=0):
		"""
		11-03-05 partition the effective_rec_array into non-interesting and interesting based on the interesting_dataset_set
			do the t-test, get the p-value
		"""
		non_interesting_rec_array = []
		interesting_rec_array = []
		for i in range(len(effective_rec_array)):
			if i+1 in interesting_dataset_set:
				interesting_rec_array.append(effective_rec_array[i])
			else:
				non_interesting_rec_array.append(effective_rec_array[i])
		"""
		if debug:
			print "interesting_rec_array",interesting_rec_array
			print "non_interesting_rec_arrray",non_interesting_rec_array
		"""
		result = r.t_test(non_interesting_rec_array, interesting_rec_array)
		return result['p.value']
	
	def computing_node_handler(self, communicator, data, parameter_list):
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		gene2enc_array, dataset_signature_set, p_value_cut_off = parameter_list
		data = cPickle.loads(data)
		good_patterns = []
		for row in data:
			id, vertex_set_string, edge_set_string, recurrence_array_string, go_no_list = row
			
			vertex_set = vertex_set_string[1:-1].split(',')
			vertex_set = map(int, vertex_set)
			edge_set = edge_set_string[2:-2].split('},{')
			for i in range(len(edge_set)):
				edge_set[i] = edge_set[i].split(',')
				edge_set[i] = map(int, edge_set[i])
			recurrence_array = recurrence_array_string[1:-1].split(',')
			recurrence_array = map(float, recurrence_array)
			effective_rec_array = self.cal_effective_rec_array(edge_set, recurrence_array, vertex_set, gene2enc_array)
			p_value = self.get_t_test_p_value(effective_rec_array, dataset_signature_set)
			if p_value<=p_value_cut_off:
				go_no_list = go_no_list[1:-1].split(',')
				go_no_list = map(int, go_no_list)
				good_patterns.append([id, vertex_set, edge_set, effective_rec_array, p_value, go_no_list])
		sys.stderr.write("Node no.%s done with %s good_patterns.\n"%(node_rank, len(good_patterns)))
		return good_patterns
		
	def output_node_handler(self, communicator, parameter_list, data):
		"""
		11-03-05
			called by common.output_node()
		11-04-05 graph output is only turned if self.report
		"""
		pic_output_dir, cluster_info_instance, gene_no2id, gene_no2go_no, writer = parameter_list
		good_patterns = cPickle.loads(data)
		for row in good_patterns:
			writer.writerow(row)
			id, vertex_set, edge_set, effective_rec_array, p_value, go_no_list = row
			if self.report:
				for go_no in go_no_list:
					graphFname = os.path.join(pic_output_dir, '%s_%s.png'%(go_no, id))
					draw_pattern(vertex_set, edge_set, go_no, gene_no2id, gene_no2go_no, cluster_info_instance, graphFname)
		
		
	def run(self):
		communicator = MPI.world.duplicate()
		node_rank = communicator.rank	
		free_computing_nodes = range(1,communicator.size-1)	#exclude the last node
		if node_rank == 0:
			(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
			schema_instance = form_schema_tables(self.fname, self.acc_cutoff, self.lm_bit)
			gene_id2no = get_gene_id2gene_no(curs)
			gene2enc_array = self.get_gene2enc_array(self.gim_inputfname, gene_id2no)
			gene2enc_array_pickle = cPickle.dumps(gene2enc_array, -1)
			
			gene_no2id = get_gene_no2gene_id(curs)
			gene_no2go_no = get_gene_no2go_no(curs)
			gene_no2id_pickle = cPickle.dumps(gene_no2id, -1)
			gene_no2go_no_pickle = cPickle.dumps(gene_no2go_no, -1)
			for node in free_computing_nodes:	#send it to the computing_node
				communicator.send(gene2enc_array_pickle, node, 0)
			
			communicator.send(gene_no2id_pickle, communicator.size-1, 0)
			communicator.send(gene_no2go_no_pickle, communicator.size-1, 0)
		elif node_rank in free_computing_nodes:
			data, source, tag = communicator.receiveString(0, 0)
			gene2enc_array = cPickle.loads(data)	#take the data
		elif node_rank==communicator.size-1:
			schema_instance = form_schema_tables(self.fname, self.acc_cutoff, self.lm_bit)
			data, source, tag = communicator.receiveString(0, 0)
			gene_no2id = cPickle.loads(data)
			data, source, tag = communicator.receiveString(0, 0)
			gene_no2go_no = cPickle.loads(data)
			
		mpi_synchronize(communicator)
		if node_rank == 0:
			curs.execute("DECLARE crs CURSOR FOR SELECT p.id, p.vertex_set, p.edge_set, p.recurrence_array,\
			g.go_no_list from %s p, %s g where g.mcl_id=p.id"%(schema_instance.pattern_table, schema_instance.good_cluster_table))
			input_node(communicator, curs, free_computing_nodes, self.message_size, self.report)
		elif node_rank in free_computing_nodes:
			parameter_list = [gene2enc_array, self.dataset_signature_set, self.p_value_cut_off]
			computing_node(communicator, parameter_list, self.computing_node_handler, report=self.report)
		elif node_rank==communicator.size-1:
			if not os.path.isdir(self.pic_output_dir):
				os.makedirs(self.pic_output_dir)
			cluster_info_instance = cluster_info()
			ofname = os.path.join(self.pic_output_dir, '%s_p%s'%(schema_instance.good_cluster_table, self.p_value_cut_off))
			writer = csv.writer(open(ofname, 'w'), delimiter='\t')
			parameter_list = [self.pic_output_dir, cluster_info_instance, gene_no2id, gene_no2go_no, writer]
			output_node(communicator, free_computing_nodes, parameter_list, self.output_node_handler, self.report)
			del writer

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:g:i:l:a:o:p:s:v:br", ["help", \
			"hostname=", "dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	gim_inputfname = None
	fname = None
	lm_bit = '00001'
	acc_cutoff = 0.6
	pic_output_dir = None
	p_value_cut_off = 0.01
	dataset_signature = []
	message_size = 10000000
	debug = 0
	report = 0
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-z", "--hostname"):
			hostname = arg
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-k", "--schema"):
			schema = arg
		elif opt in ("-g"):
			gim_inputfname = arg
		elif opt in ("-i"):
			fname = arg
		elif opt in ("-l"):
			lm_bit = arg
		elif opt in ("-a"):
			acc_cutoff = float(arg)
		elif opt in ("-o"):
			pic_output_dir = arg
		elif opt in ("-p"):
			p_value_cut_off = float(arg)
		elif opt in ("-s"):
			dataset_signature = arg.split(',')
			dataset_signature = map(int, dataset_signature)
		elif opt in ("-v"):
			message_size = int(arg)
		elif opt in ("-b"):
			debug = 1
		elif opt in ("-r"):
			report = 1
	if schema and gim_inputfname and fname and pic_output_dir and dataset_signature:
		instance = MpiDifferentialPattern(hostname, dbname, schema, gim_inputfname, fname,\
			lm_bit, acc_cutoff, p_value_cut_off, pic_output_dir, dataset_signature, message_size, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
