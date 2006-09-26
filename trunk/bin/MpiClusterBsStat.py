#!/usr/bin/env mpipython
"""
Usage: MpiClusterBsStat.py -k SCHEMA -g GOOD_CLUSTER_TABLE
	-v sig_vector_fname -l CLUSTER_BS_TABLE [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-g ...,	GOOD_CLUSTER_TABLE
	-l ...,	CLUSTER_BS_TABLE
	-s ...,	no of clusters per transmission, 2000(default)
	-a ...,	min ratio of #associated genes for one bs_no vs cluster size, 1/3(default)
	-t ...,	top_number of scores to be kept, 5(default)
	-e ...,	degree_cut_off, 0.3(default), for fuzzyDense
	-p ...,	p_value_cut_off, 0.001(default)
	-o ...,	output_file, fuzzyDense.out(default)
	-x ...,	tax_id, 4932(default, yeast), for get_gene_id2gene_symbol
	-v ...,	sig_vector_fname
	-f,	run fuzzyDense
	-n,	CLUSTER_BS_TABLE is new
	-c,	commit the database transaction
	-b,	debug version.
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	mpirun -np 20 -machinefile ~/hostfile /usr/bin/mpipython ~/script/annot/bin/MpiClusterBsStat.py
	-k mm_fim_97 -g good_clusters -l cluster_bs -v
	
Description:
	Program to do binding site enrichment analysis for a module.
	
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import sys, os, getopt, csv, math, Numeric, cPickle
from Scientific import MPI
from codense.common import mpi_synchronize, db_connect, output_node, get_gene_id2gene_symbol, \
	dict_map, get_dataset_no2desc, input_node, computing_node
from sets import Set
from TF_functions import cluster_bs_analysis
from MpiCrackSplat import MpiCrackSplat	#12-18-05	for fill_edge2encodedOccurrence()
from fuzzyDense import fuzzyDense	#12-18-05	for fuzzyDense()
if sys.version_info[:2] < (2, 3):       #python2.2 or lower needs some extra
	from python2_3 import *
	
class MpiClusterBsStat:
	"""
	09-19-05
	"""
	def __init__(self,hostname='zhoudb', dbname='graphdb', schema=None, good_cluster_table=None,\
		cluster_bs_table=None,  size=2000, ratio_cutoff=1/3.0, top_number=5, degree_cut_off=0.3, p_value_cut_off=0.001, \
		output_file='fuzzyDense.out', tax_id=4932, sig_vector_fname=None, fuzzyDense_flag=0, new_table=0, \
		commit=0, debug=0, report=0):
		"""
		12-15-05
			add p_value_cut_off
		12-20-05
			add output_file and tax_id, degree_cut_off
		2006-09-21 add fuzzyDense_flag
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.good_cluster_table = good_cluster_table
		self.cluster_bs_table = cluster_bs_table
		self.size = int(size)
		self.ratio_cutoff = float(ratio_cutoff)
		self.top_number = int(top_number)
		self.degree_cut_off = float(degree_cut_off)
		self.p_value_cut_off = float(p_value_cut_off)
		self.output_file = output_file
		self.tax_id = int(tax_id)
		self.sig_vector_fname = sig_vector_fname
		self.fuzzyDense_flag = int(fuzzyDense_flag)
		self.new_table = int(new_table)
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
	
	def get_gene_no2bs_no_block(self, curs):
		"""
		09-19-05
			input_node takes care of this
			format:
				[gene_no1, bs_no11, bs_no12, ..., -1, gene_no2, gene_no21, ..., -1]
			-1 is the dilimiter
		"""
		sys.stderr.write("Getting gene_no2bs_no_block...\n")
		gene_no2bs_no_block = []
		curs.execute("DECLARE crs0 CURSOR FOR select gene_no, trans_fac from gene")
		curs.execute("fetch 5000 from crs0")
		rows = curs.fetchall()
		while rows:
			for row in rows:
				gene_no, trans_fac_ls = row
				if trans_fac_ls:	#throw away gene_nos with no trans_fac
					trans_fac_ls = trans_fac_ls[1:-1].split(',')
					trans_fac_ls = map(int, trans_fac_ls)
					gene_no2bs_no_block.append(gene_no)
					for bs_no in trans_fac_ls:
						gene_no2bs_no_block.append(bs_no)
					gene_no2bs_no_block.append(-1)	#delimiter
			curs.execute("fetch 5000 from crs0")
			rows = curs.fetchall()
		gene_no2bs_no_block = Numeric.array(gene_no2bs_no_block)
		sys.stderr.write("Done.\n")
		return gene_no2bs_no_block
	
	def construct_two_dicts(self, node_rank, gene_no2bs_no_block):
		"""
		09-19-05
			computing nodes need to have two dicts
		"""
		sys.stderr.write("Node no.%s constructing two dicts...\n"%(node_rank))
		gene_no_and_bs_no_list = []
		gene_no2bs_no_set = {}
		bs_no2gene_no_set = {}
		for no in gene_no2bs_no_block:
			if no==-1:
				gene_no = gene_no_and_bs_no_list[0]
				for bs_no in gene_no_and_bs_no_list[1:]:
					if gene_no not in gene_no2bs_no_set:
						gene_no2bs_no_set[gene_no] = Set()
					gene_no2bs_no_set[gene_no].add(bs_no)
					if bs_no not in bs_no2gene_no_set:
						bs_no2gene_no_set[bs_no] = Set()
					bs_no2gene_no_set[bs_no].add(gene_no)
				gene_no_and_bs_no_list = []
			else:
				gene_no_and_bs_no_list.append(no)
		sys.stderr.write("Node no.%s done with constructing two dicts.\n"%(node_rank))
		return gene_no2bs_no_set, bs_no2gene_no_set
	
	def return_total_vertex_set(self, curs, good_cluster_table):
		"""
		04-06-06
		"""
		sys.stderr.write("Getting total vertex_set...\n")
		curs.execute("DECLARE crs CURSOR FOR select vertex_set from %s "%(good_cluster_table))
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		total_vertex_set = Set()
		while rows:
			for row in rows:
				vertex_set = row[0][1:-1].split(',')
				for vertex in vertex_set:
					total_vertex_set.add(int(vertex))
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		curs.execute("close crs")
		sys.stderr.write("Done with total vertex_set.\n")
		return total_vertex_set
	
	def fill_edge2encodedOccurrence(self, sig_vector_fname, min_sup, max_sup, total_vertex_set=None):
		"""
		04-04-06
		"""
		sys.stderr.write("Getting edge2encodedOccurrence...\n")
		from MpiFromDatasetSignatureToPattern import encodeOccurrenceBv
		edge2encodedOccurrence = {}
		reader = csv.reader(open(sig_vector_fname), delimiter='\t')
		no_of_datasets = 0
		counter = 0
		for row in reader:
			edge = row[:2]
			edge = map(int, edge)
			#04-06-06 any vertex of the edge doesn't appear in total_vertex_set, skip the edge
			if total_vertex_set and (edge[0] not in total_vertex_set or edge[1] not in total_vertex_set):
				continue
			edge.sort()	#04-06-06 in ascending order
			sig_vector = row[2:]
			sig_vector = map(int, sig_vector)
			if no_of_datasets==0:
				no_of_datasets = len(sig_vector)
			if sum(sig_vector)>=min_sup and sum(sig_vector)<=max_sup:
				edge2encodedOccurrence[tuple(edge)] = encodeOccurrenceBv(sig_vector)
		sys.stderr.write("Done.\n")
		del reader
		return edge2encodedOccurrence, no_of_datasets
	
	def fetch_cluster_block(self, curs, size):
		"""
		09-19-05
		12-18-05
			add recurrence_array
			
			format:
				[mcl_id1, gene_no, gene_no, ..., -1, recurrence, recurrence, ..., -2, mcl_id2, gene_no, gene_no, ...]
		"""
		curs.execute("fetch %s from crs"%(size))
		rows = curs.fetchall()
		cluster_block = []
		for row in rows:
			mcl_id, vertex_list, recurrence_array = row
			vertex_list = vertex_list[1:-1].split(',')
			vertex_list = map(int, vertex_list)
			recurrence_array = recurrence_array[1:-1].split(',')
			recurrence_array = map(int, recurrence_array)
			
			cluster_block.append(mcl_id)
			for vertex in vertex_list:
				cluster_block.append(vertex)
			cluster_block.append(-1)
			
			for recurrence in recurrence_array:
				cluster_block.append(recurrence)
			cluster_block.append(-2)
			
		cluster_block = Numeric.array(cluster_block)
		return cluster_block
	
	def computing_node_handler(self, communicator, data, parameter_list):
		"""
		12-20-05
			called common's computing_node()
		2006-09-21 add fuzzyDense_flag
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		gene_no2bs_no_set, bs_no2gene_no_set, ratio_cutoff, \
		top_number, p_value_cut_off, fuzzyDense_instance, degree_cut_off, fuzzyDense_flag = parameter_list
		result = []
		for row in data:
			id, vertex_set, recurrence_array = row
			vertex_set = vertex_set[1:-1].split(',')
			vertex_set = map(int, vertex_set)
			if fuzzyDense_flag:
				recurrence_array = recurrence_array[1:-1].split(',')
				recurrence_array = map(float, recurrence_array)
				core_vertex_list, on_dataset_index_ls =  fuzzyDense_instance.get_core_vertex_set(vertex_set, recurrence_array, degree_cut_off)
				if core_vertex_list:
					ls_to_return = cluster_bs_analysis(core_vertex_list, gene_no2bs_no_set, bs_no2gene_no_set, ratio_cutoff, \
						top_number, p_value_cut_off)
					if ls_to_return:
						result.append([id, core_vertex_list, on_dataset_index_ls, ls_to_return])
			else:
				ls_to_return = cluster_bs_analysis(vertex_set, gene_no2bs_no_set, bs_no2gene_no_set, ratio_cutoff, \
					top_number, p_value_cut_off)
				if ls_to_return:
					result.append([id, vertex_set, [0], ls_to_return])	#05-31-06, on_dataset_index_ls is faked by [0]
		sys.stderr.write("Node no.%s done.\n"%node_rank)
		return result
	
	
	def create_cluster_bs_table(self, curs, cluster_bs_table):
		"""
		04-06-06
			add core_vertex_ls and on_dataset_index_ls
		"""
		sys.stderr.write("Creating %s...\n"%cluster_bs_table)
		curs.execute("create table %s(\
			id	serial primary key,\
			mcl_id	integer,\
			core_vertex_ls	integer[],\
			on_dataset_index_ls	integer[],\
			bs_no_list	integer[],\
			gene_no_list	integer[],\
			global_ratio	float,\
			local_ratio	float,\
			expected_ratio	float,\
			unknown_ratio	float,\
			score	float,\
			score_type	integer)"%(cluster_bs_table))
		sys.stderr.write("%s created.\n"%cluster_bs_table)
	
	def submit_cluster_bs_table(self, communicator, parameter_list, data):
		"""
		09-20-05
			data is Float
		04-06-06
			the format of 'data' is changed
			add core_vertex_ls, on_dataset_index_ls
		"""
		curs, cluster_bs_table = parameter_list
		data = cPickle.loads(data)
		for row in data:
			id, core_vertex_ls, on_dataset_index_ls, ls_to_return = row
			for tfbs_row in ls_to_return:
				score, score_type, bs_no_list, target_gene_no_list, global_ratio, local_ratio, expected_ratio, unknown_ratio = tfbs_row
				curs.execute("insert into %s(mcl_id, core_vertex_ls, on_dataset_index_ls, bs_no_list, gene_no_list, global_ratio, local_ratio, \
					expected_ratio, unknown_ratio, score, score_type) values(%s, ARRAY%s, ARRAY%s, '{%s}', '{%s}', %s, %s,\
					%s, %s, %s, %s)"%(cluster_bs_table, id, repr(core_vertex_ls), repr(on_dataset_index_ls), repr(bs_no_list)[1:-1], repr(target_gene_no_list)[1:-1],\
					global_ratio, local_ratio, expected_ratio, unknown_ratio, score, score_type))
	
	def output_cluster_bs_data(self, communicator, parameter_list, data):
		"""
		12-20-05 for darwin output
		
			out:=[
			[id, [fuzzyDense gene list], { [ [TF gene list1], [TF target gene list1], p-value], [ [TF gene list2], [TF target gene list2], p-value], ...}, \
				{[dataset_no1, description], [dataset_no2, description], ... }  ],
			[...],
			...
			[]]:
			
		"""
		outf, gene_id2symbol, dataset_no2desc = parameter_list
		data = cPickle.loads(data)
		for row in data:
			id, core_vertex_ls, on_dataset_index_ls, ls_to_return = row
			#prepare the dataset_no_desc_ls
			dataset_no_desc_ls = []
			for dataset_index in on_dataset_index_ls:
				dataset_no = dataset_index +1
				dataset_no_desc_ls.append([dataset_no, dataset_no2desc[dataset_no]])
			#prepare the tfbs_row_darwin_ls
			tfbs_row_darwin_ls = []
			for tfbs_row in ls_to_return:
				score, score_type, bs_no_list, target_gene_no_list, global_ratio, local_ratio, expected_ratio, unknown_ratio = tfbs_row
				bs_no_symbol_list = dict_map(gene_id2symbol, bs_no_list)
				target_gene_no_symbol_list = dict_map(gene_id2symbol, target_gene_no_list)
				tfbs_row_darwin_ls.append([bs_no_symbol_list, target_gene_no_symbol_list, score])
			#translate the core_vertex_ls
			core_vertex_symbol_ls = dict_map(gene_id2symbol, core_vertex_ls)
			
			#output them all
			outf.write('[%s, %s,{%s},{%s}],\n'%(id, repr(core_vertex_symbol_ls), repr(tfbs_row_darwin_ls)[1:-1], repr(dataset_no_desc_ls)[1:-1]))
	
	def run(self):
		"""
		09-05-05
		2006-09-21 add fuzzyDense_flag
		
			--db_connect()
			--get_gene_no2bs_no_block()
			--construct_two_dicts()
			
			--input_node()
				--fetch_cluster_block()
			--computing_node()
				--node_fire()
					--cluster_bs_analysis()
			--create_cluster_bs_table()
			--output_node()
				--submit_cluster_bs_table()
		"""
		communicator = MPI.world.duplicate()
		node_rank = communicator.rank
		free_computing_nodes = range(1,communicator.size-1)
		
		if node_rank == 0:
			(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
			gene_no2bs_no_block = self.get_gene_no2bs_no_block(curs)
			for node in range(1, communicator.size-1):	#send it to the computing_node
				communicator.send(gene_no2bs_no_block, node, 0)
			if self.fuzzyDense_flag:	#2006-09-21 add fuzzyDense_flag
				#12-18-05 get edge2encodedOccurrence
				MpiCrackSplat_instance = MpiCrackSplat()
				edge2encodedOccurrence = {}
				min_sup = 5	#need to expose them
				max_sup = 40
				total_vertex_set = self.return_total_vertex_set(curs, self.good_cluster_table)
				edge2encodedOccurrence, no_of_datasets = self.fill_edge2encodedOccurrence(\
					self.sig_vector_fname, min_sup, max_sup, total_vertex_set)
				edge2encodedOccurrence_pickle = cPickle.dumps(edge2encodedOccurrence, -1)
				for node in free_computing_nodes:	#send it to the computing_node
					communicator.send(edge2encodedOccurrence_pickle, node, 0)
		elif node_rank>0 and node_rank<communicator.size-1:
			data, source, tag, count = communicator.receive(Numeric.Int, 0, 0)
			gene_no2bs_no_set, bs_no2gene_no_set = self.construct_two_dicts(node_rank, data)
			if self.fuzzyDense_flag:	#2006-09-21
				#12-18-05
				data, source, tag = communicator.receiveString(0, 0)
				edge2encodedOccurrence = cPickle.loads(data)
			
		elif node_rank==communicator.size-1:	#establish connection before pursuing
			(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
			
			#12-20-05 for darwin output
			gene_id2symbol = get_gene_id2gene_symbol(curs, self.tax_id)
			dataset_no2desc = get_dataset_no2desc(curs)
			
			
		mpi_synchronize(communicator)
		
		if node_rank == 0:
			curs.execute("DECLARE crs CURSOR FOR select distinct id, vertex_set, recurrence_array\
				from %s "%(self.good_cluster_table))
			input_node(communicator, curs, free_computing_nodes, self.size, self.report)
			curs.execute("close crs")
			
		elif node_rank<=communicator.size-2:	#exclude the last node
			if self.fuzzyDense_flag:	#2006-09-21
				fuzzyDense_instance = fuzzyDense(edge2encodedOccurrence)
			else:
				fuzzyDense_instance = None
			parameter_list = [gene_no2bs_no_set, bs_no2gene_no_set, self.ratio_cutoff, \
				self.top_number, self.p_value_cut_off, fuzzyDense_instance, self.degree_cut_off, self.fuzzyDense_flag]
			computing_node(communicator, parameter_list, self.computing_node_handler, report=self.report)
			
		elif node_rank==communicator.size-1:
			
			#12-20-05 comment out
			if self.new_table:
				self.create_cluster_bs_table(curs, self.cluster_bs_table)
			parameter_list = [curs, self.cluster_bs_table]
			output_node(communicator, free_computing_nodes, parameter_list, self.submit_cluster_bs_table, report=self.report)
			if self.commit:
				curs.execute("end")
			"""
			outf = open(self.output_file, 'w')
			outf.write("out:=[\n")
			parameter_list = [outf, gene_id2symbol, dataset_no2desc]
			output_node(communicator, free_computing_nodes, parameter_list, self.output_cluster_bs_data, report=self.report)
			
			outf.write('[]]:\n')
			outf.close()
			"""


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:g:l:s:a:t:e:p:o:x:v:fncbr", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	good_cluster_table = None
	cluster_bs_table = None
	size = 2000
	ratio_cutoff = 1.0/3
	top_number = 5
	degree_cut_off = 0.3
	p_value_cut_off = 0.001
	output_file = 'fuzzyDense.out'
	tax_id = 4932
	sig_vector_fname = None
	fuzzyDense_flag = 0
	new_table = 0
	commit = 0
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
		elif opt in ("-g",):
			good_cluster_table = arg
		elif opt in ("-l",):
			cluster_bs_table = arg
		elif opt in ("-s",):
			size = int(arg)
		elif opt in ("-a",):
			ratio_cutoff = float(arg)
		elif opt in ("-t",):
			top_number = int(arg)
		elif opt in ("-e",):
			degree_cut_off = float(arg)
		elif opt in ("-p",):
			p_value_cut_off = float(arg)
		elif opt in ("-f",):
			fuzzyDense_flag = 1
		elif opt in ("-n",):
			new_table = 1
		elif opt in ("-o",):
			output_file = arg
		elif opt in ("-x",):
			tax_id = int(arg)
		elif opt in ("-v",):
			sig_vector_fname = arg
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-b",):
			debug = 1
		elif opt in ("-r",):
			report = 1
	if schema and good_cluster_table and sig_vector_fname:
		if cluster_bs_table or output_file:
			instance = MpiClusterBsStat(hostname, dbname, schema, good_cluster_table,\
				cluster_bs_table, size, ratio_cutoff, top_number, degree_cut_off, p_value_cut_off, output_file, \
				tax_id, sig_vector_fname, fuzzyDense_flag, new_table, commit, debug, report)
			instance.run()
	else:
		print __doc__
		sys.exit(2)
