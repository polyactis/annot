#!/usr/bin/env mpipython
"""
Usage: MpiTightClust.py [OPTIONS] -k SCHEMA --od=OUTPUT_DIR

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	the edge_cor_vector table, edge_cor_vector(default)
	--od=...	output_dir
	--me=...	min_no_of_edges,500(default)
	
	--tp=..	top_percentage, 0.1(default)
	--nn=...	no_of_nas maximum(<=), 7(default)
	
	--tcn=...	targetClustNum, 2(default)
	--mk=...	min_k, 4(default)
	--xk=...	max_k, 8(default)
	-a ...	alpha, 0.1(default)
	-b ...	beta, 0.6(default)
	--tn=...	topNum, 1(default)
	--sn=..	seqNum, 2(default)
	--rsn=...	resampNum, 10(default)
	--ssp=...	subSampPercent, 0.7(default)
	--np=...	npass, 3(default)
	
	--debug	enable debug, test running
	-h	--help	show this help
	
Examples:
	mpirun.mpich -np 20 -machinefile ~/hostfile /usr/bin/mpipython
	~/script/annot/bin/MpiTightClust.py -k mm_79_2 --od=./edge_data/mm_79_2
	
Description:
	The program gets all the edge_data, preprocesses them into column format.
	(by module PreprocessEdgeData) and then call tightClust(from Prof Tseng).
	
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import sys, os, math, getopt, time, csv, Numeric, psycopg
from Scientific import MPI
from codense.common import system_call, mpi_schedule_jobs, mpi_synchronize, db_connect, get_go_no2depth, get_gene_no2go_no
from sets import Set
from PreprocessEdgeData import PreprocessEdgeData
from netmine_wrapper import netmine_wrapper

def callTightClust(go_no, parameter_list):
	"""
	06-03-05
		the node_function for mpi_schedule_jobs()
	"""
	import sys
	output_dir = parameter_list[0]
	no_of_nas = parameter_list[1]
	top_percentage = parameter_list[2]
	input_file = os.path.join(output_dir, 'edge_data_%s'%go_no)
	output_file = os.path.join(output_dir, 'edge_data_%s.2'%go_no)
	preprocess_instance = PreprocessEdgeData(input_file, output_file, no_of_nas, top_percentage)
	preprocess_instance.run()
	tightClust_path = os.path.join(os.path.expanduser('~/script/annot/bin/tightClust'), 'tightClust')
	job_list = [tightClust_path, '%s/'%output_dir, 'edge_data_%s.2'%go_no] + parameter_list[3:]
	exit_code = system_call("%s"%' '.join(job_list))
	
	tightClust_output_file = '%s.tightClust'%output_file
	del preprocess_instance
	return tightClust_output_file

class get_edge_data:
	"""
	06-07-05
		no_of_nas is used to control whether the edge should be outputed.
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, 
		table='edge_cor_vector', output_dir=None, min_no_of_edges=500, debug=0, no_of_nas=7):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.table = table
		self.output_dir = output_dir
		self.min_no_of_edges = int(min_no_of_edges)
		self.debug = int(debug)
		self.no_of_nas = int(no_of_nas)
		
		self.gene_no2go_no = {}
		self.go_no2edge_matrix_data = {}
		self.go_no_qualified = []
	
	def run(self):
		"""
		06-03-05
		
		--db_connect()
		--prepare_gene_no2go_no()
		--get_function_edge_matrix_data()
			--_get_function_edge_matrix_data()
				--return_common_go_no_of_edge()
				--return_edge_vector()
		--edge_data_output()
		"""		
		conn,curs = db_connect(self.hostname, self.dbname,self.schema)
		
		self.gene_no2go_no = self.prepare_gene_no2go_no(curs)
		self.get_function_edge_matrix_data(curs, self.no_of_nas, self.table)
		
		#make a directory first
		if not os.path.isdir(self.output_dir):
			os.makedirs(self.output_dir)
		
		for go_no, edge_data in self.go_no2edge_matrix_data.iteritems():
			if len(edge_data)>=self.min_no_of_edges:
				self.edge_data_output(self.output_dir, go_no, edge_data)
				self.go_no_qualified.append(go_no)
		
		
	def get_function_edge_matrix_data(self, curs, no_of_nas, edge_table='edge_cor_vector'):
		"""
		04-15-05
			
		"""
		sys.stderr.write("Getting edge matrix for all functions...\n")
		curs.execute("DECLARE crs CURSOR FOR select edge_id,edge_name,cor_vector \
			from %s"%(edge_table))
		
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		counter = 0
		while rows:
			for row in rows:
				self._get_function_edge_matrix_data(row, no_of_nas)
				counter +=1
			sys.stderr.write('%s%s'%('\x08'*20, counter))
			if self.debug:
				if counter == 10000:
					break
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		sys.stderr.write("Done\n")


	def _get_function_edge_matrix_data(self, row, no_of_nas):
		"""
		04-11-05
		"""
		edge_id = row[0]
		edge = row[1][1:-1].split(',')
		edge = map(int, edge)
		for go_no in self.return_common_go_no_of_edge(edge):
			if go_no not in self.go_no2edge_matrix_data:
				self.go_no2edge_matrix_data[go_no] = [[go_no]]	#later expanded in pop_edge_data_of_one_function()
					#to be packaged into a Numeric array
			(na_counter, edge_vector) = self.return_edge_vector(row[2])
			if na_counter<=no_of_nas:
				self.go_no2edge_matrix_data[go_no].append([edge_id]+edge_vector)
	
	def return_common_go_no_of_edge(self, edge):
		"""
		04-15-05
			return common go_no's shared by an edge.
		"""
		gene_no1, gene_no2 = edge
		common_go_no_set = self.gene_no2go_no[gene_no1] & self.gene_no2go_no[gene_no2]
		return common_go_no_set
	
	def return_edge_vector(self, edge_vector_string):
		"""
		04-16-05
			parse the edge_vector_string fetched from database,
			handle the NA issues, replace them with random numbers for Cheng2000's biclustering
		"""
		na_counter = 0
		edge_vector = []
		for item in edge_vector_string[1:-1].split(','):
			if item=='1.1':	#1.1 is 'NA', convention because of haiyan's copath
				edge_vector.append('NA')	#don't use (-800,800)
				na_counter += 1
			else:
				edge_vector.append(float(item))
		return (na_counter, edge_vector)
	
	def prepare_gene_no2go_no(self, curs):
		"""
		04-15-05
			different from get_gene_no2go_no, the value is a set.
		04-27-05
			only depth ==5
		"""
		sys.stderr.write("Preparing gene_no2go_no...")
		#from codense.common import get_gene_no2go_no, get_go_no2depth
		go_no2depth = get_go_no2depth(curs)
		gene_no2go_no = get_gene_no2go_no(curs)
		gene_no2go_no_set = {}
		for gene_no,go_no_list in gene_no2go_no.iteritems():
			gene_no2go_no_set[gene_no] = Set()
			for go_no in go_no_list:
				if go_no2depth[go_no] == 5:
					gene_no2go_no_set[gene_no].add(go_no)
		sys.stderr.write("Done.\n")
		return gene_no2go_no_set
	
	def edge_data_output(self, output_dir, go_no, edge_data):
		"""
		04-27-05
			output edge data to investigate separately
		06-03-05
			add the output_dir
		"""
		filename = 'edge_data_%s'%go_no
		filename = os.path.join(output_dir, filename)
		file = open(filename, 'w')
		writer = csv.writer(file, delimiter='\t')
		for row in edge_data:
			writer.writerow(row)
		del writer
		file.close()


class MpiTightClust:
	def __init__(self,  hostname='zhoudb', dbname='graphdb', schema=None, \
		table='edge_cor_vector', output_dir=None, min_no_of_edges=500, top_percentage=0.1, \
		no_of_nas=7, targetClustNum='2', min_k='4', max_k='8', alpha='0.1', \
		beta='0.6', topNum='1', seqNum='2', resampNum='10', subSampPercent='0.7', \
		npass='3', debug=0):
		
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.table = table
		self.output_dir = os.path.expanduser(output_dir)	#if '~' is included in it, the programs called won't recognize it.
		self.min_no_of_edges = min_no_of_edges
		self.top_percentage = top_percentage
		self.no_of_nas = no_of_nas
		self.targetClustNum = targetClustNum
		self.min_k = min_k
		self.max_k = max_k
		self.alpha = alpha
		self.beta = beta
		self.topNum = topNum
		self.seqNum = seqNum
		self.resampNum = resampNum
		self.subSampPercent = subSampPercent
		self.npass = npass
		self.debug = debug

	def run(self):
		"""
		06-03-05
			
			--<get_edge_data>
			
			--mpi_schedule_jobs()
				--callTightClust()
					--<PreprocessEdgeData>
					--tightClust
			--<netmine_wrapper>
			
		"""
		communicator = MPI.world.duplicate()
		get_edge_data_instance = get_edge_data(self.hostname, self.dbname, self.schema,\
				self.table, self.output_dir, self.min_no_of_edges, self.debug, self.no_of_nas)
		
		if communicator.rank == 0:
			sys.stderr.write("this is node %s\n"%communicator.rank)
			get_edge_data_instance.run()
		
		mpi_synchronize(communicator)
		
		job_list = get_edge_data_instance.go_no_qualified
		parameter_list =[self.output_dir, self.no_of_nas, self.top_percentage, self.targetClustNum, \
			self.min_k, self.max_k, self.alpha, self.beta, self.topNum, self.seqNum, self.resampNum,\
			self.subSampPercent, self.npass]
		if self.debug:
			sys.stderr.write("The common parameter_list is %s.\n"%repr(parameter_list))
		of_name_list = mpi_schedule_jobs(communicator, job_list, callTightClust, parameter_list, self.debug)
		
		mpi_synchronize(communicator)
		
		#collecting
		if communicator.rank==0:
			final_ofname = os.path.join(self.output_dir, 'tightClust')
			netmine_wrapper_instance = netmine_wrapper()
			netmine_wrapper_instance.collect_and_merge_output(of_name_list, final_ofname)
				#receive only from those fired nodes

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema=", "table=", "od=",\
		"me=", "tp=", "nn=", "tcn=", "mk=", "xk=", "tn=", "sn=", "rsn=", "ssp=",\
		"np=", "debug"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:a:b:", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = None
	#passed to get_edge_data
	table = 'edge_cor_vector'
	output_dir = None
	min_no_of_edges = 500
	
	#passed to PreprocessEdgeData
	top_percentage = 0.1
	no_of_nas = 7
	
	#passed to tightClust, string form
	targetClustNum = "2"
	min_k = "4"
	max_k = '8'
	alpha = '0.1'
	beta = '0.6'
	topNum = '1'
	seqNum = '2'
	resampNum = '10'
	subSampPercent = '0.7'
	npass = '3'
	
	debug=0
	
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
		elif opt in ("-t", "--table"):
			table = arg
		elif opt in ("--od"):
			output_dir = arg
		elif opt in ("--me"):
			min_no_of_edges = int(arg)
		elif opt in ("--tp"):
			top_percentage = float(arg)
		elif opt in ("--nn"):
			no_of_nas = int(arg)
		elif opt in ("--tcn"):
			targetClustNum = arg
		elif opt in ("--mk"):
			min_k = arg
		elif opt in ("--xk"):
			max_k = arg
		elif opt in ("-a"):
			alpha = arg
		elif opt in ("-b"):
			beta = arg
		elif opt in ("--tn"):
			topNum = arg
		elif opt in ("--sn"):
			seqNum = arg
		elif opt in ("--rsn"):
			resampNum = arg
		elif opt in ("--ssp"):
			subSampPercent = arg
		elif opt in ("--np"):
			npass = arg
		elif opt in ("--debug"):
			debug = 1

	if schema and output_dir:
		instance = MpiTightClust(hostname, dbname, schema, table, output_dir, min_no_of_edges,\
			top_percentage, no_of_nas, targetClustNum, min_k, max_k, alpha, beta,\
			topNum, seqNum, resampNum, subSampPercent, npass, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
