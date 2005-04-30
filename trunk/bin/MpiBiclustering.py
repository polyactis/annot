#!/usr/bin/env mpipython
"""
Usage: MpiBiclustering.py -k SCHEMA [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	the splat_result table (edge_set)
	-m ..., --mcl_table=...	the mcl_result table (vertex_set)
	-s ..., --hscore_cut_off=...	0.1 (default)
	-e ..., --min_height=...	5 (default)
	-w ..., --min_width=...		6 (default)
	-o ..., --batchThreshold=...	100 (default)
	-x ..., --cor_cut_off=...	0.9(default), for seed_grow
	-y ..., --euc_dist_cut_off=...	0.02(default) for seed_grow
	-p ..., --outfname=...	the output filename
	-c, --commit	commit this database transaction
	-r, --report	report the progress(a number)
	-b, --debug
	-h, --help	show this help
	
Examples:
	mpirun N MpiBiclustering.py -k sc_54_6661 -r
	
Description:
	Biclustering on the edge matrix. MPI version.

"""

import sys, os, psycopg, getopt, csv, re, random
from sets import Set
from codense.common import db_connect, get_vertex_edge_list_by_edge_id
from Numeric import array, Float, take, zeros, Int
from Scientific import MPI
from graph import biclustering
from graph import graph_modeling

class cluster_group:
	"""
	04-20-05
		class to hold edge cor matrix, edge_id's and all seed-clusters of one function category
	"""
	def __init__(self):
		self.edge_matrix = None
		self.edge_id_array = None
		self.edge_id_set = Set()
		self.bicluster_list = []

class bicluster:
	"""
	04-20-05
		class to hold information of a bicluster, put in cluster_group.bicluster_list
	"""
	def __init__(self):
		self.row_index_list = []
		self.column_index_list = []
		self.score = None
		self.consensus_list = []
		self.added_edge_id_list = []
		self.added_edge_matrix = []

class MpiBiclustering:
	"""
	04-11-05
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, \
		table=None, mcl_table=None, hscore_cut_off=0.1, min_height=5, min_width=6,\
		batchThreshold=100,	cor_cut_off=0.9, euc_dist_cut_off=0.02, \
		outfname=None, report=0, needcommit=0, debug=0):
		
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.table = table
		self.mcl_table = mcl_table
		self.hscore_cut_off = float(hscore_cut_off)
		self.min_height = int(min_height)
		self.min_width = int(min_width)
		self.batchThreshold = int(batchThreshold)
		self.cor_cut_off = float(cor_cut_off)
		self.euc_dist_cut_off = float(euc_dist_cut_off)
		self.outfname = outfname
		self.report = int(report)
		self.needcommit = int(needcommit)
		#debugging flag
		self.debug = int(debug)
		
		#key is go-no, value is [[edge-data],[],[],...] edge-data: edge_id,cor1,cor2...
		self.go_no2edge_matrix_data = {}	#storing the data to generate seed-clusters
		self.go_no2cluster_group = {}	#storing the clusters for each go-no.
		self.candidate_edge_array = []	#edges on which the seed-clusters will grow, later converted to Numeric array to share
		self.no_of_records = 0
		
		#tell boost::python to use Numeric array
		biclustering.set_module_and_type('Numeric', 'ArrayType')

	
	def prepare_gene_no2go_no(self, curs):
		"""
		04-15-05
			different from get_gene_no2go_no, the value is a set.
		04-27-05
			only depth ==5
		"""
		sys.stderr.write("Preparing gene_no2go_no...")
		from codense.common import get_gene_no2go_no, get_go_no2depth
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
	
	def get_function_edge_matrix_data(self, curs, edge_table='edge_cor_vector'):
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
				self._get_function_edge_matrix_data(row)	
				counter +=1
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, counter))
			if self.debug and counter>10000:
				break
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		sys.stderr.write("Done\n")		

	def _get_function_edge_matrix_data(self, row):
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
			self.go_no2edge_matrix_data[go_no].append([edge_id]+self.return_edge_vector(row[2]))
	
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
		edge_vector = []
		for item in edge_vector_string[1:-1].split(','):
			if item=='1.1':	#1.1 is 'NA', convention because of haiyan's copath
				edge_vector.append(random.uniform(-1,1))	#don't use (-800,800)
			else:
				edge_vector.append(float(item))
		return edge_vector
		
	def get_candidate_edge_matrix(self, curs, edge_table='edge_cor_vector'):
		"""
		04-15-05
			fetch all edges on which seed-clusters will grow, (named candidate edges)
			
			NOTE: name another cursor, crs1
			--_get_candidate_edge_matrix_data()
		"""
		sys.stderr.write("Getting candidate edge matrix...\n")
		curs.execute("DECLARE crs1 CURSOR FOR select edge_id,edge_name,cor_vector \
			from %s"%(edge_table))
		
		curs.execute("fetch 5000 from crs1")
		rows = curs.fetchall()
		no_of_records = 0
		while rows:
			for row in rows:
				self._get_candidate_edge_matrix_data(row)
				no_of_records+=1
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, no_of_records))
			if self.debug and no_of_records>10000:
				break
			curs.execute("fetch 5000 from crs1")
			rows = curs.fetchall()
		sys.stderr.write("Done\n")

	def _get_candidate_edge_matrix_data(self, row):
		"""
		04-20-05
			--edge_contain_unknown_gene()
		"""
		edge_id = row[0]
		edge = row[1][1:-1].split(',')
		edge = map(int, edge)
		if self.edge_contain_unknown_gene(edge):
			self.candidate_edge_array.append([edge_id]+self.return_edge_vector(row[2]))
	
	def edge_contain_unknown_gene(self, edge):
		"""
		04-15-05
			return True if the edge contains one unknown gene, False otherwise
		04-27-05
			unknown gene's go_no set is nothing, change the judging rule
		"""
		gene_no1, gene_no2 = edge
		if len(self.gene_no2go_no[gene_no1]) ==0 or len(self.gene_no2go_no[gene_no2])==0:
			return True
		else:
			return False
		
	def seed_grow(self, node_rank, cor_cut_off, euc_dist_cut_off):
		"""
		04-20-05
			add candidate edge based on its correlation with consensus_list, (>=0.8)
		"""
		sys.stderr.write("Node %s, seed_growing...\n"%node_rank)
		for i in range(self.candidate_edge_array.shape[0]):
			candidate_edge_vector = self.candidate_edge_array[i,:]
			edge_id = int(candidate_edge_vector[0])	#first grab the edge_id
			candidate_edge_vector = candidate_edge_vector[1:]	#then grab it's correlation vector
			for go_no, cluster_group in self.go_no2cluster_group.iteritems():
				if edge_id in cluster_group.edge_id_set:
					continue	#this edge is among the function group
				for j in range(len(cluster_group.bicluster_list)):
					bicluster = cluster_group.bicluster_list[j]
					selected_candidate_edge_vector = list(take(candidate_edge_vector, bicluster.column_index_list))
					edge_data = graph_modeling.ind_cor(selected_candidate_edge_vector, \
						bicluster.consensus_list, -1)	#leave_one_out = -1, means no leave_one_out
					euc_edge_data = graph_modeling.euc_dist(selected_candidate_edge_vector,\
						bicluster.consensus_list)
					if edge_data.value>=cor_cut_off and euc_edge_data.value/(euc_edge_data.degree+2)<=euc_dist_cut_off:	#average euclidean distance
						bicluster.added_edge_id_list.append(edge_id)
						bicluster.added_edge_matrix.append(selected_candidate_edge_vector)
						cluster_group.bicluster_list[j] = bicluster	#update the list values, different from dictionary
		sys.stderr.write("Node %s, Done.\n"%(node_rank))
	
	def generate_seeds(self, node_rank, go_no, edge_id_array, edge_matrix):
		"""
		04-20-05
			stop when the hscore_cut_off is not met
		04-29-05
			try 4 times more after hscore_cut_off is not met.
		"""
		sys.stderr.write("Node %s, generate seeds on function %s...\n"%(node_rank, go_no))
		if edge_matrix.shape[0] < self.min_height:
			sys.stderr.write("node %s, go_no %s, no enough rows\n"%(node_rank, go_no))
			return
		self.go_no2cluster_group[go_no] = cluster_group()
		self.go_no2cluster_group[go_no].matrix = edge_matrix
		self.go_no2cluster_group[go_no].edge_id_array = edge_id_array
		self.go_no2cluster_group[go_no].edge_id_set = Set(edge_id_array)
		biclustering_instance = biclustering.biclustering(self.hscore_cut_off, self.min_height, \
			self.min_width, self.batchThreshold)
		biclustering_instance.data_read_in(edge_matrix)
		result = biclustering_instance.getbicluster()
		no_of_retries = 0	#try to get more clusters.
		while result:
			if result.score<=self.hscore_cut_off:
				bicluster_unit = bicluster()
				bicluster_unit.row_index_list = result.row_index_list
				bicluster_unit.column_index_list = result.column_index_list
				bicluster_unit.score = result.score
				bicluster_unit.consensus_list = result.consensus_list
				self.go_no2cluster_group[go_no].bicluster_list.append(bicluster_unit)
			else:
				no_of_retries += 1
				if no_of_retries == 4:
					break
			result = biclustering_instance.getbicluster()
		del biclustering_instance
		sys.stderr.write('Node %s go_no %s Done.\n'%(node_rank, go_no))

	def output_in_copath_format(self, outfname, node_rank):
		"""
		04-20-05
			output go_no2cluster_group
		04-25-05
			cluster_id redefined
		"""
		(conn, curs) = db_connect(self.hostname, self.dbname, self.schema)
		
		outf = open(outfname, 'a')
		writer = csv.writer(outf, delimiter='\t')
		for go_no, cluster_group in self.go_no2cluster_group.iteritems():
			counter = 0
			for bicluster in cluster_group.bicluster_list:
				seed_edge_id_list = list(take(cluster_group.edge_id_array, bicluster.row_index_list))
				edge_id_list = seed_edge_id_list + bicluster.added_edge_id_list
				vertex_list , edge_list = get_vertex_edge_list_by_edge_id(curs, edge_id_list)
				no_of_nodes = len(vertex_list)
				connectivity = len(edge_list)*2.0/(no_of_nodes*(no_of_nodes-1))
				vertex_string = '{' + ';'.join(vertex_list) + ';}'
				edge_string  = self.edge_string_from_edge_list(edge_list)
				cluster_id = "%s.%s"%(go_no, counter)
				writer.writerow([cluster_id, connectivity, vertex_string, edge_string])
				counter += 1
		del writer
		outf.close()
		
	def output(self, outfname, node_rank):
		"""
		04-26-05
			output the information about the bicluster, easy to check
		04-25-05
			cluster_id redefined
		"""
		outf = open(outfname, 'a')
		writer = csv.writer(outf, delimiter='\t')
		for go_no, cluster_group in self.go_no2cluster_group.iteritems():
			counter = 0
			for bicluster in cluster_group.bicluster_list:
				cluster_id = "%s.%s"%(go_no, counter)
				seed_edge_id_list = list(take(cluster_group.edge_id_array, bicluster.row_index_list))
				edge_id_list = seed_edge_id_list + bicluster.added_edge_id_list
				writer.writerow([cluster_id, bicluster.score, repr(edge_id_list), repr(bicluster.column_index_list)])
				counter += 1
	
	def edge_string_from_edge_list(self, edge_list):
		"""
		04-25-05
		"""
		edge_string  = ''
		for edge in edge_list:
			edge_string += "(%s,%s );"%(edge[0], edge[1])
		edge_string = '{'+edge_string+'}'
		return edge_string
	
	def pop_edge_data_of_one_function(self, go_no2edge_matrix_data):
		"""
		04-16-05
			pop edge_data of one function out of the dictionary
			package both go-no and edge_data in numeric array
		04-29-05
			if self.report, output the edge_data.
		"""
		if len(go_no2edge_matrix_data)>0:
			go_no, edge_data = go_no2edge_matrix_data.popitem()
			edge_data[0] = edge_data[0]*len(edge_data[1])	#go-no in the first row is duplicated to be packaged in a Numeric array
			edge_data[0][1] = len(edge_data[1])	#no_of_columns is encoded in the 2nd entry of the first row.
			if self.report:
				self.edge_data_output(go_no, edge_data)
			return array(edge_data)
		else:
			return None

	
	def edge_data_output(self, go_no, edge_data):
		"""
		04-27-05
			output edge data to investigate separately
		"""
		filename = 'edge_data_%s'%go_no
		file = open(filename, 'w')
		writer = csv.writer(file, delimiter='\t')
		for row in edge_data:
			writer.writerow(row)
		del writer
		file.close()
	
	def mpi_synchronize(self, communicator):
		"""
		04-16-05
			synchronize all processes in the communicator
		"""
		sys.stdout.flush()
		sys.stderr.flush()
		communicator.barrier()
	
	def run(self):
		"""
		04-15-05
			if node_rank == 0
				--db_connect()
				--prepare_gene_no2go_no()
				--get_function_edge_matrix_data()
					--_get_function_edge_matrix_data()
			
			--mpi_synchronize()
			if node_rank == 0
				--pop_edge_data_of_one_function()
			else:
				--generate_seeds()
			
			--mpi_synchronize()
			
			if node_rank == 0
				--get_candidate_edge_matrix()
					--_get_candidate_edge_matrix_data()
			
			(first broadcast candidate_edge_array's shape)
			
			(second broadcast candidate_edge_array)
			
			--mpi_synchronize()
			
			--seed_grow()
			
			--mpi_synchronize()
			
			if node_rank == 0:
				(tell other nodes to output)
			else:
				output()
		"""
		communicator = MPI.world.duplicate()
		node_rank = communicator.rank
		
		if node_rank == 0:
			(conn, curs) = db_connect(self.hostname, self.dbname, self.schema)
			self.gene_no2go_no = self.prepare_gene_no2go_no(curs)
			self.get_function_edge_matrix_data(curs)
		
		self.mpi_synchronize(communicator)
		
		#schedule generate_seeds()
		if node_rank == 0:
			seed_utilized = Set()
			for node in range(1, communicator.size):
				edge_data = self.pop_edge_data_of_one_function(self.go_no2edge_matrix_data)
				if edge_data:
					communicator.send(edge_data, node, 0)
					sys.stderr.write("Node %s schedule job to %s\n"%(node_rank, node))
					seed_utilized.add(node)
				else:
					break
			seed_to_tell_output = seed_utilized.copy()	#seed_utilized will be emptied soon, seed_to_tell_output is used later.
			
			received_value, source, tag = communicator.receiveString(None, None)	#listen
			while received_value:		#??check what the received_value is
				edge_data = self.pop_edge_data_of_one_function(self.go_no2edge_matrix_data)
				if edge_data:
					communicator.send(edge_data, source, 0)	#more jobs
					sys.stderr.write("Node %s get one more job\n"%source)
				else:
					stop_signal = array([0.0])
					communicator.send(stop_signal, source, 0)	#no more jobs, stop that node,
						#array([0.0]) is the stopping signal cause the receiver requires the data type to be Numeric.Float
					sys.stderr.write("node %s stopped.\n"%source)
					seed_utilized.remove(source)
					if len(seed_utilized) == 0:	#all seed used have finished their jobs
						break
				received_value, source, tag = communicator.receiveString(None, None)	#listen
		else:
			received_data, source, tag, count = communicator.receive(Float, 0, None)	#get data from node 0,
				#04-24-05 the array is one-dimension no matter what dimension the original array is
			while received_data:
				if count==1:	#it's array([0.0]), stopping signal, don't use 'received_data==array([0.0])' to judge.
					sys.stderr.write("node %s breaked"%node_rank)
					break
				else:
					#handle the received_data
					go_no = int(received_data[0])	#integer form
					no_of_columns = int(received_data[1])	#no_of_columns
					no_of_rows = count/no_of_columns
					received_data.shape = (no_of_rows, no_of_columns)
					edge_id_array = map(int, received_data[1:,0])	#integer form
					edge_matrix = received_data[1:,1:]
					self.generate_seeds(node_rank, go_no, edge_id_array, edge_matrix)	#generate seeds
					communicator.send("finished", 0, node_rank)
					
				received_data, source, tag, count = communicator.receive(Float, 0, None)	#get data from node 0
		
		self.mpi_synchronize(communicator)
		
		
		candidate_edge_array_shape = zeros((2), Int)
		
		#get the candidate edges, first broadcast its shape
		if node_rank == 0:
			self.get_candidate_edge_matrix(curs)
			self.candidate_edge_array = array(self.candidate_edge_array)	#convert to Numeric array
			candidate_edge_array_shape = array(self.candidate_edge_array.shape)
		
		communicator.broadcast(candidate_edge_array_shape,0)
		
		self.mpi_synchronize(communicator)
		
		#broadcast candidate edges
		if node_rank!=0:
			self.candidate_edge_array = zeros(candidate_edge_array_shape, Float)
		communicator.broadcast(self.candidate_edge_array, 0)
		
		self.mpi_synchronize(communicator)
		
		#seed_grow
		if node_rank!=0:
			self.seed_grow(node_rank, self.cor_cut_off, self.euc_dist_cut_off)
		
		self.mpi_synchronize(communicator)
		
		if node_rank == 0:
			for node in seed_to_tell_output:
				communicator.send("output", node, 0)
				return_value, source, tag = communicator.receiveString(node, None)
				if return_value=="finished":
					print "%s finished its output"%(source)
				else:
					print "%s encounted error: %s in its output"%(source, return_value)
		else:
			return_value, source, tag = communicator.receiveString(0, None)	#get 'output' signal from node 0
			if return_value=='output':
				self.output_in_copath_format(self.outfname,node_rank)
				self.output(self.outfname+'.1', node_rank)
				communicator.send("finished", 0, node_rank)	#tell node 0 after done.

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema=", "table=", "mcl_table=",\
		"hscore_cut_off=", "min_height=", "min_width=", "batchThreshold=", "cor_cut_off=",\
		"euc_dist_cut_off=", "outfname=", "report", "commit", "debug"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:m:s:e:w:o:x:y:p:rcb", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = 'cluster_stat'
	mcl_table = 'mcl_result'
	hscore_cut_off= 0.1
	min_height = 5
	min_width = 6
	batchThreshold = 100
	cor_cut_off = 0.9
	euc_dist_cut_off = 0.02
	outfname = None
	report = 0
	commit = 0
	debug = 0
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
		elif opt in ("-m", "--mcl_table"):
			mcl_table = arg
		elif opt in ("-s", "--hscore_cut_off"):
			hscore_cut_off = float(arg)
		elif opt in ("-e", "--min_height"):
			min_height = int(arg)
		elif opt in ("-w", "--min_width"):
			min_width = int(arg)
		elif opt in ("-o", "--batchThreshold"):
			batchThreshold = int(arg)
		elif opt in ("-x", "--cor_cut_off"):
			cor_cut_off = float(arg)
		elif opt in ("-y", "--euc_dist_cut_off"):
			euc_dist_cut_off = float(arg)
		elif opt in ("-p", "--outfname"):
			outfname = arg
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-b", "--debug"):
			debug = 1

	if schema and outfname:
		instance = MpiBiclustering(hostname, dbname, schema, table, mcl_table, hscore_cut_off,\
			min_height, min_width, batchThreshold, cor_cut_off, euc_dist_cut_off, outfname, \
			report, commit, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
