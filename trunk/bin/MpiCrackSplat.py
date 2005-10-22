#!/usr/bin/env mpipython
"""
Usage: MpiCrackSplat.py -k SCHEMA -i INPUTFILE -o OUTPUTFILE [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-i ...,	--inputfile=...	the input file, dataset signature output by fim
	-o ..., --outputfile=...	the output file
	-m ..., --min_sup=...	minimum support of an edge in the database, 0(default)
	-x ..., --max_sup=...	maximum support of an edge, 200(default)
	-s ..., --min_size=...	minimum size of a cluster(5 default)
	-c ..., --min_con=...	minimum connectivity(0.5, default)
	-y ..., --cluster_block_size=...	min no of clusters in cluster_block_matrix(1500, default)IGNORE
	-w ..., --cluster_block_edges=...	min no of edges in cluster_block_matrix(25,000, default)
	-b, --debug	debug version.
	-r, --report	enable report flag
	-h, --help	Display the usage infomation.
	
Examples:
	mpirun -np 20 -machinefile ~/hostfile /usr/bin/mpipython ~/script/annot/bin/MpiCrackSplat.py
		-k rn_fim_33 -i ~/bin/hhu_clustering/data/output/netmine/rn_fim_33m3x100
		-o ~/bin/hhu_clustering/data/output/netmine/rn_fim_33m3x100.dense -m 3
	
Description:
	Counterpart of CrackSplat.py, but uses ClusterByEBC and MPI.
	Watch -w, control the transmission block size.
	
"""

import sys, os, getopt, csv, math, Numeric
sys.path += [os.path.expanduser('~/script/annot/bin')]
from Scientific import MPI
from graph.cc_from_edge_list import ClusterByEBC
from codense.common import system_call, mpi_schedule_jobs, mpi_synchronize, db_connect, output_node
from netmine_wrapper import netmine_wrapper
from codense.codense2db import codense2db
from sets import Set
from MpiFromDatasetSignatureToPattern import decodeOccurrenceToBv, encodeOccurrenceBv, MpiFromDatasetSignatureToPattern

class MpiCrackSplat:
	def __init__(self,hostname='zhoudb', dbname='graphdb', schema=None, inputfile=None,\
		outputfile=None, min_sup=0, max_sup=200,  min_size=5, min_con=0.5, \
		cluster_block_size=1500, cluster_block_edges=25000, debug=0, report=0):
		"""
		08-07-05
			Program to handle the fim_closed output over the edge-transaction mining
			Restore the patterns(vertex_set, edge_set) based on the dataset signature.
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.inputfile = inputfile
		self.outputfile = outputfile
		self.min_sup = int(min_sup)
		self.max_sup = int(max_sup)
		self.min_size = int(min_size)
		self.min_con = float(min_con)
		self.cluster_block_size = int(cluster_block_size)
		self.cluster_block_edges = int(cluster_block_edges)
		self.debug = int(debug)
		self.report = int(report)
	
	def fill_edge2encodedOccurrence(self, hostname, dbname, schema, edge2encodedOccurrence, min_sup, max_sup, edge_table='edge_cor_vector'):
		"""
		09-05-05
			get the edge2encodedOccurrence from the database
		"""
		sys.stderr.write("Getting edges...\n")
		(conn, curs) = db_connect(hostname, dbname, schema)
		curs.execute("DECLARE crs CURSOR FOR select edge_name,sig_vector \
			from %s"%(edge_table))
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		no_of_datasets = 0
		counter = 0
		while rows:
			for row in rows:
				edge = row[0][1:-1].split(',')
				edge = map(int, edge)
				sig_vector = row[1][1:-1].split(',')
				sig_vector = map(int, sig_vector)
				if no_of_datasets==0:
					no_of_datasets = len(sig_vector)
				if sum(sig_vector)>=min_sup and sum(sig_vector)<=max_sup:
					edge2encodedOccurrence[tuple(edge)] = encodeOccurrenceBv(sig_vector)
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		sys.stderr.write("Done.\n")
		return no_of_datasets
		
	def get_cluster_block(self, reader, min_size, max_no_of_clusters=1000, edge_limit=20000):
		"""
		09-05-05
		
			Format of cluster_block_matrix:
				[vertex1, vertex2, cluster1, cluster2, ...]
			cluster1, cluster2, ... are occurrence flags
			Dictionay edge2index keeps track of the pointer of each edge in the cluster_block_matrix.
			
			Numeric resize is very expensive, try to avoid it.
		09-07-05
			re-encode the cluster_block_matrix
			Format:
				[
				[v1_of_e1, v2_of_e1]
				[v2_of_e1, v2_of_e2]
				...
				[-2,-2]
				[v1_of_e1, v2_of_e1]
				[v2_of_e1, v2_of_e2]
				...
				[-2,-2]
				...
				]
		"""
		cluster_block_matrix = []	#Numeric.zeros((edge_limit,max_no_of_clusters+2), Numeric.Int)
		numb_edge = [-2, -2]	#separator of clusters
		counter = 0
		for row in reader:
			vertex_list = row[0][1:-1].split(',')
			if len(vertex_list)<min_size:	#skip too small clusters
				continue
			edge_list = row[1][2:-2].split('], [')
			for edge in edge_list:
				edge = edge.split(',')
				edge = map(int, edge)
				edge.sort()	#ascending order
				cluster_block_matrix.append(edge)
			cluster_block_matrix.append(numb_edge)	#attach the separator
			counter += 1	#posterior plus
			if len(cluster_block_matrix)>=edge_limit:	#enough, revisit later
				if self.debug:
					sys.stderr.write("No of clustes is %s\n"%counter)
					sys.stderr.write("No of edges is %s.\n"%len(cluster_block_matrix))
				break
		#transform it into Numeric.array
		cluster_block_matrix = Numeric.array(cluster_block_matrix)
		if self.debug:
			sys.stderr.write("The shape of cluster_block_matrix is %s.\n"%(repr(cluster_block_matrix.shape)))
		return cluster_block_matrix
	
	def input_node(self, communicator, inputfile, min_size, max_no_of_clusters, cluster_block_edges):
		"""
		09-05-05
			strategy:
				node 0 is for reading:	read max_no_of_clusters clusters as a block and send it to the sib-nodes
				node size-1 is for writing
				others for computing
				if all computing nodes are too busy to receive the block, the reading will pause
				
			--get_cluster_block()
		10-19-05
			ask output_node for a free_computing_node
		"""
		node_rank = communicator.rank
		sys.stderr.write("Reading clusters from %s...\n"%inputfile)
		reader = csv.reader(open(inputfile, 'r'), delimiter='\t')
		cluster_block_matrix = self.get_cluster_block(reader, min_size, max_no_of_clusters, cluster_block_edges)
		counter = 0	#10-19-05
		while cluster_block_matrix:
			communicator.send(Numeric.array([-2]), communicator.size-1, 1)	#10-19-05 WATCH: tag is 1, to the output_node. Message is array,Float.
			free_computing_node, source, tag = communicator.receiveString(communicator.size-1, 2)	#10-19-05
				#WATCH: tag is 2, from the output_node
			communicator.send(cluster_block_matrix, int(free_computing_node), 0)	#10-19-05	#WATCH: int()
			if self.report:
				sys.stderr.write("block %s sent to %s.\n"%(counter, free_computing_node))	#10-19-05
			cluster_block_matrix = self.get_cluster_block(reader,min_size, max_no_of_clusters, cluster_block_edges)
			counter += 1	#10-19-05
		del reader
		#tell computing_node to exit the loop
		stop_signal = Numeric.zeros((1,1), Numeric.Int)
		stop_signal[0,0] = -1
		for node in range(1, communicator.size-1):	#send it to the computing_node
			communicator.send(stop_signal, node, 0)
		sys.stderr.write("%s reading done\n"%inputfile)
	
	def computing_node(self, communicator, max_no_of_clusters, min_size, min_con):
		"""
		09-05-05
			what the computing_node does
		09-07-05
			the cluster_block_matrix's format was changed in get_cluster_block(), follow it.
			--node_fire()
		"""
		node_rank = communicator.rank
		data, source, tag, count = communicator.receive(Numeric.Int, 0, 0)	#get data from node 0
		while data:
			if data[0]==-1:
				if self.debug:
					sys.stderr.write("node %s breaked.\n"%node_rank)
				break
			else:
				shape_x = count/2
				data.shape = (shape_x, 2)
				self.node_fire(communicator, data, min_size, min_con)
			data, source, tag, count = communicator.receive(Numeric.Int, 0, 0)	#get data from node 0
		#tell the last node to stop
		stop_signal = Numeric.zeros((1,1), Numeric.Int)
		stop_signal[0,0] = -1
		communicator.send(stop_signal, communicator.size-1, 1)
		
	def node_fire(self, communicator, cluster_block_matrix, min_size, min_con):
		"""
		09-05-05
			
			call ClusterByEBC
			pass the results to the last node for output
			format of the passing block
				[
				[shape_x, shape_y,0,0,...]
				[vertex1, vertex2, vertex3, ...]
				[v1_of_e1, v2_of_e1,0,0,...]
				[v1_of_e2, v2_of_e2,0,0,...]
				]
		09-07-05
			the cluster_block_matrix's format was changed in get_cluster_block(), follow it.
			Also change the format of edge_block:
				[
				[no_of_vertices, no_of_edges]
				[v1,0]
				[v2,0]
				...
				[vn,0]
				[v1_of_e1, v2_of_e1]
				...
				]
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		edge_list = []
		for edge in cluster_block_matrix:
			if edge==[-2,-2]:	#cluster separator
				if self.debug:
					sys.stderr.write("The raw edge_list is %s.\n"%repr(edge_list))
				ClusterByEBC_instance = ClusterByEBC(edge_list, min_size, min_con)
				ClusterByEBC_instance.run()
				for k in range(len(ClusterByEBC_instance.cc_list)):
					edge_list = ClusterByEBC_instance.cc_list[k]
					vertex_list = ClusterByEBC_instance.cc_vertex_list[k]
					if self.debug:
						sys.stderr.write("The edge_list is %s.(from node %s)\n"%(repr(edge_list), node_rank))
						sys.stderr.write("The vertex_list is %s.(from node %s)\n"%(repr(vertex_list), node_rank))
					no_of_vertices = len(vertex_list)
					no_of_edges = len(edge_list)
					edge_block = Numeric.zeros((no_of_vertices+no_of_edges+1, 2), Numeric.Int)	#Watch: +2
					edge_block[0,0], edge_block[0,1] = no_of_vertices, no_of_edges
					vertex_list.sort()
					for i in range(no_of_vertices):
						edge_block[i+1][0] = vertex_list[i]
					for l in range(no_of_edges):
						edge = edge_list[l]
						if edge[0]>edge[1]:
							edge.reverse()
						edge_block[no_of_vertices+l+1,0], edge_block[no_of_vertices+l+1,1] = edge
					if self.debug:
						sys.stderr.write("The edge_block is %s.(from node %s)\n"%(repr(edge_block), node_rank))
					communicator.send(edge_block, communicator.size-1, 1)
				del ClusterByEBC_instance
				edge_list= []	#clean up
			else:
				edge_list.append(list(edge))	#Watch: convert to list
		sys.stderr.write("Node no.%s done.\n"%node_rank)
	
	def output_cluster(self, communicator, parameter_list, edge_block):
		"""
		09-05-05
			output the cluster
		09-07-05
			node_fire() changed the format of edge_block, follow it
		10-21-05
			called by common.output_node()
		"""
		writer, codense2db_instance, edge2encodedOccurrence, no_of_datasets = parameter_list
		shape_x = len(edge_block)/2
		edge_block.shape = (shape_x, 2)
		no_of_vertices, no_of_edges = edge_block[0]
		vertex_list = []
		for i in range(1, no_of_vertices+1):
			vertex_list.append(edge_block[i,0])
		
		edge_list = []
		combined_vector = []
		for i in range(no_of_vertices+1, edge_block.shape[0]):
			edge = [edge_block[i,0], edge_block[i,1]]
			edge_list.append(edge)
			edge_tuple = tuple(edge)
			if edge_tuple not in edge2encodedOccurrence:
				sys.stderr.write("edge %s not in edge2encodedOccurrence.\n"%(repr(edge)))
			combined_vector.append(decodeOccurrenceToBv(edge2encodedOccurrence[edge_tuple], no_of_datasets))
		recurrence_array = codense2db_instance.parse_recurrence(combined_vector)
		edge_list.sort()
		writer.writerow([repr(vertex_list), repr(edge_list), repr(recurrence_array)] )
	
	def run(self):
		"""
		09-05-05
			Watch: when sending via MPI, tag 0 means from node 0,  tag 1 means goes to the last node.
		10-21-05
			replace output_node() with the one from codense.common for better scheduling
			
			--fill_edge2encodedOccurrence()
			
			--input_node()
				--get_cluster_block()
			--computing_node()
				--node_fire()
			--output_node()
				--output_cluster()
			
			--uniqueSort()
		"""
		communicator = MPI.world.duplicate()
		node_rank = communicator.rank
		intermediateFile = '%s.unsorted'%self.outputfile	#intermediateFile to store concatenated results
		if communicator.rank==(communicator.size-1):
			edge2encodedOccurrence = {}
			no_of_datasets = self.fill_edge2encodedOccurrence(self.hostname, self.dbname, self.schema, edge2encodedOccurrence, self.min_sup, self.max_sup)
		
		mpi_synchronize(communicator)
		
		if node_rank == 0:
			self.input_node(communicator, self.inputfile, self.min_size, self.cluster_block_size, self.cluster_block_edges)
		elif node_rank<=communicator.size-2:	#exclude the last node
			self.computing_node(communicator, self.cluster_block_size, self.min_size, self.min_con)
		elif node_rank==communicator.size-1:
			codense2db_instance = codense2db()
			free_computing_nodes = range(1,communicator.size-1)
			writer = csv.writer(open(intermediateFile,'w'), delimiter='\t')
			parameter_list = [writer, codense2db_instance, edge2encodedOccurrence, no_of_datasets]
			output_node(communicator, free_computing_nodes, parameter_list, self.output_cluster, report=self.report, type=Numeric.Int)
			del writer
			#10-21-05self.output_node(communicator, intermediateFile, codense2db_instance, edge2encodedOccurrence, no_of_datasets)
		
		mpi_synchronize(communicator)
		#collecting
		if node_rank==0:
			MpiFromDatasetSignatureToPattern_instance = MpiFromDatasetSignatureToPattern()
			MpiFromDatasetSignatureToPattern_instance.uniqueSort(intermediateFile, self.outputfile)
		

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:o:m:x:s:c:y:w:br", ["help", "hostname=", \
			"dbname=", "schema=", "inputfile=", "outputfile=", "min_sup", "max_sup=", "min_size=", \
			"min_con=", "cluster_block_size=", "cluster_block_edges=", "debug", "report"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	inputfile = None
	outputfile = None
	min_sup = 0
	max_sup = 200
	min_size = 5
	min_con = 0.5
	cluster_block_size = 1500
	cluster_block_edges = 25000
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
		elif opt in ("-i", "--inputfile"):
			inputfile = arg
		elif opt in ("-o", "--outputfile"):
			outputfile = arg
		elif opt in ("-m", "--min_sup"):
			min_sup = int(arg)
		elif opt in ("-x", "--max_sup"):
			max_sup = int(arg)
		elif opt in ("-s", "--min_size"):
			min_size = int(arg)
		elif opt in ("-c", "--min_con"):
			min_con = float(arg)
		elif opt in ("-y", "--cluster_block_size"):
			cluster_block_size = int(arg)
		elif opt in ("-w", "--cluster_block_edges"):
			cluster_block_edges = int(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	if schema and inputfile and outputfile:
		instance = MpiCrackSplat(hostname, dbname, schema, inputfile, outputfile, \
			min_sup, max_sup, min_size, min_con, cluster_block_size, cluster_block_edges, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
