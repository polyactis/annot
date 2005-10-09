#!/usr/bin/env mpipython
"""
Usage: MpiBFSCluster.py -k SCHEMA -i SPLAT_TABLE -o D_MATRIX_TABLE [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-i ...,	SPLAT_TABLE, input_table
	-o ...,	D_MATRIX_TABLE, output_table
	-s ...,	no of clusters per transmission, 2000(default)
	-n,	D_MATRIX_TABLE is new
	-c,	commit the database transaction
	-b,	debug version.
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	mpirun -np 20 -machinefile ~/hostfile /usr/bin/mpipython ~/script/annot/bin/MpiBFSCluster.py
	-k mm_fim_97 -i -o -n -r -c
	
Description:
	Program to do BFS search for each vertex and record down the layer.
	It's actually all to all shortest path. But BGL's python binding don't offer
	all_pairs shortest path, so use BFS.
	
"""

import sys, os, getopt, csv, math, Numeric
sys.path += [os.path.expanduser('~/script/annot/bin'), os.path.expanduser('~/lib/boost/lib')]
from Scientific import MPI
from codense.common import mpi_synchronize, db_connect
from sets import Set
import bgl

class cal_distance_bfs_visitor(bgl.Graph.BFSVisitor):
	def __init__(self, name_map, vertex_no2index):
		bgl.Graph.BFSVisitor.__init__(self)
		self.name_map = name_map
		self.vertex_no2index = vertex_no2index
		self.source = -1	#None is not good to judge, index 0 could be mixed with it.
		self.distance_list = [-1]*len(vertex_no2index)	#-1 means infiniti
	
	def discover_vertex(self, u, g):
		vertex_no = int(self.name_map[u])
		vertex_index = self.vertex_no2index[vertex_no]
		if self.source!=-1:	#not just "if self.source:", the latter confuses with 0 and None
			self.distance_list[vertex_index] = self.distance_list[self.source]+1
		else:
			self.distance_list[vertex_index] = 0
	
	def examine_vertex(self, u, g):
		vertex_no = int(self.name_map[u])
		vertex_index = self.vertex_no2index[vertex_no]
		self.source = vertex_index

class MpiBFSCluster:
	"""
	09-19-05
	"""
	def __init__(self,hostname='zhoudb', dbname='graphdb', schema=None, input_table=None,\
		output_table=None,  size=2000, new_table=0, commit=0, debug=0, report=0):
		"""
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_table = input_table
		self.output_table = output_table
		self.size = int(size)
		self.new_table = int(new_table)
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
	
	def fetch_cluster_block(self, curs, size):
		"""
		10-07-05
			format:
				[
				[splat_id1, 0],
				[g1,g2],
				[g3,g4],
				[-1, -1]
				...
				[splat_id2, 0],
				[g1, g2],
				...
				[-1, -1]
				
				]
			[-1, -1] is cluster separator
		"""
		curs.execute("fetch %s from crs"%(size))
		rows = curs.fetchall()
		cluster_block = []
		for row in rows:
			splat_id, edge_set = row
			cluster_block.append([splat_id, 0])
			edge_set = edge_set[2:-2].split('},{')
			for edge in edge_set:
				edge = edge.split(',')
				edge = map(int, edge)
				cluster_block.append(edge)
			cluster_block.append([-1,-1])
		cluster_block = Numeric.array(cluster_block, Numeric.Int)
		return cluster_block
	
	def input_node(self, communicator, curs, input_table, size):
		sys.stderr.write("Reading clusters from tables...\n")
		node_rank = communicator.rank
		curs.execute("DECLARE crs CURSOR FOR select splat_id, edge_set\
			from %s "%(input_table))
		cluster_block = self.fetch_cluster_block(curs, size)
		which_node = 0
		while cluster_block:
			node_to_receive_block = which_node%(communicator.size-2)+1	#it's -2, work like this. 
				#"which_node%(communicator.size-2)" is in range (0,size-3). Regard node 1 to size-2 as 0 to size-3. So need +1.
			communicator.send(cluster_block, node_to_receive_block, 0)
			if self.debug:
				sys.stderr.write("block %s sent to %s.\n"%(which_node, node_to_receive_block))
			cluster_block = self.fetch_cluster_block(curs, size)
			which_node += 1
		#tell computing_node to exit the loop
		stop_signal = Numeric.zeros((1,1), Numeric.Int)
		stop_signal[0,0] = -1
		for node in range(1, communicator.size-1):	#send it to the computing_node
			communicator.send(stop_signal, node, 0)
		sys.stderr.write("Node no.%s reading done\n"%(node_rank))
	
	def bfs_cluster(self, splat_id_edge_block, communicator):
		"""
		10-07-05
			result_block format:
				[
				[no_of_genes, splat_id, 0, 0,...]
				[g1, g2, ..., gn]
				[g1-g1, g1-g2, g1-g3, ..., g1-gn]
				...
				[gn-g1, gn-g2, ..., gn-gn]
				]
		"""
		edge_set = splat_id_edge_block[1:]
		g = bgl.Graph(edge_set, "nm")
		name_map = g.get_vertex_string_map("nm")	#WATCH: bgl transforms the integers in edge_set into strings, so int() later
		vertex_set = []
		vertex_no2index = {}
		#construct vertex_no2index
		for v in g.vertices:
			vertex_set.append(int(name_map[v]))
		vertex_set.sort()
		for i in range(len(vertex_set)):
			vertex_no2index[vertex_set[i]] = i
		
		#result
		try:
			result_block = Numeric.zeros((len(vertex_set)+2, len(vertex_set)), Numeric.Int)
			result_block[0][0] = len(vertex_set)
			result_block[0][1] = splat_id_edge_block[0][0]
			result_block[1] = vertex_set
			for v in g.vertices:
				bfs_visitor_instance = cal_distance_bfs_visitor(name_map, vertex_no2index)
				bgl.breadth_first_search(g, v, visitor=bfs_visitor_instance)
				vertex_index = vertex_no2index[int(name_map[v])]
				result_block[vertex_index+2] = bfs_visitor_instance.distance_list
		except:
			common_string = "Node %s error"%communicator.rank
			print common_string, sys.exc_info()
			print common_string, "splat_id_edge_block:",splat_id_edge_block[:50]
			print common_string, "result_block:", result_block[:50]
		communicator.send(result_block, communicator.size-1, 1)
	
	def node_fire(self, communicator, data):
		"""
		10-07-05
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		splat_id_edge_block =[]
		for i in range(len(data)/2):
			if data[2*i]==-1:	#another new edge_set
				self.bfs_cluster(splat_id_edge_block, communicator)
				splat_id_edge_block = []	#re-initialize
			else:
				splat_id_edge_block.append([data[2*i], data[2*i+1]])
		sys.stderr.write("Node no.%s done.\n"%node_rank)
	
	def computing_node(self, communicator):
		"""
		10-07-05
			
		"""
		node_rank = communicator.rank
		data, source, tag, count = communicator.receive(Numeric.Int, 0, 0)	#get data from node 0
		while 1:
			if data[0]==-1:
				if self.debug:
					sys.stderr.write("node %s breaked.\n"%node_rank)
				break
			else:
				self.node_fire(communicator, data)
			data, source, tag, count = communicator.receive(Numeric.Int, 0, 0)	#get data from node 0
		#tell the last node to stop
		stop_signal = Numeric.zeros((1,1), Numeric.Int)
		stop_signal[0,0] = -1
		communicator.send(stop_signal, communicator.size-1, 1)
	
	def create_output_table(self, curs, output_table):
		sys.stderr.write("Creating %s...\n"%output_table)
		curs.execute("create table %s(\
			id	serial primary key,\
			splat_id	integer,\
			vertex_set	integer[],\
			d_matrix	integer[][])"%(output_table))
		sys.stderr.write("%s created.\n"%output_table)
	
	def submit2output_table(self, curs, output_table, data):
		"""
		10-07-05
		"""
		no_of_genes = int(data[0])
		splat_id = int(data[1])
		data = Numeric.resize(data, (no_of_genes+2, no_of_genes))
		vertex_set = data[1].tolist()
		d_matrix = data[2:].tolist()
		d_matrix_string = repr(d_matrix)
		d_matrix_string = d_matrix_string.replace('[', '{')
		d_matrix_string = d_matrix_string.replace(']', '}')
		curs.execute("insert into %s(splat_id, vertex_set, d_matrix) values(%s, '{%s}', '%s')"%\
			(output_table, splat_id, repr(vertex_set)[1:-1], d_matrix_string))
			
	
	def output_node(self, communicator, curs, output_table):
		"""
		10-07-05
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s ready to accept output...\n"%node_rank)
		data, source, tag, count = communicator.receive(Numeric.Int, None, 1)
		no_of_resting_nodes = 0	#to keep track how many computing_nodes have rested
		while 1:
			if data[0]==-1.0:
				no_of_resting_nodes += 1
				if self.debug:
					sys.stderr.write("node %s(%s-th) rested.\n"%(source, no_of_resting_nodes))
				if no_of_resting_nodes==communicator.size-2:	#all computing_nodes have stopped, i'm done.
					break
					if self.debug:
						sys.stderr.write("node %s output finished.\n"%node_rank)
			else:
				self.submit2output_table(curs, output_table, data)
			data, source, tag, count = communicator.receive(Numeric.Int, None, 1)
		sys.stderr.write("Node no.%s output done.\n"%node_rank)
	
	def run(self):
		"""
		10-07-05
			
			--db_connect()
			
			--input_node()
				--fetch_cluster_block()
			--computing_node()
				--node_fire()
					--bfs_cluster()
			--create_output_table()
			--output_node()
				--submit2output_table()
		"""
		communicator = MPI.world.duplicate()
		node_rank = communicator.rank
		if node_rank == 0:
			(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		elif node_rank==communicator.size-1:	#establish connection before pursuing
			(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		mpi_synchronize(communicator)
		
		if node_rank == 0:
			self.input_node(communicator, curs, self.input_table, self.size)
		elif node_rank<=communicator.size-2:	#exclude the last node
			self.computing_node(communicator)
		elif node_rank==communicator.size-1:
			if self.new_table:
				self.create_output_table(curs, self.output_table)
			self.output_node(communicator, curs, self.output_table)
			if self.commit:
				curs.execute("end")

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:o:s:ncbr", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	input_table = None
	output_table = None
	size = 2000
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
		elif opt in ("-i"):
			input_table = arg
		elif opt in ("-o"):
			output_table = arg
		elif opt in ("-s"):
			size = int(arg)
		elif opt in ("-n"):
			new_table = 1
		elif opt in ("-c"):
			commit = 1
		elif opt in ("-b"):
			debug = 1
		elif opt in ("-r"):
			report = 1
	if schema and input_table and output_table:
		instance = MpiBFSCluster(hostname, dbname, schema, input_table,\
			output_table, size, new_table, commit, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
