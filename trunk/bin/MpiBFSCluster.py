#!/usr/bin/env mpipython
"""
Usage: MpiBFSCluster.py -i -o [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database(IGNORE)
	-i ...,	inputfile
	-o ...,	outputfile
	-s ...,	no of clusters per transmission, 1000000(default)
	-n,	output_table is new(IGNORE)
	-c,	commit the database transaction(IGNORE)
	-b,	debug version.
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	mpirun -np 20 -machinefile ~/hostfile /usr/bin/mpipython ~/script/annot/bin/MpiBFSCluster.py
	-i -o -r
	
Description:
	Program to do BFS search for each vertex and record down the layer.
	It's actually all to all shortest path. But BGL's python binding don't offer
	all_pairs shortest path, so use BFS.
	
"""

import sys, os, getopt, csv, Numeric, cPickle
sys.path += [os.path.expanduser('~/script/annot/bin')]
from Scientific import MPI
from codense.common import mpi_synchronize, db_connect, output_node, \
	form_schema_tables, input_node, computing_node
from sets import Set
from graph.johnson_sp import johnson_sp

"""
class cal_distance_bfs_visitor(bgl.Graph.BFSVisitor):
	'''
	10-09-05
		code call this function has been deleted from bfs_cluster()
		see 1.1 version
	'''
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
"""

class MpiBFSCluster:
	"""
	09-19-05
	"""
	def __init__(self,hostname='zhoudb', dbname='graphdb', schema=None, inputfile=None,\
		outputfile=None, size=2000, new_table=0, commit=0, debug=0, report=0):
		"""
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.inputfile = inputfile
		self.outputfile = outputfile
		self.size = int(size)
		self.new_table = int(new_table)
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
	
	def input_handler(self, parameter_list, message_size, report=0):
		"""
		10-28-05 the input is from MpiFromDatasetSignatureToPattern.py
		"""
		if report:
			sys.stderr.write("Fetching stuff...\n")
		inf = parameter_list[0]
		block = []
		string_length = 0
		for row in inf:
			block.append(row)
			string_length += len(repr(row))	#the length to control MPI message size
			if string_length>=message_size:
				break
		if report:
			sys.stderr.write("Fetching done.\n")
		return block
	
	def node_fire(self, communicator, data, parameter_list):
		"""
		10-07-05
		10-21-05
			called by common.computing_node()
		10-28-05
			the input is from MpiFromDatasetSignatureToPattern.py
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		result = []
		for row in data:
			try:
				vertex_set, edge_set, recurrence_array = row
				vertex_set = vertex_set[1:-1].split(',')
				vertex_set = map(int, vertex_set)
				vertex_set.sort()	#10-28-05
				edge_set = edge_set[2:-2].split('], [')	#10-28-05
				for i in range(len(edge_set)):
					edge_set[i] = edge_set[i].split(',')
					edge_set[i] = map(int, edge_set[i])
					edge_set[i].sort()	#10-28-05
				j_instance = johnson_sp()
				j_instance.run(vertex_set, edge_set)
				j_instance.D
				edge_set.sort()	#10-28-05 to ease codense2db.py
				row = [repr(vertex_set), repr(edge_set), row[2], repr(j_instance.D)]	#10-28-05
				result.append(row)
			except:
				common_string = "Node %s error"%communicator.rank
				print common_string, sys.exc_info()
				print common_string, "vertex_set",vertex_set,"edge_set",edge_set
				print common_string, "d_matrix",j_instance.D
		sys.stderr.write("Node no.%s done.\n"%node_rank)
		return result
	
	def cleanup_handler(self, communicator):
		#tell the last node to stop
		communicator.send("-1", communicator.size-1, 1)	#tag is 1
	
	def create_output_table(self, curs, output_table):
		"""
		10-10-05
			splat_id becomes the primary key
		10-20-05
			output_table is just the pattern_table
		10-28-05 DEFUNCT
		"""
		sys.stderr.write("Creating %s...\n"%output_table)
		curs.execute("create table %s(\
			id	serial primary key,\
			vertex_set	integer[],\
			edge_set	integer[][],\
			no_of_vertices	integer,\
			no_of_edges	integer,\
			repr_p_value	float[],\
			repr_go_no	integer[],\
			connectivity	float,\
			unknown_gene_ratio	float,\
			recurrence_array	float[],\
			recurrence	float,\
			d_matrix	integer[][],\
			comment	varchar)"%output_table)
		sys.stderr.write("%s created.\n"%output_table)
	
	def submit2output_table(self, communicator, parameter_list, data):
		"""
		10-07-05
		10-21-05
			called by common.output_node_int()
		10-28-05 DEFUNCT
		"""
		curs, output_table = parameter_list
		data = cPickle.loads(data)
		for row in data:
			id, vertex_set, edge_set, no_of_vertices, no_of_edges, connectivity, \
			unknown_gene_ratio, recurrence_array, recurrence, d_matrix = row
			curs.execute("insert into %s(id, vertex_set, edge_set, no_of_vertices, no_of_edges, \
				connectivity, unknown_gene_ratio, recurrence_array, recurrence, d_matrix) \
				values(%s, '%s', '%s', %s, %s, %s, %s, '%s', %s, ARRAY%s)"%\
				(output_table, id, vertex_set, edge_set, no_of_vertices, no_of_edges, \
				connectivity, unknown_gene_ratio, recurrence_array, recurrence, repr(d_matrix)))
	
	def output_handler(self, communicator, parameter_list, data):
		"""
		10-28-05 write it to outputfile
		"""
		writer = parameter_list[0]
		data = cPickle.loads(data)
		for row in data:
			writer.writerow(row)
	
	def run(self):
		"""
		10-07-05
		10-09-05 input_node() add mcl_table
		10-24-05 create new views for splat_table and mcl_table
		10-28-05 no views, no new pattern_table, read from inputfile, write to outputfile
			
			--input_node()
				--input_handler()
			--computing_node()
				--node_fire()
				--cleanup_handler()
			--output_node()
				--output_handler()
		"""
		communicator = MPI.world.duplicate()
		node_rank = communicator.rank		
		free_computing_nodes = range(1,communicator.size-1)
		if node_rank == 0:
			inf = csv.reader(open(self.inputfile,'r'), delimiter='\t')
			parameter_list = [inf]
			input_node(communicator, parameter_list, free_computing_nodes, self.size, self.report, input_handler=self.input_handler)
			del inf
		elif node_rank<=communicator.size-2:	#exclude the last node
			parameter_list = []
			computing_node(communicator, parameter_list, self.node_fire, self.cleanup_handler, self.report)
		elif node_rank==communicator.size-1:
			writer = csv.writer(open(self.outputfile, 'w'), delimiter='\t')
			parameter_list = [writer]
			output_node(communicator, free_computing_nodes, parameter_list, self.output_handler, self.report)
			del writer

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
	inputfile = None
	outputfile = None
	size = 1000000
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
			inputfile = arg
		elif opt in ("-o"):
			outputfile = arg
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
	if inputfile and outputfile:
		instance = MpiBFSCluster(hostname, dbname, schema, inputfile, outputfile,\
			size, new_table, commit, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
