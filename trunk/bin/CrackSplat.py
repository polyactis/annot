#!/usr/bin/env python
"""
Usage: CrackSplat.py -k SCHEMA [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	the splat_result table (edge_set)
	-m ..., --mcl_table=...	the mcl_result table (vertex_set)
	-f ,,,, --dir_files=...	the temporary directory to store temporary files,
		/tmp/CrackSplat(default).
	-a ..., --parameter_a=...	the first parameter(connectivity), 0.2(default)
	-b ..., --parameter_b=...	the second parameter(size), 4(default)
	-y ..., --type=...	the type, 1(modes, default), 2
	-n,	--new_table	target_table is new. (create it first)
	-c, --commit	commit this database transaction
	-r, --report	report flag
	-b, --debug debug flag
	-h, --help              show this help

Examples:
	CrackSplat.py -k sc_54 -t splat_result_p3g5e6d4q5n80 -m mcl_result_p3g5e6d4q5n80
		-r
	
Description:
	02-24-05
	Break copath results, which are stored in splat_result table, where the name
	comes, into dense parts via calling other clustering algorithm.
	The results will be stored in a new mcl_result-like table.
"""

import sys, os, psycopg, getopt, csv
from graphlib import Graph
from sets import Set
from visualize.clustering_test import clustering_test
from codense.codense2db import codense2db

class MclResult:
	"""
	data structure for mcl_result
	"""
	def __init__(self):
		self.splat_id = -1
		self.vertex_set = []
		self.connectivity = 0
		self.recurrence_array = []
		
class CrackSplat:
	"""
	02-24-05
		Break copath results, which are stored in splat_result table, where the name
		comes, into dense parts via calling other clustering algorithm.
		The results will be stored in a new mcl_result-like table.
	02-25-05
		a debug flag	
		--run
	"""
	def __init__(self, hostname=None, dbname=None, schema=None, table=None, mcl_table=None, dir_files='/tmp/CrackSplat',\
		parameter_a=0.2, parameter_b=4, type=1, new_table=0, needcommit=0, report=0, debug=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.table = table
		self.mcl_table = mcl_table
		self.dir_files = dir_files
		self.parameter_a = float(parameter_a)
		self.parameter_b = int(parameter_b)
		self.type = int(type)
		self.new_table = int(new_table)
		self.needcommit = int(needcommit)
		self.report = int(report)
		self.debug = int(debug)
	
	def init(self):
		"""
		02-24-05
			instantiate a class, create the temp directory if necessary,
		"""
		self.clustering_test_instance = clustering_test()
		self.codense2db_instance = codense2db()
		
		if not os.path.isdir(self.dir_files):
			os.makedirs(self.dir_files)
		else:
			sys.stderr.write("Warning, directory %s already exists.\n"%(self.dir_files))
		self.tmpinfname = os.path.join(self.dir_files, 'input')
		self.tmpoutfname = os.path.join(self.dir_files, 'output')
		
		self.crack_dict = {1: self.crack_by_modes}
		
	def db_connect(self, hostname, dbname, schema):
		"""
		02-24-05
			establish database connection, return (conn, curs).
		"""
		conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		curs = conn.cursor()
		curs.execute("set search_path to %s"%schema)
		return (conn, curs)
		
	def splat2graph_dict(self, edge_set):
		"""
		02-24-05
			convert a splat pattern, an edge_list to a graph_dict.
		"""
		graph_dict = {}
		edge_list = edge_set[2:-2].split('},{')
		for edge in edge_list:
			vertex_list = edge.split(',')
			vertex1 = int(vertex_list[0])
			vertex2 = int(vertex_list[1])
			if vertex1 < vertex2:
				graph_dict[(vertex1, vertex2)] = 1
			else:
				graph_dict[(vertex2, vertex1)] = 1
		return graph_dict
	
	def crack_by_modes(self, splat_id, graph_dict, infname, outfname, parameter_a, parameter_b, \
		argument1=None, argument2=None, argument3=None):
		"""
		02-24-05
			--graph_dict2modes_input
			--clustering_test_instance.reformat
			--call_modes
			--parse_modes_results
		02-25-05
			check the exit code of modes
		02-25-05
			debug flag used here to visualize the modes cracking result
		#argument1 and argument2 are for future purpose
		
		"""
		(index2no, graph) = self.graph_dict2graph(graph_dict)
		clustering_test_instance = argument1
		clustering_test_instance.reformat(graph, infname, len(index2no))
		return_code = self.call_modes(infname, outfname, len(index2no), parameter_a, parameter_b)
		if return_code!=1:
			#modes' normal exit code is 1
			print 'call modes failed'
			sys.exit(1)
		if self.debug:
			clustering_test_instance.visualize_clusters(outfname, graph, index2no, '/tmp/test.R')
		codense2db_instance = argument2
		curs = argument3
		return self.parse_modes_results(splat_id, outfname, index2no, graph, codense2db_instance, curs)
		
	def graph_dict2graph(self, graph_dict):
		"""
		02-24-05
			return (index2no, graph)
			the vertices in the graph are reordered from 0 to len(index2no)-1
		"""
		no2index = {}		#used to keep track of whether one node has appeared or not
		index2no = {}
		graph = Graph.Graph()
		no_of_genes = 0
		for (edge, weight) in graph_dict.iteritems():
			if edge[0] not in no2index:
				index1 = no_of_genes
				no2index[edge[0]] = index1				
				index2no[index1] = edge[0]
				no_of_genes += 1
			else:
				index1 = no2index[edge[0]]
				
			if edge[1] not in no2index:
				index2 = no_of_genes
				no2index[edge[1]] = index2
				index2no[index2] = edge[1]
				no_of_genes += 1
			else:
				index2 = no2index[edge[1]]
			if index1<index2:
				graph.add_edge(index1, index2, weight)
			else:
				graph.add_edge(index2, index1, weight)
		
		return (index2no, graph)
	
	def call_modes(self, infname, outfname, no_of_genes, parameter_a=0.2, parameter_b=4):
		"""
		02-24-05
			make the real call to modes.
		
		modes' parameter order
		inputfile gene_num outputfile min_graph_size min_edge_weight densitycutoff_order*10 cutStopSize

		"""
		modes_path = os.path.join(os.path.expanduser('~'), 'bin/hhu_clustering/bin/modes')
		if not os.path.isfile(modes_path):
			print 'modes program not available'
			sys.exit(1)
		###!!!, the argument list must comprise of all strings, not integer, or float.
		wl = ['modes', infname, repr(no_of_genes), outfname, repr(parameter_b), '1', repr(int(parameter_a*10)), '80']
		return os.spawnvp(os.P_WAIT, modes_path, wl)

	
	def parse_modes_results(self, splat_id, outfname, index2no, graph, codense2db_instance, curs):
		"""
		02-24-05
			parse the modes output file
		"""
		listOfMclResult = []
		reader = csv.reader(file(outfname), delimiter='\t')
		for row in reader:
			unit = MclResult()
			unit.splat_id = splat_id
			mcl_id = int(row[0])
			#no_of_nodes = int(row[1])
			#no_of_edges = float(row[2])		#02-25-05	I saw connectivity =0 in the database with real connectivity=0.5. Some errors in this.
			node_list = row[3:]
			node_list = map(int, node_list)
			subgraph = graph.subgraph_from_node_list(node_list)
			no_of_edges = float(len(subgraph.edges))
			no_of_nodes = float(len(subgraph.nodes))
			unit.connectivity = 2*no_of_edges/(no_of_nodes*(no_of_nodes-1))
			#map the index back to gene_no
			for i in range(len(node_list)):
				node_list[i] = index2no[node_list[i]]
			unit.vertex_set = node_list
			edge_set = []
			#same mapping for the edges
			for edge_id in subgraph.edges:
				edge = subgraph.edge_by_id(edge_id)
				index1 = index2no[edge[0]]
				index2 = index2no[edge[1]]
				if index1 < index2:
					edge_set.append([index1, index2])
				else:
					edge_set.append([index2, index1])
			#calculate the recurrence_array via codense2db_instance's functions
			combined_cor_vector = codense2db_instance.get_combined_cor_vector(curs, edge_set)
			cor_cut_off = 0		#0 means no cut off for those edges.
			unit.recurrence_array = codense2db_instance.parse_recurrence(combined_cor_vector, \
				len(edge_set), cor_cut_off)
			listOfMclResult.append(unit)
		del reader
		return listOfMclResult
		
	def mcl_table_create(self, curs, mcl_table):
		try:
			curs.execute("create table %s(\
				mcl_id	serial primary key,\
				splat_id	integer,\
				vertex_set	integer[],\
				parameter	varchar,\
				connectivity	float,\
				p_value_min	float,\
				go_no_vector	integer[],\
				unknown_gene_ratio	float,\
				recurrence_array	float[])"%mcl_table)
		except:
			sys.stderr.write("Error occurred when creating table %s\n"%mcl_table)	

	def submit(self, curs, mcl_table, listOfMclResult):
		"""
		02-24-05
			submit the cluster to a mcl_result-like table
		"""
		
		try:
			for mclResult in listOfMclResult:
				curs.execute("insert into %s(splat_id, vertex_set, connectivity, recurrence_array)\
						values (%d, ARRAY%s, %f, ARRAY%s)"%\
						(mcl_table, mclResult.splat_id, repr(mclResult.vertex_set),\
						mclResult.connectivity, repr(mclResult.recurrence_array)) )
		except:
			sys.stderr.write('Error occurred when inserting pattern. Aborted.\n')
			sys.exit(1)
		
	def run(self):
		"""
		02-24-05
			--init
			--db_connect
			--splat2graph_dict
			--crack_dict
			--submit
		
		"""
		#some additional initialization
		self.init()
		
		(conn, curs) = self.db_connect(self.hostname, self.dbname, self.schema)
		curs.execute("begin")
		if self.new_table:
			self.mcl_table_create(curs, self.mcl_table)
		curs.execute("DECLARE crs CURSOR FOR select splat_id,edge_set from %s"%self.table)
		curs.execute("fetch 2000 from crs")
		rows = curs.fetchall()
		no_of_splats = 0
		no_of_mcl_results = 0
		while rows:
			for row in rows:
				splat_id = row[0]
				graph_dict = self.splat2graph_dict(row[1])
				#cracking here
				listOfMclResult = self.crack_dict[self.type](splat_id, graph_dict, self.tmpinfname, \
					self.tmpoutfname, self.parameter_a, self.parameter_b, \
					argument1=self.clustering_test_instance, argument2=self.codense2db_instance, argument3=curs)
				no_of_splats += 1
				no_of_mcl_results += len(listOfMclResult)
				if self.needcommit:
					self.submit(curs, self.mcl_table, listOfMclResult)
			if self.report:
				string = 'Splat_patterns: %d, Mcl_results: %d'%(no_of_splats, no_of_mcl_results)
				sys.stderr.write('%s%s'%('\x08'*80,string))
			curs.execute("fetch 2000 from crs")
			rows = curs.fetchall()
		if self.needcommit:
			curs.execute("end")
		#02-25-05 cleanup the temp stuff
		os.remove(self.tmpinfname)
		os.remove(self.tmpoutfname)
		os.rmdir(self.dir_files)
		
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:m:f:a:b:y:ncr", ["help", "hostname=", \
			"dbname=", "schema=", "table=", "mcl_table=", "dir_files=", "parameter_a=",\
			"parameter_b=", "type=", "new_table", "commit", "report"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = 'splat_result'
	mcl_table = 'mcl_result'
	dir_files = '/tmp/CrackSplat'
	parameter_a = 0.2
	parameter_b = 4
	type = 1
	new_table = 0
	commit = 0
	report = 0
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
		elif opt in ("-f", "--dir_files"):
			dir_files = arg
		elif opt in ("-a", "--parameter_a"):
			parameter_a = float(arg)
		elif opt in ("-b", "--parameter_b"):
			parameter_b = int(arg)
		elif opt in ("-y", "--type"):
			type = int(arg)
		elif opt in ("-n", "--new_table"):
			new_table = 1
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-b", "--debug"):
			debug = 1
	if schema:
		instance = CrackSplat(hostname, dbname, schema, table, mcl_table, dir_files, parameter_a,\
			parameter_b, type, new_table, commit, report, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
