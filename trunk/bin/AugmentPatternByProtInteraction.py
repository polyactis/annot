#!/usr/bin/env python
"""
Usage: AugmentPatternByProtInteraction.py -k SCHEMA -i xxx -g xx -c xxx -e xxx -p xxx
	-o xxx [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-i ..., input table, either good_cluster or pattern_table
	-o ...,	output table, eg, aug_pi_xxx
	-n ...,	prot_interaction_table, 'mrinal_pi.intact_interaction'(default)
	-y ...,	type of input table, 1(default, good_cluster) or 2(pattern_table)
	-x ...,	tax_id, 9606(default)
	-c,	commit
	-b,	debug version.
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	python2.4 ./script/annot/bin/AugmentPatternByProtInteraction.py -k hs_fim_65
	-i good_cl_hs_fim_65_n2s175_m5x65s4l5_ft2_e5_000001a60
	-o aug_pi_hs_fim_65_n2s175_m5x65s4l5_ft2_e5_000001a60 -c -r
	
Description:
	To find all protein interaction among nodes in the pattern via shortest_path.
"""


import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import sys, os, getopt, csv, math
from codense.common import db_connect, pg_1d_array2python_ls, return_vertex_set_string, return_edge_set_string
from sets import Set
import networkx as nx
from DrawPredGOExptTFCompTF_Patterns import DrawPredGOExptTFCompTF_Patterns

class AugmentPatternByProtInteraction:
	def __init__(self,hostname='zhoudb', dbname='graphdb', schema=None, \
		input_table=None, output_table=None, prot_interaction_table=None, \
		input_type=1, tax_id=9606, need_commit=0, debug=0, report=0):
		"""
		2006-12-18
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_table = input_table
		self.output_table = output_table
		self.prot_interaction_table = prot_interaction_table
		self.input_type = int(input_type)
		self.tax_id = int(tax_id)
		self.need_commit = int(need_commit)
		self.debug = int(debug)
		self.report = int(report)
	
	def augment_pattern_by_prot_interaction(self, vertex_set, prot_interaction_graph):
		g = nx.Graph()
		no_of_nodes = len(vertex_set)
		for m in range(no_of_nodes):
			for n in range(m+1, no_of_nodes):
				u = vertex_set[m]
				v = vertex_set[n]
				if prot_interaction_graph.has_node(u) and prot_interaction_graph.has_node(v):
					shortest_path_list = nx.shortest_path(prot_interaction_graph, u, v)
					if shortest_path_list:	#check the whole shortest path
						for i in range(len(shortest_path_list)-1):
							if not g.has_edge(shortest_path_list[i], shortest_path_list[i+1]):
								g.add_edge(shortest_path_list[i], shortest_path_list[i+1])
		new_vertex_set = []
		for v in g:
			new_vertex_set.append(v)
		new_vertex_set.sort()
		new_edge_set = []
		for (u, v) in g.edges():
			if u<=v:
				edge_tuple = [u, v]
			else:
				edge_tuple = [v, u]
			new_edge_set.append(edge_tuple)
		new_edge_set.sort()
		return new_vertex_set, new_edge_set
	
	def create_aug_pi_table(self, curs, table):
		sys.stderr.write("Creating table %s..."%table)
		curs.execute("create table %s(\
			id	serial,\
			mcl_id	integer primary key,\
			vertex_set	integer[],\
			edge_set	integer[][])"%table)
		sys.stderr.write("done.\n")
	
	def submit_aug_pi_graph(self, curs, table, mcl_id, vertex_set, edge_set):
		if self.debug:
			sys.stderr.write("\nSubmitting augmented graph to %s..."%table)
		curs_sentence = "insert into %s(mcl_id, vertex_set, edge_set) values(%s, '%s', '%s')"%\
			(table, mcl_id, return_vertex_set_string(vertex_set), return_edge_set_string(edge_set))
		curs.execute(curs_sentence)
		if self.debug:
			sys.stderr.write("submission done.\n")
	
	def batch_augment_all_patterns(self, curs, output_table, prot_interaction_graph):
		sys.stderr.write("Augment patterns from %s..."%input_table)
		curs.execute("fetch 1000 from crs")
		rows = curs.fetchall()
		counter = 0
		while rows:
			for row in rows:
				mcl_id, vertex_set = row
				vertex_set = vertex_set[1:-1].split(',')
				vertex_set = map(int, vertex_set)
				new_vertex_set, new_edge_set = self.augment_pattern_by_prot_interaction(vertex_set, prot_interaction_graph)
				self.submit_aug_pi_graph(curs, output_table, mcl_id, new_vertex_set, new_edge_set)
				counter += 1
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, counter))
			curs.execute("fetch 1000 from crs")
			rows = curs.fetchall()
		
		sys.stderr.write("done.\n")
	
	def run(self):
		"""
		2006-12-25
			
		"""
		(conn, curs) = db_connect(self.hostname, self.dbname, self.schema)
		self.create_aug_pi_table(curs, self.output_table)
		if self.input_type==1:
			curs.execute("DECLARE crs CURSOR FOR select mcl_id, vertex_set from %s"%self.input_table)
		elif self.input_type==2:
			curs.execute("DECLARE crs CURSOR FOR select id, vertex_set from %s"%self.input_table)
		DrawPredGOExptTFCompTF_Patterns_instance = DrawPredGOExptTFCompTF_Patterns()
		prot_interaction_graph = DrawPredGOExptTFCompTF_Patterns_instance.get_prot_interaction_graph(\
			curs, self.prot_interaction_table, self.tax_id)
		self.batch_augment_all_patterns(curs, self.output_table, prot_interaction_graph)
		curs.execute("close crs")
		if self.need_commit:
			curs.execute("end")

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:o:n:y:x:cbr", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	input_table = None
	output_table = None
	prot_interaction_table = 'mrinal_pi.intact_interaction'
	input_type = 1
	tax_id = 9606
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
		elif opt in ("-i",):
			input_table = arg
		elif opt in ("-o",):
			output_table = arg
		elif opt in ("-n",):
			prot_interaction_table = arg
		elif opt in ("-y",):
			input_type = int(arg)
		elif opt in ("-x",):
			tax_id = int(arg)
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-b",):
			debug = 1
		elif opt in ("-r",):
			report = 1
	if schema and input_table and output_table:
			instance = AugmentPatternByProtInteraction(hostname, dbname, schema, \
				input_table, output_table, prot_interaction_table, input_type, \
				tax_id, commit, debug, report)
			instance.run()
	else:
		print __doc__
		sys.exit(2)