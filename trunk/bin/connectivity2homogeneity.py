#!/usr/bin/env python
"""
Usage: connectivity2homogeneity.py -k SCHEMA -i INPUT -o OUTPUT [OPTIONS]

Option:
	STAT_TABLE_FILE is the file where the output goes.
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-s ..., --splat_table=...
	-m ..., --mcl_table=...
	-i ..., --input_file=...	input_file, summary_graph in gspan format
	-o ..., --output_fname=...	where column data goes.
	-n ..., --no_of_random_subgraphs=...	1e5 (default)
	-e ..., --random_subgraph_size=...	random generate(default)
	-y ..., --type=...	running type, 1(default)
	-r, --report	report the progress(a number)
	-b, --debug	enable debugging, no debug by default
	-h, --help              show this help

Examples:
	connectivity2homogeneity.py -k sc_54_6661 -i ~/bin/hhu_clustering/data/input/sc_54_6661_7 -e 6 -o /tmp/yh/data

Description:
	This class investigates the relationship between connectivity and function homogeneity.
"""

import sys, os, psycopg, getopt, csv, math, random
from sets import Set
from codense.common import db_connect

class connectivity2homogeneity:
	"""
	04-03-05
		This class investigates the relationship between connectivity and function homogeneity.
	
	--run()
		--db_connect
		--get_summary_graph()
		--get_gene_no2go_no_list()
		--get_random_subgraph()
			--_connectivity2homogeneity()
	"""
	def __init__(self, input_file=None, output_fname=None, hostname='zhoudb', \
		dbname='graphdb', schema=None, splat_table=None, mcl_table=None, \
		no_of_random_subgraphs=1e5, random_subgraph_size=None,\
		type=1, report=0, debug=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema		
		self.splat_table = splat_table
		self.mcl_table = mcl_table
		self.no_of_random_subgraphs = int(no_of_random_subgraphs)
		self.random_subgraph_size = random_subgraph_size
		self.report = int(report)
		self.debug = int(debug)
		self.type = int(type)
		
		self.input_file = input_file
		self.output_fname = output_fname
		

		
	def get_summary_graph(self, input_file, min_weight=None):
		"""
		04-03-05
			read a summary graph from a file
		"""
		from codense.common import get_gspan_graph
		return get_gspan_graph(input_file, min_weight)
	
	def get_gene_no2go_no_list(self, curs, depth=5):
		"""
		04-03-05
			get the mapping between gene_no and its associated functions
			given the depth.(unknown's depth is 2, so all these genes are known genes.
			the fact is used in _connectivity2homogeneity())
		"""
		from codense.common import get_go_no2depth
		go_no2depth = get_go_no2depth(curs)
		if self.debug:
			print "length of go_no2depth is %s"%len(go_no2depth)
			
		#codes below similar to get_gene_no2go_no of codense.common, but different.
		sys.stderr.write("Getting gene_no2go_no (go_no depth:%s) ..."%depth)
		gene_no2go_no = {}
		curs.execute("select gene_no,go_functions from gene")
		rows = curs.fetchall()
		for row in rows:
			go_functions_list = row[1][1:-1].split(',')
			#don't forget to transform the data type to integer.
			go_functions_list = map(int, go_functions_list)
			if self.debug:
				print "gene is %s"%row[0]
				print "go_functions_list is %s"%row[1]
			for go_no in go_functions_list:
				go_no_depth = go_no2depth.get(go_no)
				if self.debug:
					print "go_no %s depth: %s"%(go_no, go_no_depth)
					raw_input("pause:")
				if go_no_depth==depth:
					if row[0] not in gene_no2go_no:
						gene_no2go_no[row[0]] = []
					gene_no2go_no[row[0]].append(go_no)
		sys.stderr.write("Done\n")
		return gene_no2go_no
		
	def get_random_subgraph(self, summary_graph, size=None):
		"""
		04-03-05
			get a random subgraph out of a summary_graph with same size
			or different sizes.
		"""
		sys.stderr.write("Generating random subgraphs and analyzing...\n")
		node_list = summary_graph.node_list()
		for i in range(self.no_of_random_subgraphs):
			if size == None:
				#the cluster size is random too.
				size = random.randint(1, len(node_list))
			node_sublist_index = random.sample(range(len(node_list)), size)
			node_sublist = []
			for index in node_sublist_index:
				node_sublist.append(node_list[index])
			subgraph = summary_graph.subgraph_from_node_list(node_sublist)
			edge_list = subgraph.edge_list()
			edge_set = map(subgraph.edge_by_id, edge_list)
			self._connectivity2homogeneity(i, node_sublist, edge_set, self.gene_no2go_no, self.writer)
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*10, i))
		sys.stderr.write("Done.\n")
	
	def subgraph_fetch(self, curs, mcl_table, splat_table, size=None):
		"""
		04-03-05
			counterpart of get_random_subgraph, get subgraph from database
		"""
	
	def _connectivity2homogeneity(self, subgraph_id, node_list, edge_list, gene_no2go_no, writer):
		"""
		04-03-05
			core function.
			carry out several homogeneity measures on the cluster and output
		"""
		no_of_nodes = len(node_list)
		connectivity = 2*len(edge_list)/(no_of_nodes*(no_of_nodes-1.0))
		
		#codes below copied from cluster_stat.py
		local_go_no_dict = {}
		no_of_knowns = 0
		for gene_no in node_list:
			if gene_no in gene_no2go_no:
				#it's a known gene first
				no_of_knowns += 1
				go_no_list = gene_no2go_no[gene_no]
				for go_no in go_no_list:
					if go_no in local_go_no_dict:
						local_go_no_dict[go_no] += 1
					else:
						local_go_no_dict[go_no] = 1
		if self.debug:
			print "local_go_no_dict is %s"%repr(local_go_no_dict)
		#get the go_no with maximum associations
		association_go_no_list = []
		for go_no in local_go_no_dict:
			association_go_no_list.append([local_go_no_dict[go_no], go_no])
		association_go_no_list.sort()
		if self.debug:
			print "association_go_no_list is %s"%repr(association_go_no_list)
		if len(association_go_no_list)>0:
			#has some known genes
			go_no_with_max_associations = association_go_no_list[-1][1]
			max_associations = association_go_no_list[-1][0]
		else:
			return
			
		
		#two percentages,
		perc2all = max_associations/float(no_of_nodes)
		perc2known = max_associations/float(no_of_knowns)
		#output
		writer.writerow([perc2all,perc2known, connectivity, no_of_nodes, no_of_knowns, \
			go_no_with_max_associations, subgraph_id])
		
	def run(self):
		"""
		04-03-05
		"""
		(conn, curs) = db_connect(self.hostname, self.dbname, self.schema)
		self.summary_graph = self.get_summary_graph(self.input_file, min_weight=None)
		self.gene_no2go_no = self.get_gene_no2go_no_list(curs, depth=5)
		self.writer = csv.writer(open(self.output_fname, 'w'), delimiter='\t')
		#self.gene_no2go_no and self.writer should be ready before calling get_random_subgraph()
		#because _connectivity2homogeneity() needs them.
		self.get_random_subgraph(self.summary_graph, size=self.random_subgraph_size)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema=", "splat_table=", "mcl_table=", \
		"input_file=", "output_fname=", "no_of_random_subgraphs=", "random_subgraph_size=",\
		"type=", "report", "debug"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:s:m:i:o:n:e:y:rb", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	splat_table = None
	mcl_table = None
	input_file = None
	output_fname = None
	no_of_random_subgraphs = 1e5
	random_subgraph_size = None
	type = 1
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
		elif opt in ("-s", "--splat_table"):
			table = arg
		elif opt in ("-m", "--mcl_table"):
			mcl_table = arg
		elif opt in ("-i", "--input_file"):
			input_file = arg
		elif opt in ("-o", "--output_fname"):
			output_fname = arg
		elif opt in ("-n", "--no_of_random_subgraphs="):
			no_of_random_subgraphs = int(arg)
		elif opt in ("-e", "--random_subgraph_size="):
			random_subgraph_size = int(arg)
		elif opt in ("-y", "--type="):
			type = int(arg)
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-b", "--debug"):
			debug = 1
			
	if schema and output_fname:
		instance = connectivity2homogeneity(input_file, output_fname, hostname, dbname, schema, \
			splat_table, mcl_table, no_of_random_subgraphs, random_subgraph_size,\
			type, report, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
