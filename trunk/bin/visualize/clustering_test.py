#!/usr/bin/env python
"""
Usage: clustering_test.py -k SCHEMA [OPTION] R_FILE

Option:
	R_FILE is the file to store the R code.
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	splat_result table(default)
	-m ..., --mcl_table=...	mcl_result(default), mcl_result corresponding to table above
	-g ..., --gene_table=...	table storing the stat results, p_gene(default)
	-i ..., --input_file=...	gspan format
	-s ..., --cluster_file=...	modes' output
	-o ..., --output_file=...	haiyan's matrix format
	-b ..., --label=...	use gene 1(index, default) or 2(no) or 3(id) to label
	-p ..., --plot_type=...	dot(default), neato or twopi
	-f ..., --function=...	a go_no
	-c ..., --functioncolor=...	'green'(default) or 'red'
	-n ..., --centralnode=...	a gene_no
	-l ..., --mcl_id=...	the id corresponding to a mcl_cluster
	-y ..., --type=...	the type, 1(reformat), 2 (visualize)
	-h, --help              show this help

Examples:
	clustering_test.py -i sc_gspan -y 1 -o matrix /tmp/test.R
	
	clustering_test.py -i sc_gspan -y 2 /tmp/test.R
	
	clustering_test.py -i sc_gspan -y 2 -c myoutput /tmp/test.R
	
Description:
	This program is used for two purposes:
	1. reformat gspan's format into haiyan's matrix format
	2. visualize the input graph, even further label the
	nodes from clusters outputted by haiyan's MODES
	
"""

import sys, os, psycopg, getopt, csv
from graphlib import Graph
from sets import Set
from rpy import r
		
class clustering_test:
	"""
	run	--dstruc_loadin
		--[type==1]
		--reformat	--get_weight
		
		--[type==2]
		--draw_graph
		--visualize_clusters		--draw_graph
	"""
	def __init__(self, hostname, dbname, schema, table, mcl_table,\
			gene_table, input_file, cluster_file, output_file, label, plot_type,\
			function, functioncolor, centralnode, mcl_id, type, r_fname):
		"""
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		"""
		self.table = table
		self.mcl_table = mcl_table
		self.gene_table = gene_table
		self.input_file = input_file
		self.cluster_file = cluster_file
		self.output_file = output_file
		self.label = int(label)
		self.plot_type = plot_type
		self.function = int(function)
		self.functioncolor = functioncolor
		self.centralnode = int(centralnode)
		self.mcl_id = int(mcl_id)
		self.type = int(type)
		self.r_fname = r_fname
		
		#graph got from input_file
		self.input_graph = Graph.Graph()
		#the table for edge_correlation_vector
		self.edge_table = 'edge_cor_vector'
		#mapping between go_no and go_id
		self.go_no2go_id = {}
		#mapping between go_no and go's name
		self.go_no2go_name = {}
		#mapping between gene_no and gene_id
		self.gene_index2no = {}
		self.gene_no2index = {}
		self.gene_index2id = {}
		self.global_gene_to_go_dict = {}
		
		self.label_dict = {}
		self.no_of_genes = 0
		
	def dstruc_loadin(self):
		'''
		'''
		sys.stderr.write("Loading Data STructure...")
		reader = csv.reader(file(self.input_file), delimiter=' ')
		self.no_of_genes = 0
		for row in reader:
			if row[0] == 'v':
				gene_no = int(row[1])
				gene_id = row[2]
				self.gene_index2no[self.no_of_genes] = gene_no
				self.gene_no2index[gene_no] = self.no_of_genes
				self.gene_index2id[self.no_of_genes] = gene_id
				self.input_graph.add_node(self.no_of_genes)
				self.no_of_genes += 1
			elif row[0] == 'e':
				gene_no1 = int (row[1])
				gene_no2 = int(row[2])
				weight = row[3]
				gene_index1 = self.gene_no2index[gene_no1]
				gene_index2 = self.gene_no2index[gene_no2]
				if gene_index1 < gene_index2:
					self.input_graph.add_edge(gene_index1, gene_index2, weight)
				else:
					self.input_graph.add_edge(gene_index2, gene_index1, weight)

		gene_index_list = range(self.no_of_genes)
		self.label_dict = {1: gene_index_list,
			2: self.gene_index2no,
			3: self.gene_index2id}
		del reader	
		sys.stderr.write("Done\n")
	
	def get_weight(self, gene_index1, gene_index2):
		if gene_index1<gene_index2:
			edge_id = self.input_graph.edge_by_node(gene_index1, gene_index2)
		else:
			edge_id = self.input_graph.edge_by_node(gene_index2, gene_index1)
		weight = self.input_graph.edges[edge_id][2]
		return weight
	
	def reformat(self):
		"""
		reformat a graph file in gspan formt to a haiyan's matrix format
		"""
		writer = csv.writer(open(self.output_file, 'w'), delimiter = '\t')
		for i in range(self.no_of_genes):
			nbrs = Set(self.input_graph.all_nbrs(i))
			row = []
			for j in range(self.no_of_genes):
				if j in nbrs:
					weight = self.get_weight(i, j)
					row.append(weight)
				else:
					#this includes i itself and other unlinked nodes
					row.append(0)
			writer.writerow(row)
		del writer
	
	def visualize_clusters(self, cluster_file):
		reader = csv.reader(file(cluster_file), delimiter='\t')
		for row in reader:
			cluster_id = int(row[0])
			node_list = row[3:]
			node_list = map(int, node_list)
			sys.stdout.write("Cluster %s \n"%(cluster_id))
			self.draw_graph(node_list)
			r.source(self.r_fname)
			#asking  to continue or not
			no_stop = raw_input("Continue? Y/n:\t")
			if no_stop == 'n' or no_stop == 'N':
			    sys.exit(3)
		del reader
	
	def draw_graph(self, node_list = []):
		'''
		write the R script to draw an unweighted subgraph
		'''
		self.r_f = open(r_fname, 'w')
		self.r_f.write('library("Rgraphviz")\n')
		vertex_set = self.input_graph.node_list()
		vertex_labels = []
		for vertex in vertex_set:
			vertex_labels.append('"%s"'%self.label_dict[self.label][vertex])
		self.r_f.write('V <- c(%s)\n'%(','.join(vertex_labels)))
		self.r_f.write('edL2 <- vector("list", length=%d)\n'%(len(vertex_set)))
		self.r_f.write("names(edL2) <- V\n")
		for i in range(len(vertex_set)):
			vertex = vertex_set[i]
			nbrs = self.input_graph.all_nbrs(vertex)
			nbr_list = []
			for neighbor in nbrs:
				nbr_list.append(vertex_set.index(neighbor)+1)
			self.r_f.write('edL2[[%d]] <- list(edges=c(%s))\n'%((i+1), ','.join(map(repr,nbr_list))))
		
		self.r_f.write('gR2 <- new("graphNEL", nodes=V, edgeL=edL2, edgemode="undirected")\n')
		self.r_f.write('nAttrs = list()\n')
		#self.r_f.write('eAttrs = list()\n')
		self.r_f.write("defAttrs = getDefaultAttrs()\n")
		self.r_f.write('defAttrs$node$color <- "black"\n')
		self.r_f.write('defAttrs$node$fillcolor <- "transparent"\n')
		self.r_f.write('defAttrs$node$shape <- "ellipse"\n')
		color_list = []
		for vertex in node_list:
			label = self.label_dict[self.label][vertex]
			color_list.append('"%s"="%s"'%(label, self.functioncolor))
		self.r_f.write("nAttrs$color <- c(%s)\n"%(','.join(color_list)))
		self.r_f.write('plot(gR2, attrs=defAttrs, nodeAttrs=nAttrs, "%s")\n'%self.plot_type)
		self.r_f.close()	

	def run(self):
		self.dstruc_loadin()
		if self.type == 1:
			if self.output_file:
				self.reformat()
			else:
				sys.stderr.write("Please the output file.\n")
				sys.exit(2)
		elif self.type == 2:
			if self.r_fname:
				if self.cluster_file:
					self.visualize_clusters(self.cluster_file)				
		
				else:
					self.draw_graph()
					r.source(self.r_fname)
			else:
				sys.stderr.write("Please give the R_FILE.\n")

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema=", "table=", "mcl_table=", \
		"gene_table=",  "input_file=", " cluster_file=", "output_file=", "label=", "plot_type=", \
		"function=", "functioncolor=", "centralnode=", "mcl_id=", "type="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:m:g:i:s:o:b:p:f:c:n:l:y:", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = 'splat_result'
	mcl_table = 'mcl_result'
	gene_table = ''
	input_file = None
	cluster_file = None
	output_file = None
	label = 1
	plot_type = "dot"
	function = 0
	functioncolor = 'green'
	centralnode = 0
	mcl_id = 0
	type = 1
	r_fname = None
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
		elif opt in ("-g", "--gene_table"):
			gene_table = arg
		elif opt in ("-i", "--input_file"):
			input_file = arg
		elif opt in ("-s", "--cluster_file"):
			cluster_file = arg
		elif opt in ("-o", "--output_file"):
			output_file = arg
		elif opt in ("-b", "--label"):
			label = int(arg)
		elif opt in ("-p", "--plot_type"):
			plot_type = arg
		elif opt in ("-f", "--function"):
			function = int(arg)
		elif opt in ("-c", "--functioncolor"):
			functioncolor = arg
		elif opt in ("-n", "--centralnode"):
			centralnode = int(arg)
		elif opt in ("-l", "--mcl_id"):
			mcl_id = int(arg)
		elif opt in ("-y", "--type"):
			type = int(arg)
	
	if len(args) == 1:
		r_fname = args[0]	
	if input_file:
		instance = clustering_test(hostname, dbname, schema, table, mcl_table,\
			gene_table, input_file, cluster_file, output_file, label, plot_type,\
			function, functioncolor, centralnode, mcl_id, type, r_fname)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
