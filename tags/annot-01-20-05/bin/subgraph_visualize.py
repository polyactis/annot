#!/usr/bin/env python
"""
Usage: subgraph_visualize.py -k SCHEMA [OPTION] R_FILE

Option:
	R_FILE is the file to store the R code.
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	splat_result table(default)
	-m ..., --mcl_table=...	mcl_result(default), mcl_result corresponding to table above
	-g ..., --gene_table=...	table storing the stat results, p_gene(default)
	-f ..., --function=...	a go_no
	-c ..., --functioncolor=...	'green' or 'red'
	-n ..., --centralnode=...	a gene_no
	-l ..., --mcl_id=...	the id corresponding to a mcl_cluster
	-y ..., --type=...	the type, 1(single graph), 2 (meta_graph)
	-h, --help              show this help

Examples:
	subgraph_visualize.py -k sc_yh60_splat -t splat_result_sup_3 -m 
		mcl_result_sup_3_2 -f 45 -c 'red' -n 4096 -l 27 test.R
	
	subgraph_visualize.py -k sc_yh60_splat -t splat_result_sup_3 -m 
		mcl_result_sup_3_2 -g p_gene_cluster_stat_sup_3_2_3
		-f 45 -c 'red' -n 4096 test.R
	
Description:
	This program is used to visualize the MCL clusters via R(bioconductor's
	Rgraphviz). Two types of visualization:
	1.
	single graph  only needs the table, mcl_table, mcl_id and centralnode.
	2.
	meta_graph needs table, mcl_table, gene_table, function and centralnode.
	
"""

import sys, os, psycopg, getopt
from graphlib import Graph
from sets import Set
		
class subgraph_visualize:
	def __init__(self, hostname, dbname, schema, table, mcl_table,\
			gene_table, function, functioncolor, centralnode, mcl_id, type, r_fname):
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.table = table
		self.mcl_table = mcl_table
		self.gene_table = gene_table
		self.function = int(function)
		self.functioncolor = functioncolor
		self.centralnode = int(centralnode)
		self.mcl_id = int(mcl_id)
		self.type = int(type)
		self.r_f = open(r_fname, 'w')
		self.r_f.write('library("Rgraphviz")\n')

		#mapping between go_no and go_id
		self.go_no2go_id = {}
		#mapping between go_no and go's name
		self.go_no2go_name = {}
		#mapping between gene_no and gene_id
		self.gene_no2gene_id = {}
		self.gene_id2gene_no = {}
		self.global_gene_to_go_dict = {}
		
	def dstruc_loadin(self):
		'''
		'''
		sys.stderr.write("Loading Data STructure...")
		
		#setup self.go_no2go_id
		self.curs.execute("select go_no, go_id from go")
		rows = self.curs.fetchall()
		for row in rows:
			self.go_no2go_id[row[0]] = row[1]
		
		#setup self.go_no2go_name
		self.curs.execute("select g.go_no,t.name from go g, go.term t where g.go_id=t.acc")
		rows = self.curs.fetchall()
		for row in rows:
			self.go_no2go_name[row[0]] = row[1]
				
		#setup self.gene_no2gene_id
		self.curs.execute("select gene_no, gene_id from gene")
		rows = self.curs.fetchall()
		for row in rows:
			self.gene_no2gene_id[row[0]] = row[1]
			self.gene_id2gene_no[row[1]] = row[0]

		self.curs.execute("select gene_no,go_functions from gene")
			
		rows = self.curs.fetchall()
		for row in rows:
			self.global_gene_to_go_dict[row[0]] = []
			go_functions_list = row[1][1:-1].split(',')
			for go_no in go_functions_list:
				self.global_gene_to_go_dict[row[0]].append(int(go_no))
		
		sys.stderr.write("Done\n")
	
	def get_subgraph(self, mcl_id):
		#first get the mcl cluster information
		self.curs.execute("select splat_id, vertex_set, p_value_min, \
			array_upper(recurrence_array,1) from %s where mcl_id=%d"%(self.mcl_table, mcl_id))
		rows = self.curs.fetchall()
		for row in rows:
			splat_id = row[0]
			vertex_set = row[1][1:-1]
			vertex_set = vertex_set.split(',')
			vertex_set = map(int, vertex_set)
			p_value_min = row[2]
			recurrence = row[3]
		
		#get the splat pattern from splat_result based on splat_id
		self.curs.execute("select edge_set from %s where splat_id=%d"%(self.table, splat_id))
		rows = self.curs.fetchall()
		for row in rows:
			edge_set = row[0]
		
		big_graph = Graph.Graph()
		edge_list = edge_set[2:-2].split('},{')
		for edge in edge_list:
			vertex_list = edge.split(',')
			big_graph.add_edge(int(vertex_list[0]), int(vertex_list[1]), 1)
			
		subgraph = big_graph.subgraph_from_node_list(vertex_set)
		return subgraph
		
	def one_graph(self, gene_no, mcl_id):
		'''
		write the R script to draw an unweighted subgraph
		'''
		subgraph = self.get_subgraph(mcl_id)
		vertex_set = subgraph.node_list()
		vertex_labels = []
		for vertex in vertex_set:
			vertex_labels.append('"%s"'%self.gene_no2gene_id[vertex])
		self.r_f.write('V <- c(%s)\n'%(','.join(vertex_labels)))
		self.r_f.write('edL2 <- vector("list", length=%d)\n'%(len(vertex_set)))
		self.r_f.write("names(edL2) <- V\n")
		for i in range(len(vertex_set)):
			vertex = vertex_set[i]
			nbrs = subgraph.inc_nbrs(vertex) + subgraph.out_nbrs(vertex)
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
		for vertex in vertex_set:
			gene_id = self.gene_no2gene_id[vertex]
			if self.function in self.global_gene_to_go_dict[vertex]:
				color_list.append('%s="%s"'%(gene_id, self.functioncolor))
		self.r_f.write("nAttrs$color <- c(%s)\n"%(','.join(color_list)))
		gene_id = self.gene_no2gene_id[centralnode]
		self.r_f.write('nAttrs$fillcolor <- c(%s="yellow")\n'%gene_id)
		self.r_f.write('nAttrs$shape <- c(%s="box")\n'%gene_id)
		self.r_f.write('plot(gR2, attrs=defAttrs, nodeAttrs=nAttrs)\n')
	
	def weighted_subgraph(self, subgraph):
		'''
		write the R script to draw a weighted subgraph
		'''
		vertex_set = subgraph.node_list()
		vertex_labels = []
		for vertex in vertex_set:
			vertex_labels.append('"%s"'%self.gene_no2gene_id[vertex])
		self.r_f.write('V <- c(%s)\n'%(','.join(vertex_labels)))
		self.r_f.write('edL2 <- vector("list", length=%d)\n'%(len(vertex_set)))
		self.r_f.write("names(edL2) <- V\n")
		for i in range(len(vertex_set)):
			vertex = vertex_set[i]
			nbrs = subgraph.inc_nbrs(vertex) + subgraph.out_nbrs(vertex)
			nbr_list = []
			weight_list = []
			for neighbor in nbrs:
				#first order
				edge = subgraph.edge_by_node(vertex, neighbor)
				if edge == None:
					#second order
					edge = subgraph.edge_by_node(neighbor, vertex)
				weight_list.append(subgraph.edges[edge][2])
				#Note +1 to index
				nbr_list.append(vertex_set.index(neighbor)+1)
			#Note: +1 to i
			self.r_f.write('edL2[[%d]] <- list(edges=c(%s), weights=c(%s))\n'%((i+1), ','.join(map(repr,nbr_list)), ','.join(map(repr,weight_list)) ))
		
		self.r_f.write('gR2 <- new("graphNEL", nodes=V, edgeL=edL2, edgemode="undirected")\n')
		self.r_f.write('nAttrs = list()\n')
		self.r_f.write('eAttrs = list()\n')
		#this block adds weights to the edges, copied from 'Rgraphviz.pdf'
		self.r_f.write('ew <- edgeWeights(gR2)\n')
		self.r_f.write('lw <- unlist(unlist(ew))\n')
		self.r_f.write('toRemove <- removedEdges(gR2)\n')
		self.r_f.write('if (length(toRemove) > 0) lw <- lw[-toRemove]\n')
		self.r_f.write('names(lw) <- edgeNames(gR2)\n')
		self.r_f.write('eAttrs$label <- lw\n')
		
		self.r_f.write("defAttrs = getDefaultAttrs()\n")
		self.r_f.write('defAttrs$node$color <- "black"\n')
		self.r_f.write('defAttrs$node$fillcolor <- "transparent"\n')
		self.r_f.write('defAttrs$node$shape <- "ellipse"\n')
		color_list = []
		for vertex in vertex_set:
			gene_id = self.gene_no2gene_id[vertex]
			if self.function in self.global_gene_to_go_dict[vertex]:
				color_list.append('%s="%s"'%(gene_id, self.functioncolor))
		self.r_f.write("nAttrs$color <- c(%s)\n"%(','.join(color_list)))
		gene_id = self.gene_no2gene_id[centralnode]
		self.r_f.write('nAttrs$fillcolor <- c(%s="yellow")\n'%gene_id)
		self.r_f.write('nAttrs$shape <- c(%s="box")\n'%gene_id)
		self.r_f.write('plot(gR2, attrs=defAttrs, nodeAttrs=nAttrs, edgeAttrs=eAttrs)\n')
		
	def context_subgraph(self, gene_no, go_no):
		self.curs.execute("select cluster_context, cluster_array from %s where gene_no=%d \
			and go_no=%d"%(self.gene_table, gene_no, go_no))
		rows = self.curs.fetchall()
		mcl_id_list = rows[0][1][1:-1]
		mcl_id_list = mcl_id_list.split(',')
		mcl_id_list = map(int, mcl_id_list)
		
		if mcl_id_list == []:
			sys.stderr.write("No MCL clusters associated.\n")
			sys.exit(1)
		meta_graph = Graph.Graph()
		for mcl_id in mcl_id_list:
			subgraph = self.get_subgraph(mcl_id)
			for edge in subgraph.edge_list():
				node1 = subgraph.edges[edge][0]
				node2 = subgraph.edges[edge][1]
				if (node1 in meta_graph.nodes) and (node2 in meta_graph.nodes):
					#first check the first order.
					meta_edge = meta_graph.edge_by_node(node1, node2)
					if meta_edge == None:
						#then second order.
						meta_edge = meta_graph.edge_by_node(node2, node1)
					if meta_edge:
						#increase the weight by 1
						meta_graph.edges[meta_edge][2] += 1
					else:
						#add a new edge
						meta_graph.add_edge(node1, node2, 1)
				else:
					#add a new edge
					meta_graph.add_edge(node1, node2, 1)
		return meta_graph
		
	
	def run(self):
		if self.type == 1:
			#subgraph = self.get_subgraph(self.mcl_id)
			#self.weighted_subgraph(subgraph)
			self.one_graph(self.centralnode, self.mcl_id)
		elif self.type == 2:
		
			subgraph = self.context_subgraph(self.centralnode, self.function)
			self.weighted_subgraph(subgraph)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema=", "table=", "mcl_table=", \
		"gene_table=", "function=", "functioncolor=", "centralnode=", "mcl_id=", "type="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:m:g:f:c:n:l:y:", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = 'splat_result'
	mcl_table = 'mcl_result'
	gene_table = ''
	function = 0
	functioncolor = 'green'
	centralnode = 0
	mcl_id = 0
	type = 0
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
			
	if schema and centralnode and len(args)==1 and gene_table and function and type:
		instance = subgraph_visualize(hostname, dbname, schema, table, mcl_table,\
			gene_table, function, functioncolor, centralnode, mcl_id, type, args[0])
		instance.dstruc_loadin()
		instance.run()
	else:
		print __doc__
		sys.exit(2)
