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
	-n ..., --gene_p_table=...	the table which store the id of real predictions
	-g ..., --gene_table=...	table storing the stat results, p_gene(default)
	-f ..., --function=...	a go_no
	-o ..., --functioncolor=...	'green' or 'red'
	-c ..., --centralnode=...	a gene_no
	-l ..., --mcl_id=...	the id corresponding to a mcl_cluster
	-y ..., --type=...	the type, 1(single graph), 2 (meta_graph), 3
	-p ..., --plot_type=...	dot(default), neato or twopi
	-h, --help              show this help

Examples:
	subgraph_visualize.py -k sc_yh60_splat -t splat_result_sup_3 -m 
		mcl_result_sup_3_2
		-f 45 -c 4096 -l 27 -y 1 test.R
	
	subgraph_visualize.py -y 2 -k sc_54 -t splat_result_p3g5e6d4q5n80
		-m mcl_result_repos_2 -g  p_gene_repos_2_e5 -n gene_p_repos_2_e5
		-f 202 -c 923 test.R
	
	subgraph_visualize.py -k sc_54 -t splat_result_1 -m mcl_result_1
		-l 0 -y 3 /tmp/test.R
	
Description:
	This program is used to visualize the MCL clusters via R(bioconductor's
	Rgraphviz). Two types of visualization:
	1.
	single graph needs the table, mcl_table, mcl_id
	2.
	meta_graph needs table, mcl_table, gene_p_table, gene_table, centralnode and function
	3.
	display codense sub_subgraph dataset by dataset, need table, mcl_table,
	mcl_id, edge_table(hidden parameter)
	
"""

import sys, os, psycopg, getopt
from graphlib import Graph
from sets import Set
from rpy import r
from codense.common import db_connect

class subgraph_visualize:
	"""
	run
		--dstruc_loadin
		--[type==1]
		--get_subgraph
		--subgraph_output
		
		--[type==2]
		--context_subgraph
			--get_subgraph
		--subgraph_output

		--[type==3]
		--subgraph_in_one_dataset
			--get_subgraph
		--subgraph_output
	
	02-20-05
	unweighted_subgraph_output()
	subgraph_output()
		offer choices on how to draw clusters, "dot", "neato" or "twopi"
	
	02-20-05
	context_subgraph()
		change the way to mcl_id_list according to the change in gene_stat_plot.py's submit()
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, table=None, mcl_table=None,\
			gene_p_table=None, gene_table=None, function=0, functioncolor='green', centralnode=1, mcl_id=1, \
			type=1, r_fname=None, plot_type="dot"):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.table = table
		self.mcl_table = mcl_table
		self.gene_p_table = gene_p_table
		self.gene_table = gene_table
		self.function = int(function)
		self.functioncolor = functioncolor
		self.centralnode = int(centralnode)
		self.mcl_id = int(mcl_id)
		self.type = int(type)
		self.r_fname = r_fname
		self.plot_type = plot_type
		
		#the table for edge_correlation_vector
		self.edge_table = 'edge_cor_vector'
		#mapping between go_no and go_id
		self.go_no2go_id = {}
		#mapping between go_no and go's name
		self.go_no2go_name = {}
		#mapping between gene_no and gene_id
		self.gene_no2gene_id = {}
		self.gene_id2gene_no = {}
		self.global_gene_to_go_dict = {}
		
	def dstruc_loadin(self, curs):
		'''
		'''
		sys.stderr.write("Loading Data STructure...\n")
		
		from codense.common import get_go_no2go_id, get_gene_no2gene_id, get_go_no2name, get_gene_id2gene_no, get_gene_no2go_no
		self.go_no2go_id = get_go_no2go_id(curs)
		self.go_no2go_name = get_go_no2name(curs)
		self.gene_no2gene_id = get_gene_no2gene_id(curs)
		self.gene_id2gene_no = get_gene_id2gene_no(curs)
		self.global_gene_to_go_dict = get_gene_no2go_no(curs)

		curs.execute("select gene_no,go_functions from gene")
		
		if self.type == 3:
			curs.execute("select array_upper(recurrence_array,1) from %s limit 1"%self.table)
			rows = curs.fetchall()
			self.no_of_datasets = int(rows[0][0])
			
		sys.stderr.write("Done\n")
	
	def get_subgraph(self, curs, splat_table, mcl_table, mcl_id, edge_table=None):
		"""
		03-09-05
			make it class independent
			
			input: a mcl_id
			output: a graph with weight
		"""
		sys.stderr.write("Getting subgraph for mcl_id %s..."%mcl_id)
		#first get the mcl cluster information (node_list)
		curs.execute("select splat_id, vertex_set, p_value_min, \
			array_upper(recurrence_array,1) from %s where mcl_id=%d"%(mcl_table, mcl_id))
		rows = curs.fetchall()
		for row in rows:
			splat_id = row[0]
			vertex_set = row[1][1:-1]
			vertex_set = vertex_set.split(',')
			vertex_set = map(int, vertex_set)
			p_value_min = row[2]
			recurrence = row[3]
		
		#second, get the splat pattern from splat_result based on splat_id (big_graph)
		curs.execute("select edge_set from %s where splat_id=%d"%(splat_table, splat_id))
		rows = curs.fetchall()
		for row in rows:
			edge_set = row[0]

		big_graph = Graph.Graph()
		edge_list = edge_set[2:-2].split('},{')
		for edge in edge_list:
			vertex_list = edge.split(',')
			#order it
			vertex_list = map(int, vertex_list)
			vertex_list.sort()
			#set the default edge_data to be 1
			edge_data = 1
			#if edge_table available, get the correlation vector for an edge
			if edge_table:
				#03-31-05 sig_vector replaces cor_vector
				curs.execute("select sig_vector from %s where edge_name='{%s,%s}'"%(edge_table, vertex_list[0], vertex_list[1]))
				rows = curs.fetchall()
				edge_data = rows[0][0][1:-1]
				edge_data = edge_data.split(',')
				edge_data = map(float, edge_data)
			big_graph.add_edge(vertex_list[0], vertex_list[1], edge_data)
			
		subgraph = big_graph.subgraph_from_node_list(vertex_set)
		sys.stderr.write("Done\n")
		return subgraph
	
	def subgraph_in_one_dataset(self, curs, splat_table, mcl_table, edge_table, mcl_id, dataset_no):
		"""
		for type 3
		"""
		sys.stderr.write("Getting subgraph for mcl_id %s in dataset %s..."%(mcl_id, dataset_no))
		subgraph = self.get_subgraph(curs, splat_table, mcl_table, mcl_id, edge_table)
		sub_subgraph = Graph.Graph()
		for edge_id in subgraph.edge_list():
			edge = subgraph.edge_by_id(edge_id)
			edge_data = subgraph.edges[edge_id][2]
			sub_subgraph.add_edge(edge[0], edge[1], edge_data[dataset_no])
		sys.stderr.write("Done\n")
		return sub_subgraph
	
	
	def unweighted_subgraph_output(self, output_f, subgraph, label_dict, gene_no2go_no, centralnode=1, function=0, functioncolor='green', plot_type='dot', weighted=0):
		'''
		write the R script to draw an unweighted subgraph
		03-09-05
			defunct, use subgraph_output() instead.
		
		'''
		output_f.write('library("Rgraphviz")\n')
		vertex_set = subgraph.node_list()
		vertex_labels = []
		for vertex in vertex_set:
			vertex_labels.append('"%s"'%label_dict[vertex])
		output_f.write('V <- c(%s)\n'%(','.join(vertex_labels)))
		output_f.write('edL2 <- vector("list", length=%d)\n'%(len(vertex_set)))
		output_f.write("names(edL2) <- V\n")
		for i in range(len(vertex_set)):
			vertex = vertex_set[i]
			nbrs = subgraph.inc_nbrs(vertex) + subgraph.out_nbrs(vertex)
			nbr_list = []
			for neighbor in nbrs:
				nbr_list.append(vertex_set.index(neighbor)+1)
			output_f.write('edL2[[%d]] <- list(edges=c(%s))\n'%((i+1), ','.join(map(repr,nbr_list))))
		
		output_f.write('gR2 <- new("graphNEL", nodes=V, edgeL=edL2, edgemode="undirected")\n')
		output_f.write('nAttrs = list()\n')
		#output_f.write('eAttrs = list()\n')
		output_f.write("defAttrs = getDefaultAttrs()\n")
		output_f.write('defAttrs$node$color <- "black"\n')
		output_f.write('defAttrs$node$fillcolor <- "transparent"\n')
		output_f.write('defAttrs$node$shape <- "ellipse"\n')
		color_list = []
		central_label = None
		for vertex in vertex_set:
			gene_id = label_dict[vertex]
			if function in gene_no2go_no[vertex]:
				color_list.append('"%s"="%s"'%(gene_id, functioncolor))
			if vertex==centralnode:
				central_label = label_dict[centralnode]
		output_f.write("nAttrs$color <- c(%s)\n"%(','.join(color_list)))
		if central_label:
			output_f.write('nAttrs$fillcolor <- c("%s"="yellow")\n'%central_label)
			output_f.write('nAttrs$shape <- c("%s"="box")\n'%central_label)
		output_f.write('plot(gR2, attrs=defAttrs, nodeAttrs=nAttrs, "%s")\n'%plot_type)
	
	def subgraph_output(self, output_f, subgraph, label_dict, gene_no2go_no, centralnode=1, function=0, functioncolor='green', plot_type='dot', weighted=1):
		'''
		write the R script to draw a weighted subgraph
		
		03-09-05
			make it class-independent
			
			Not giving the centralnode or function is ok.
			
			change name from weighted_subgraph_output to subgraph_output
			
			add a parameter to specify weighted or not
		'''
		sys.stderr.write("Outputing subgraph...")
		output_f.write('library("Rgraphviz")\n')
		vertex_set = subgraph.node_list()
		vertex_labels = []
		for vertex in vertex_set:
			vertex_labels.append('"%s"'%label_dict[vertex])
		output_f.write('V <- c(%s)\n'%(','.join(vertex_labels)))
		output_f.write('edL2 <- vector("list", length=%d)\n'%(len(vertex_set)))
		output_f.write("names(edL2) <- V\n")
		for i in range(len(vertex_set)):
			vertex = vertex_set[i]
			nbrs = subgraph.inc_nbrs(vertex) + subgraph.out_nbrs(vertex)
			nbr_list = []
			weight_list = []
			for neighbor in nbrs:
				if vertex<neighbor:
					edge = subgraph.edge_by_node(vertex, neighbor)
				else:
					edge = subgraph.edge_by_node(neighbor, vertex)
				weight_list.append(subgraph.edges[edge][2])
				#Note +1 to index
				nbr_list.append(vertex_set.index(neighbor)+1)
			if weighted:
				#Note: +1 to i
				output_f.write('edL2[[%d]] <- list(edges=c(%s), weights=c(%s))\n'%((i+1), ','.join(map(repr,nbr_list)), ','.join(map(repr,weight_list)) ))
			else:
				output_f.write('edL2[[%d]] <- list(edges=c(%s))\n'%((i+1), ','.join(map(repr,nbr_list))))
		
		output_f.write('gR2 <- new("graphNEL", nodes=V, edgeL=edL2, edgemode="undirected")\n')
		output_f.write('nAttrs = list()\n')
		output_f.write('eAttrs = list()\n')
		if weighted:
			#this block adds weights to the edges, copied from 'Rgraphviz.pdf'
			output_f.write('ew <- edgeWeights(gR2)\n')
			output_f.write('lw <- unlist(unlist(ew))\n')
			output_f.write('toRemove <- removedEdges(gR2)\n')
			output_f.write('if (length(toRemove) > 0) lw <- lw[-toRemove]\n')
			output_f.write('names(lw) <- edgeNames(gR2)\n')
			output_f.write('eAttrs$label <- lw\n')
			
		output_f.write("defAttrs = getDefaultAttrs()\n")
		output_f.write('defAttrs$node$color <- "black"\n')
		output_f.write('defAttrs$node$fillcolor <- "transparent"\n')
		output_f.write('defAttrs$node$shape <- "ellipse"\n')
		color_list = []
		central_label = None
		for vertex in vertex_set:
			label = label_dict[vertex]
			if vertex in gene_no2go_no:
				if function in gene_no2go_no[vertex]:
					color_list.append('"%s"="%s"'%(label, functioncolor))
			if vertex==centralnode:
				central_label = label_dict[centralnode]
		output_f.write("nAttrs$color <- c(%s)\n"%(','.join(color_list)))
		if central_label:
			output_f.write('nAttrs$fillcolor <- c("%s"="yellow")\n'%central_label)
			output_f.write('nAttrs$shape <- c("%s"="box")\n'%central_label)
		output_f.write('plot(gR2, attrs=defAttrs, nodeAttrs=nAttrs, edgeAttrs=eAttrs, "%s")\n'%plot_type)
		sys.stderr.write("Done\n")
		
	def context_subgraph(self, curs, splat_table, mcl_table, gene_p_table, gene_table, gene_no, go_no):
		"""
		03-09-05
			make it class independent
			
			input: a prediction determined by gene_no and go_no
			output: a meta_graph of all clusters supporting this prediction
		
		03-09-05
			use the gene_p_table to check the real predictions
		"""
		sys.stderr.write("Getting context_subgraph for gene %s and go %s...\n"%(gene_no, go_no))
		curs.execute("select p.mcl_id from %s g, %s p where p.gene_no=%d and \
			p.go_no=%d and g.p_gene_id=p.p_gene_id"%(gene_p_table, gene_table, gene_no, go_no))
		rows = curs.fetchall()
		mcl_id_list = []
		for row in rows:
			mcl_id_list.append(row[0])
		
		if mcl_id_list == []:
			sys.stderr.write("No MCL clusters associated with gene %s and go %s\n"%(gene_no, go_no))
			sys.exit(1)
		meta_graph = Graph.Graph()
		for mcl_id in mcl_id_list:
			subgraph = self.get_subgraph(curs, splat_table, mcl_table, mcl_id)
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
		sys.stderr.write("Done\n")
		return meta_graph
		
	
	def run(self):
		"""
		03-09-05
		"""
		(conn, curs) = db_connect(self.hostname, self.dbname, self.schema)
		self.dstruc_loadin(curs)
		if self.r_fname:
			self.r_f = open(self.r_fname, 'w')
		if self.type == 1:
			subgraph = self.get_subgraph(curs, self.table, self.mcl_table, self.mcl_id)
			#unweighted
			weighted=0
			self.subgraph_output(self.r_f, subgraph, self.gene_no2gene_id, self.global_gene_to_go_dict, \
				self.centralnode, self.function, self.functioncolor, self.plot_type, weighted)
			self.r_f.close()
			r.source(self.r_fname)
			raw_input("Pause:\t")
		elif self.type == 2:
			if self.gene_table==None or self.gene_p_table==None:
				sys.stderr.write("Error: Please specify both the gene_p_table and gene_table.\n")
				sys.exit(2)
			subgraph = self.context_subgraph(curs, self.table, self.mcl_table, self.gene_p_table, self.gene_table, \
				self.centralnode, self.function)
			self.subgraph_output(self.r_f, subgraph, self.gene_no2gene_id, self.global_gene_to_go_dict, \
				self.centralnode, self.function, self.functioncolor, self.plot_type)
			self.r_f.close()
			r.source(self.r_fname)
			raw_input("Pause:\t")
		elif self.type == 3:
			for i in range(self.no_of_datasets):
				sys.stdout.write("Dataset %s\n"%(i+1))
				if self.edge_table==None:
					sys.stderr.write("Error: Please specify both the edge_table.\n")
					sys.exit(2)
					
				sub_subgraph = self.subgraph_in_one_dataset(curs, self.table, self.mcl_table, self.edge_table, self.mcl_id, i)
				self.subgraph_output(self.r_f, sub_subgraph, self.gene_no2gene_id, self.global_gene_to_go_dict, \
					self.centralnode, self.function, self.functioncolor, self.plot_type)
				
				self.r_f.close()
				r.source(self.r_fname)
				#asking  to continue or not
				no_stop = raw_input("Continue? Y/n:\t")
				if no_stop == 'n' or no_stop == 'N':
					sys.exit(3)
				#open it again for the next dataset
				self.r_f = open(self.r_fname, 'w')

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema=", "table=", "mcl_table=", "gene_p_table=", \
		"gene_table=", "function=", "functioncolor=", "centralnode=", "mcl_id=", "type=", "plot_type="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:m:n:g:f:o:c:l:y:p:", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = None
	table = None
	mcl_table = None
	gene_p_table = None
	gene_table = None
	function = 0
	functioncolor = 'green'
	centralnode = 0
	mcl_id = 0
	type = 0
	plot_type = "dot"
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
		elif opt in ("-n", "--gene_p_table="):
			gene_p_table = arg
		elif opt in ("-g", "--gene_table"):
			gene_table = arg
		elif opt in ("-f", "--function"):
			function = int(arg)
		elif opt in ("-o", "--functioncolor"):
			functioncolor = arg
		elif opt in ("-c", "--centralnode"):
			centralnode = int(arg)
		elif opt in ("-l", "--mcl_id"):
			mcl_id = int(arg)
		elif opt in ("-y", "--type"):
			type = int(arg)
		elif opt in ("-p", "--plot_type"):
			plot_type = arg		
	if schema and table and mcl_table and len(args)==1 and type:
		instance = subgraph_visualize(hostname, dbname, schema, table, mcl_table, gene_p_table,\
			gene_table, function, functioncolor, centralnode, mcl_id, type, args[0], plot_type)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
