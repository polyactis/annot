#!/usr/bin/env python
"""
Usage: function_mapping.py -f FUNCTION_MAPPING FILE [OPTION]  INPUT_FILE OUTPUT_FILE

Option:
	INPUT_FILE is the graph file in gspan format.
	OUTPUT_FILE is the file to hold all clusters.
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	the splat_result table (edge_set)
	-f ..., --function_mapping_file=...	the file containing GO:Gene, output by go_informative_node.py
	-s ..., --seed_size_cutoff=...	the minimum size of a seed, 4(default)
	-l ..., --linking_cutoff=...	the minimum number of links between a seed and a node, 2(default)
	-a ..., --association_cut_off=...	the cutoff for a node to be included in the cluster, 0.6(default)
		overlapping/seed_size
	-o ..., --out_association_cut_off=...	the cutoff supplementary to association_cut_off, 0(default)
		overlapping/len(all_nbrs)
	-n ..., --connectivity_cut_off=...	the connectivity_cut_off of a cluster, not used now
	-b, --debug	debug
	-r, --report	report the progress(a number)
	-h, --help		show this help
	
Examples:
	function_mapping.py -f sc_54_info_nodes sc_54_merge_6 sc_54_merge_6.clusters

Description:
	This program implements a graph-clustering algorithm totally biology-driven, named
	function_mapping.

"""
import sys, os, getopt, csv, math, psycopg
from graphlib import Graph
from sets import Set

class function_mapping:
	'''
	A clustering algorithm based on function_mapping.
	'''
	def __init__(self, input_file, output_file, hostname, dbname, schema, table, function_mapping_file, \
			seed_size_cutoff, linking_cutoff,  association_cut_off, out_association_cut_off,\
			connectivity_cut_off, debug, report):
		self.input_file = input_file
		self.output_file = output_file
		"""
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		"""
		self.table = table
		self.function_mapping_file = function_mapping_file
		#the cutoff for a candidate cluster(seed)
		self.seed_size_cutoff = int(seed_size_cutoff)
		#the minimum number of links between a seed and a node, apart from the association_cut_off
		self.linking_cutoff = int(linking_cutoff)
		self.association_cut_off = float(association_cut_off)
		self.out_association_cut_off = float(out_association_cut_off)
		self.connectivity_cut_off = float(connectivity_cut_off)
		self.debug = int(debug)
		self.report = int(report)
		
		#to output the clusters
		self.outf = csv.writer(open(self.output_file, 'w'), delimiter='\t')
		#the initial graph to be clustered.
		self.input_graph = Graph.Graph()
		#the gene_id2no dictionary based on the input_graph
		self.gene_id2no = {}
		#the function association dictionary
		self.go_id2gene_no_list = {}
		#the number of clusters
		self.no_of_clusters = 0
		
	def dstruc_loadin(self):
		sys.stderr.write("Loading Data STructure...")
		reader = csv.reader(file(self.input_file), delimiter=' ')
		for row in reader:
			if row[0] == 'v':
				gene_no = int(row[1])
				gene_id = row[2]
				self.gene_id2no[gene_id] = gene_no
				self.input_graph.add_node(gene_no)
			elif row[0] == 'e':
				gene_no1 = int (row[1])
				gene_no2 = int(row[2])
				weight = row[3]
				if gene_no1 < gene_no2:
					self.input_graph.add_edge(gene_no1, gene_no2, weight)
				else:
					self.input_graph.add_edge(gene_no2, gene_no1, weight)
		del reader
		
		reader = csv.reader(file(self.function_mapping_file), delimiter='\t')
		for row in reader:
			go_id = row[0]
			gene_id = row[2]
			if gene_id in self.gene_id2no:
				gene_no = self.gene_id2no[gene_id]
				if go_id in self.go_id2gene_no_list:
					self.go_id2gene_no_list[go_id].append(gene_no)
				else:
					self.go_id2gene_no_list[go_id] = [gene_no]
		del reader
		
		sys.stderr.write("Done\n")
	
	def component_grow(self, component):
		if self.debug:
			print "\t\t## in function component_grow()"
		#the invariant set is used to judge the association_cut_off
		component_set_invariant = component.copy()
		component_size = len(component)
		for node in component_set_invariant:
			nbrs = self.input_graph.all_nbrs(node)
			for nbr in nbrs:
				if nbr not in component:
					#not in the seed and not included in the component
					#
					#get all the neighbors of the nbr
					nbrs_of_nbr = Set(self.input_graph.all_nbrs(nbr))
					overlapping_nodes = nbrs_of_nbr & component_set_invariant
					if self.debug:
						print "node %s has %s nbrs with %s overlapping_nodes, with component size %s"%\
							(nbr, len(nbrs_of_nbr), len(overlapping_nodes), component_size)
						raw_input("pause:")
					if len(overlapping_nodes) >= self.linking_cutoff:
						#first the number of overlapping_nodes should pass cutoff
						association_ratio = len(overlapping_nodes)/float(component_size)
						out_association_ratio = len(overlapping_nodes)/float(len(nbrs_of_nbr))
						if self.debug:
							print "\tnode %s passed the linking cutoff %s, with association_ratio %s(%s)"%\
								(nbr, self.linking_cutoff, association_ratio, self.association_cut_off)
							print "\tout_association_ratio: %s , cutoff: %s"%(out_association_ratio, self.out_association_cut_off)
						if association_ratio>=self.association_cut_off and out_association_ratio>=self.out_association_cut_off:
							#second, the in-cluster/total ratio should exceed the cutoff
							component.add(nbr)
							if self.debug:
								print "\tnode %s is added to the component, new_size: %s"%(nbr, len(component))
		if self.debug:
			if len(component) == len(component_set_invariant):
				print "this component doesn't grow, size %d"%component_size
			else:
				print "this component grows from %d to %d"%(component_size, len(component))
			print "\t\t##leave function component_grow()"
		return component
	
	
	def get_connected_component(self, gene_no_list):
		if self.debug:
			print "\t\t## in function get_connected_component()"
		subgraph = self.input_graph.subgraph_from_node_list(gene_no_list)
		connected_components = []
		for gene_no in gene_no_list:
			if gene_no in subgraph.nodes:
				#keep it in a set form, easy for grow
				component = Set(subgraph.all_nbrs(gene_no))
				#add gene_no itself
				component.add(gene_no)
				#hide them in the subgraph
				for node in component:
					subgraph.hide_node(node)
				if self.debug:
					print "gene_no: %s"%gene_no
					print "component with size %s: %s"%(len(component), repr(component))
				if len(component)>= self.seed_size_cutoff:
					#not too small
					connected_components.append(component)
					if self.debug:
						print "\t***pass the seed_size_cutoff %d"%self.seed_size_cutoff
					
					cluster = self.component_grow(component)
					self.no_of_clusters += 1
					self.cluster_output(cluster)
		if self.debug:
			print "\t\t## leave function get_connected_component()"
		return connected_components
		
	def cluster_output(self, cluster):
		'''
		output in haiyan's codense format
		'''
		#1st column
		row = [self.no_of_clusters]
		subgraph = self.input_graph.subgraph_from_node_list(cluster)
		connectivity = len(cluster)/float(len(subgraph.edges))
		#2nd column
		row.append(connectivity)
		#3rd column
		node_string = '{'
		for node in cluster:
			node_string += '%s;'%node
		node_string += '}'
		row.append(node_string)
		#4th column
		edge_string = '{'
		for edge_id in subgraph.edges:
			edge = subgraph.edge_by_id(edge_id)
			edge_string += '(%s,%s );'%(edge[0], edge[1])
		edge_string += '}'
		row.append(edge_string)
		self.outf.writerow(row)
		
	def run(self):
		self.dstruc_loadin()
		no_of_functions = 0
		for go_id in self.go_id2gene_no_list:
			no_of_functions += 1
			if self.report:
				print "function: %s, %s/%s"%(go_id, no_of_functions, len(self.go_id2gene_no_list))
			gene_no_list = self.go_id2gene_no_list[go_id]
			connected_components = self.get_connected_component(gene_no_list)
			if self.report:
				print "no of clusters: %s"%self.no_of_clusters
			"""
			#do it in get_connected_component()
			for component in connected_components:
				cluster = self.component_grow(component)
				self.no_of_clusters += 1
				self.cluster_output(cluster)
			"""

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
			
	long_options_list = ["help", "hostname=", "dbname=", "schema=", "table=", "function_mapping_file=", \
		"seed_size_cutoff=",  "linking_cutoff=", "association_cut_off=", "out_association_cut_off=", \
		"connectivity_cut_off=", "debug", "report"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:f:s:l:a:o:n:br", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = 'splat_result'
	function_mapping_file = None
	seed_size_cutoff = 4
	linking_cutoff = 2
	association_cut_off = 0.6
	out_association_cut_off = 0
	connectivity_cut_off = 0.4
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
		elif opt in ("-t", "--table"):
			table = arg
		elif opt in ("-f", "--function_mapping_file"):
			function_mapping_file = arg
		elif opt in ("-s", "--seed_size_cutoff"):
			seed_size_cutoff = int(arg)
		elif opt in ("-l", "--linking_cutoff"):
			linking_cutoff = int(arg)
		elif opt in ("-a", "--association_cut_off"):
			association_cut_off = float(arg)
		elif opt in ("-o", "--out_association_cut_off"):
			out_association_cut_off = float(arg)
		elif opt in ("-n", "--connectivity_cut_off"):
			connectivity_cut_off = float(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	if len(args) == 2 and function_mapping_file:
		instance = function_mapping(args[0], args[1], hostname, dbname, schema, table, function_mapping_file,\
			seed_size_cutoff, linking_cutoff, association_cut_off, out_association_cut_off, connectivity_cut_off, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
