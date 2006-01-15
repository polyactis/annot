#!/usr/bin/env python
"""
Usage: fuzzyDense.py -i pattern_id [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database, 'scfim30' (default)
	-i ..., 	input pattern_id
	-j ...,		input pattern_table, pattern_scfim30m3x30bfs(default)
	-e ...,	degree_cut_off, 0.3(default)
	-p ...,	p_value_cut_off, 0.001 (default)
	-o ...,	output_file, fuzzyDense.out(default)
	-a ...,	min ratio of #associated genes for one bs_no vs cluster size, 0(default)
	-t ...,	top_number of scores to be kept, 5(default)
	-x ...,	tax_id, 4932(default, yeast), for get_gene_id2gene_symbol
	-b,	debug version.
	-r,	enable report flag
	-h,	Display the usage infomation.

12-16-05
	counterpart of haiyan's fuzzyDense
12-18-05
	This standalone program gets a pattern from the pattern_table, 
	get the dense part and do 
	hypergeometric p-value calculation. Output in darwin format.
	
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import boost.graph as bgl
from MpiFromDatasetSignatureToPattern import decodeOccurrence
from sets import Set


class fuzzyDense:
	"""
	12-16-05
		unsolved problem: when >=2 vertices have same degree, the choice is arbitrary
	"""
	def __init__(self, edge2encodedOccurrence, debug=0):
		self.edge2encodedOccurrence = edge2encodedOccurrence
		self.debug=int(debug)
	
	def init_graph_from_vertex_set(self, vertex_list):
		g = bgl.Graph()
		vertex_id=g.vertex_property_map('integer')
		anti_vertex_id = {}	#record the reverse mapping
		for v in vertex_list:
			v_descriptor = g.add_vertex()
			vertex_id[v_descriptor] = v
			anti_vertex_id[v] = v_descriptor
		g.vertex_properties['vertex_id'] = vertex_id
		return g, anti_vertex_id
	
	def remove_singleton_vertices(self, graph):
		vertices_to_remove = []
		for v in graph.vertices:
			if graph.in_degree(v) == 0:
				vertices_to_remove.append(v)
		
		if self.debug:
			print "no of singleton vertices to remove", len(vertices_to_remove)
		
		for v in vertices_to_remove:
			graph.remove_vertex(v)
	
	def get_vertex_min_degree(self, graph):
		min_degree = graph.num_vertices()
		vertex_min_degree = None
		for v in graph.vertices:
			if graph.in_degree(v)<min_degree:
				min_degree = graph.in_degree(v)
				vertex_min_degree = v
		return vertex_min_degree, min_degree
	
	def remove_loose_part_of_graph(self, graph, degree_cut_off):
		self.remove_singleton_vertices(graph)
		vertex_min_degree, min_degree = self.get_vertex_min_degree(graph)
		degree_percentage = float(min_degree)/(graph.num_vertices()-1)
		while degree_percentage<degree_cut_off and vertex_min_degree:
			
			if self.debug:
				print "vertex %s removed with min_degree %s, degree_percentage %s"%(graph.vertex_properties['vertex_id'][vertex_min_degree], min_degree, degree_percentage)
			
			graph.clear_vertex(vertex_min_degree)
			graph.remove_vertex(vertex_min_degree)
			vertex_min_degree, min_degree = self.get_vertex_min_degree(graph)
			degree_percentage = float(min_degree)/(graph.num_vertices()-1)
		return degree_percentage
	
	
	def get_core_vertex_set(self, vertex_list, recurrence_array, degree_cut_off):
		"""
		12-16-05
			global structures used:
				--self.edge2encodedOccurrence
		12-18-05
			expand to all datasets
			
			--init_graph_from_vertex_set()
			--decodeOccurrence()
			--remove_loose_part_of_graph()
				--remove_singleton_vertices()
				--get_vertex_min_degree()
		"""
		no_of_datasets = len(recurrence_array)
		#initialize all graphs
		graph_list = [None]*no_of_datasets
		anti_vertex_id_list = [None]*no_of_datasets
		recurrence_set = Set()
		for i in range(no_of_datasets):
			graph_list[i], anti_vertex_id_list[i] = self.init_graph_from_vertex_set(vertex_list)
			if recurrence_array[i] == 1:
				recurrence_set.add(i)
		
		no_of_vertices = len(vertex_list)
		#vertex_list.sort()	#presorted
		#construct graphs for each 'on' dataset
		for i in range(no_of_vertices):
			for j in range(i+1, no_of_vertices):
				edge_tuple = (vertex_list[i], vertex_list[j])
				"""
				if self.debug:
					print "checking", edge_tuple
				"""
				if edge_tuple in self.edge2encodedOccurrence:
					edge_recurrence = decodeOccurrence(self.edge2encodedOccurrence[edge_tuple])	#starting from 1
					"""
					if self.debug:
						print "edge_recurrence", edge_recurrence
					"""
					for k in edge_recurrence:
						index = k-1
						v_descriptor1 = anti_vertex_id_list[index][vertex_list[i]]
						v_descriptor2 = anti_vertex_id_list[index][vertex_list[j]]
						graph_list[index].add_edge(v_descriptor1, v_descriptor2)
		#remove loose part for each graph
		on_dataset_index_ls = [0]*no_of_datasets
		for i in range(no_of_datasets):
			if graph_list[i].num_edges()>1:	#at least the graph has two edges
				degree_percentage = self.remove_loose_part_of_graph(graph_list[i], degree_cut_off)
				if graph_list[i].num_vertices()>=4:	#min graph size
					"""
					if self.debug:
						print "graph %s has %s vertices remaining with degree_percentage: %s."%(i, graph_list[i].num_vertices(), degree_percentage)
					"""
					on_dataset_index_ls[i] = 1	#this dataset should be counted as 'on'
				else:
					on_dataset_index_ls[i] = 0
		
		#find core vertex_list only in those recurrent 'on' datasets
		vertex_id2occurrence = {}
		recurrent_and_on_datasets_ls = []
		on_but_not_recurrent_dataset2vertex_set = {}
		for i in range(no_of_datasets):
			if recurrence_array[i] == 1 and on_dataset_index_ls[i] == 1:
				recurrent_and_on_datasets_ls.append(i)
				for v in graph_list[i].vertices:
					vertex_id = graph_list[i].vertex_properties['vertex_id'][v]
					if vertex_id not in vertex_id2occurrence:
						vertex_id2occurrence[vertex_id]	= 0
					vertex_id2occurrence[vertex_id] += 1
			if recurrence_array[i] == 0 and on_dataset_index_ls[i] == 1:
				on_but_not_recurrent_dataset2vertex_set[i] = Set()
				for v in graph_list[i].vertices:
					vertex_id = graph_list[i].vertex_properties['vertex_id'][v]
					on_but_not_recurrent_dataset2vertex_set[i].add(vertex_id)
		#only vertices in recurrent and 'on' datasets go into core_vertex_set
		core_vertex_set = Set()
		for vertex_id in vertex_id2occurrence:
			if vertex_id2occurrence[vertex_id] == len(recurrent_and_on_datasets_ls):
				core_vertex_set.add(vertex_id)
		
		#find other on datasets from on_but_not_recurrent_dataset2vertex_set
		for dataset_no, vertex_set in on_but_not_recurrent_dataset2vertex_set.iteritems():
			intersection_set = core_vertex_set & vertex_set
			if len(intersection_set)==len(core_vertex_set):
				recurrent_and_on_datasets_ls.append(dataset_no)
		
		core_vertex_ls  = list(core_vertex_set)
		core_vertex_ls.sort()
		recurrent_and_on_datasets_ls.sort()
		
		return core_vertex_ls, recurrent_and_on_datasets_ls


if __name__ == '__main__':
	"""
	12-16-05 give it a test
	"""
	import sys,os,getopt
	sys.path += [os.path.expanduser('~/script/annot/bin')]
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:j:e:p:o:a:t:x:br", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = 'scfim30'
	pattern_id = None
	pattern_table = 'pattern_scfim30m3x30bfs'
	degree_cut_off = 0.3
	p_value_cut_off = 0.001
	output_file = 'fuzzyDense.out'
	ratio_cutoff = 0.0
	top_number = 5
	tax_id = 4932
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
			pattern_id = int(arg)
		elif opt in ("-j",):
			pattern_table = arg
		elif opt in ("-e",):
			degree_cut_off = float(arg)
		elif opt in ("-p",):
			p_value_cut_off = float(arg)
		elif opt in ("-o",):
			output_file = arg
		elif opt in ("-a",):
			ratio_cutoff = float(arg)
		elif opt in ("-t",):
			top_number = int(arg)
		elif opt in ("-x",):
			tax_id = int(arg)
		elif opt in ("-b",):
			debug = 1
		elif opt in ("-r",):
			report = 1
	
	if schema and pattern_id:
		from codense.common import db_connect, get_gene_id2gene_symbol, dict_map, get_dataset_no2desc
		from MpiCrackSplat import MpiCrackSplat
		conn, curs = db_connect(hostname, dbname, schema)
		MpiCrackSplat_instance = MpiCrackSplat()
		edge2encodedOccurrence = {}
		min_sup = 3
		max_sup = 30
		no_of_datasets = MpiCrackSplat_instance.fill_edge2encodedOccurrence(hostname, dbname, \
		schema, edge2encodedOccurrence, min_sup, max_sup)
		
		"""
		vertex_list = [854078, 854317, 854381, 854469, 855539, 855547, 855613, 855730, 855862, 856019, 856226]
		recurrence_array =  [0.0, 0.0, 0.0, 1.0, 0.53333333333333333, 0.0, 0.066666666666666666, \
		0.20000000000000001, 0.33333333333333331, 0.0, 0.26666666666666666, 0.20000000000000001, \
		0.0, 0.0, 0.33333333333333331, 0.0, 0.66666666666666663, 1.0, 0.0, 0.66666666666666663, 1.0, \
		0.0, 0.066666666666666666, 0.0, 0.46666666666666667, 0.0, 0.13333333333333333, \
		0.20000000000000001, 0.26666666666666666, 1.0]
		degree_cut_off = 0.3
		"""
		
		curs.execute("select vertex_set, recurrence_array from %s where id=%s"%(pattern_table, pattern_id))
		rows = curs.fetchall()
		vertex_list, recurrence_array = rows[0]
		vertex_list = vertex_list[1:-1].split(',')
		vertex_list = map(int, vertex_list)
		recurrence_array = recurrence_array[1:-1].split(',')
		recurrence_array = map(float, recurrence_array)
		
		fuzzyDense_instance = fuzzyDense(edge2encodedOccurrence, debug)
		core_vertex_ls, recurrent_and_on_datasets_ls = fuzzyDense_instance.get_core_vertex_set(vertex_list, recurrence_array, degree_cut_off)
		
		from MpiClusterBsStat import MpiClusterBsStat
		MpiClusterBsStat_instance = MpiClusterBsStat()
		gene_no2bs_no_block = MpiClusterBsStat_instance.get_gene_no2bs_no_block(curs)
		gene_no2bs_no_set, bs_no2gene_no_set = MpiClusterBsStat_instance.construct_two_dicts(0, gene_no2bs_no_block)
		from TF_functions import cluster_bs_analysis
		ls_to_return = cluster_bs_analysis(core_vertex_ls, gene_no2bs_no_set, bs_no2gene_no_set, ratio_cutoff, \
			top_number, p_value_cut_off)
		
		gene_id2symbol = get_gene_id2gene_symbol(curs, tax_id)
		dataset_no2desc = get_dataset_no2desc(curs)
		
		dataset_no_desc_ls = []
		for dataset_index in recurrent_and_on_datasets_ls:
			dataset_no = dataset_index +1
			dataset_no_desc_ls.append([dataset_no, dataset_no2desc[dataset_no]])
		
		
		outf = open(output_file, 'w')
		outf.write("out:=[\n")
		for i in range(len(ls_to_return)):
			row = ls_to_return[i]
			score, score_type, bs_no_list, target_gene_no_list, global_ratio, local_ratio, expected_ratio, unknown_ratio = row
			core_vertex_symbol_ls = dict_map(gene_id2symbol, core_vertex_ls)
			bs_no_symbol_list = dict_map(gene_id2symbol, bs_no_list)
			if i == len(ls_to_return)-1:
				outf.write('[{%s},{%s},{%s}]\n'%(repr(core_vertex_symbol_ls)[1:-1], repr(bs_no_symbol_list)[1:-1], repr(dataset_no_desc_ls)[1:-1]))
			else:
				outf.write('[{%s},{%s},{%s}],\n'%(repr(core_vertex_symbol_ls)[1:-1], repr(bs_no_symbol_list)[1:-1], repr(dataset_no_desc_ls)[1:-1]))
		
		outf.write(']:\n')
		
	else:
		print __doc__
		sys.exit(2)
