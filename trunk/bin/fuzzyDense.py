#!/usr/bin/env python
"""
12-16-05
	counterpart of haiyan's fuzzyDense
"""
import boost.graph as bgl
from MpiFromDatasetSignatureToPattern import decodeOccurrence

class fuzzyDense:
	"""
	12-16-05
		unsolved problem: when >=2 vertices have same degree, the choice is arbitrary
	"""
	def __init__(self, debug=0):
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
		"""
		if self.debug:
			print "no of singleton vertices to remove", len(vertices_to_remove)
		"""
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
		while degree_percentage<degree_cut_off:
			"""
			if self.debug:
				print "vertex %s removed with min_degree %s, degree_percentage %s"%(graph.vertex_properties['vertex_id'][vertex_min_degree], min_degree, degree_percentage)
			"""
			graph.clear_vertex(vertex_min_degree)
			graph.remove_vertex(vertex_min_degree)
			vertex_min_degree, min_degree = self.get_vertex_min_degree(graph)
			degree_percentage = float(min_degree)/(graph.num_vertices()-1)
		return degree_percentage
	
	
	def get_core_vertex_set(self, edge2encodedOccurrence, vertex_list, recurrence_array, degree_cut_off):
		"""
		12-16-05
			--init_graph_from_vertex_set()
			--decodeOccurrence()
			--remove_loose_part_of_graph()
				--remove_singleton_vertices()
				--get_vertex_min_degree()
		"""
		no_of_datasets = len(recurrence_array)
		graph_list = [None]*no_of_datasets
		anti_vertex_id_list = [None]*no_of_datasets
		for i in range(no_of_datasets):
			if recurrence_array[i] == 1:
				graph_list[i], anti_vertex_id_list[i] = self.init_graph_from_vertex_set(vertex_list)
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
				if edge_tuple in edge2encodedOccurrence:
					edge_recurrence = decodeOccurrence(edge2encodedOccurrence[edge_tuple])	#starting from 1
					"""
					if self.debug:
						print "edge_recurrence", edge_recurrence
					"""
					for k in edge_recurrence:
						index = k-1
						if recurrence_array[index] == 1:	#just those on datasets
							"""
							if self.debug:
								print "checking dataset", index
							"""
							v_descriptor1 = anti_vertex_id_list[index][vertex_list[i]]
							v_descriptor2 = anti_vertex_id_list[index][vertex_list[j]]
							graph_list[index].add_edge(v_descriptor1, v_descriptor2)
		#remove loose part for each graph
		vertex_id2occurrence = {}
		on_dataset_index_ls = []
		for i in range(no_of_datasets):
			if graph_list[i]:
				"""
				if self.debug:
					print "checking graph", i
				"""
				degree_percentage = self.remove_loose_part_of_graph(graph_list[i], degree_cut_off)
				if graph_list[i].num_vertices()>=4:	#min graph size
					"""
					if self.debug:
						print "graph %s has %s vertices remaining with degree_percentage: %s."%(i, graph_list[i].num_vertices(), degree_percentage)
					"""
					on_dataset_index_ls.append(i)	#this dataset should be counted as 'on'
					for v in graph_list[i].vertices:
						vertex_id = graph_list[i].vertex_properties['vertex_id'][v]
						if vertex_id not in vertex_id2occurrence:
							vertex_id2occurrence[vertex_id]	= 0
						vertex_id2occurrence[vertex_id] += 1
		core_vertex_set = []
		for vertex_id in vertex_id2occurrence:
			if vertex_id2occurrence[vertex_id] == len(on_dataset_index_ls):
				core_vertex_set.append(vertex_id)
		core_vertex_set.sort()
		return core_vertex_set, on_dataset_index_ls


if __name__ == '__main__':
	"""
	12-16-05 give it a test
	"""
	import sys,os
	sys.path += [os.path.expanduser('~/script/annot/bin')]
	from codense.common import db_connect
	from MpiCrackSplat import MpiCrackSplat
	hostname='zhoudb'
	dbname='graphdb'
	schema = 'scfim30'
	conn, curs = db_connect(hostname, dbname, schema)
	MpiCrackSplat_instance = MpiCrackSplat()
	edge2encodedOccurrence = {}
	min_sup = 3
	max_sup = 30
	no_of_datasets = MpiCrackSplat_instance.fill_edge2encodedOccurrence(hostname, dbname, \
	schema, edge2encodedOccurrence, min_sup, max_sup)
	vertex_list = [854078, 854317, 854381, 854469, 855539, 855547, 855613, 855730, 855862, 856019, 856226]
	recurrence_array =  [0.0, 0.0, 0.0, 1.0, 0.53333333333333333, 0.0, 0.066666666666666666, \
	0.20000000000000001, 0.33333333333333331, 0.0, 0.26666666666666666, 0.20000000000000001, \
	0.0, 0.0, 0.33333333333333331, 0.0, 0.66666666666666663, 1.0, 0.0, 0.66666666666666663, 1.0, \
	0.0, 0.066666666666666666, 0.0, 0.46666666666666667, 0.0, 0.13333333333333333, \
	0.20000000000000001, 0.26666666666666666, 1.0]
	degree_cut_off = 0.3
	fuzzyDense_instance = fuzzyDense(debug=1)
	print fuzzyDense_instance.get_core_vertex_set(edge2encodedOccurrence, vertex_list, recurrence_array, degree_cut_off)
