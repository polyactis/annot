#!/usr/bin/env python

import sys, pickle,os

class vertex_attributes:
	"a class wrapping the attributes of a vertex"
	def __init__(self, no = None):
		self.freq = 1
		self.no = no
	
	def inc_freq(self):
		self.freq += 1

class graph_attributes:
	"a class wrapping the attributes of a graph"
	def __init__(self, no=None):
		self.no = no
		self.graph_dict = {}
		self.vertex_set = {} #key stores vertex's label, value is a list which stores the adjacent vertices.
	def vertex_set_init(self):
		for i in self.graph_dict:
			if self.vertex_set.has_key(i[0]):
				self.vertex_set[i[0]].append(i[1])
			if i[0] not in self.vertex_set:
				self.vertex_set[i[0]] = []
			if self.vertex_set.has_key(i[1]):
				self.vertex_set[i[1]].append(i[0])
			if i[1] not in self.vertex_set:
				self.vertex_set[i[1]] = []
				
class graph_reorganize_old:
	'''
	Read the output from graph_construct class(above), convert it to gSpan or Splat input format.
	Store the mapping of vertex_label v.s. number and graph_label v.s. number in file no_label_map in user's home directory.
	Based on this class, optimization of all the graphs is possible.
	
	!!!Defuncted. This class requires too much memory.
	Replaced by the class graph_reorganize.
	'''
	def __init__(self):
		self.vertex_list_dict = {}
		self.graph_list_dict = {}
		self.vertex_block=''		#a block of vertices in the output file. "v M L" means the Mth vertex with L label.
		
	def init(self):
		self.vertex_list_dict = {}
		self.graph_list_dict = {}
		self.vertex_block = ''
		
	def parse(self, inf):
		self.init()

		#loop below initilizes graph_list_dict.
		line = inf.readline()
		while line:
			list = line[:-1].split('\t')
			if line[0] == 't':
				graph_label = list[2]
				no_of_graphs = len(self.graph_list_dict) + 1
				self.graph_list_dict[graph_label] = graph_attributes(no=no_of_graphs)
			if line[0] == 'e':
				vertex1 = list[1]
				vertex2 = list[2]
				self.graph_list_dict[graph_label].graph_dict[(vertex1,vertex2)] = float(list[3])
			line = inf.readline()
		
		#loop below initilizes vertex_list_dict.
		for graph_label in self.graph_list_dict:
			graph_attr = self.graph_list_dict[graph_label]
			graph_attr.vertex_set_init()
			#initlize vertex_set of each graph

			for vertex in graph_attr.vertex_set:
				if self.vertex_list_dict.has_key(vertex):
					self.vertex_list_dict[vertex].inc_freq()
				else:
					no_of_vertices = len(self.vertex_list_dict) + 1
					self.vertex_list_dict[vertex] = vertex_attributes(no=no_of_vertices)
					self.vertex_block+='v %d %s\n'%(no_of_vertices,vertex,)
		self.output()
		
	def output(self):
		import os
		no_label_map_filename = os.path.join(os.path.expanduser('~'),'no_label_map')
		map_fhandler = open(no_label_map_filename, 'w')
		map_fhandler.write('>no_graphlabel_mapping\n')
		for graph_label in self.graph_list_dict:
			outf = open('splat_'+graph_label, 'w')
			graph_attr = self.graph_list_dict[graph_label]
			outf.write('t # %s\n'%graph_label)
			
			map_fhandler.write('%d\t%s\n'%(graph_attr.no, graph_label,))
			
			outf.write(self.vertex_block)

			for edge in graph_attr.graph_dict:
				vertex1_no = self.vertex_list_dict[edge[0]].no
				vertex2_no = self.vertex_list_dict[edge[1]].no
				outf.write('e %d %d %f\n'%(vertex1_no, vertex2_no, graph_attr.graph_dict[edge],))
			outf.close()

		map_fhandler.write('>no_vertexlabel_mapping\n')
		for vertex_label in self.vertex_list_dict:
			vertex_attr = self.vertex_list_dict[vertex_label]
			map_fhandler.write('%d\t%s\n'%(vertex_attr.no, vertex_label,))
			
class graph_reorganize:
	'''
	This class plays two major roles.
	Collect labels.
	Transform graph.cc results into gspan input.
	'''
	def __init__(self):
		self.global_vertex_list = []
		self.global_graph_list = []
		self.global_vertex_dict = {}
		self.global_graph_dict = {}
		self.pickle_fname = os.path.join(os.path.expanduser('~'), 'pickle/yeast_global_struc')
		self.vertex_block = ''
		
	def global_mapping_construct(self, inf):
		'''
		global_mapping_construct() collects the graph label and vertex label from inf.
		mapping_batch() will iteratively call this function to collect all labels pertaining
		to one organism.
		The collection will be stored in a file through pickle.
		'''
		line = inf.readline()
		while line:
			list = line[:-1].split()
			if line[0] == 't':
				graph_label = list[2]
				if graph_label in self.global_graph_dict:
					sys.stderr.write('Error. graph with this name, "%s" appears twice\n'%graph_label)
				self.global_graph_dict[graph_label]=1
			if line[0] == 'e':
				vertex1_label = list[1]
				vertex2_label = list[2]
				if vertex1_label not in self.global_vertex_dict:
					self.global_vertex_dict[vertex1_label] = 1
				if vertex2_label not in self.global_vertex_dict:
					self.global_vertex_dict[vertex2_label] = 1
			line = inf.readline()

	def dstruc_loadin(self):
		'''
		This method loads in the data structures from a pre-stored file through pickle.
		'''
		if not os.path.isfile(self.pickle_fname):
			sys.stderr.write('Error, file: %s not existent\n'%self.pickle_fname)
			return	1
		global_struc = pickle.load(open(self.pickle_fname, 'r'))
		self.global_graph_list = global_struc['graph_list']
		self.global_vertex_dict = global_struc['vertex_dict']
		self.global_vertex_list = global_struc['vertex_list']
		for i in self.global_vertex_list:
			self.vertex_block += 'v %d %s\n'%(self.global_vertex_list.index(i)+1, i)
		return 0

	def transform(self, inf, outf):
		'''
		This method transforms a graph.cc output file(inf) into a gspan input file(outf).
		It requires self.global_vertex_list and other data structures to be filled in first.
		So self.dstruc_loadin() must be called before this method.
		transform_batch() will iteratively call this method to transform all graph.cc output files.
		'''
		line = inf.readline()
		while line:
			list = line[:-1].split()
			if line[0] == 't':
				outf.write(line.replace('\t', ' '))
				outf.write(self.vertex_block)
			if line[0] == 'e':
				vertex1_label = list[1]
				vertex2_label = list[2]
				outf.write('e %d %d %s\n'% \
								(self.global_vertex_dict[vertex1_label], \
								self.global_vertex_dict[vertex2_label], \
								list[3],)
								)
			line = inf.readline()

def transform_batch(dir, output_dir):
	'''
	See comments of graph_reorganize.transform().
	'''
	files = os.listdir(dir)
	sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
	if not os.path.isdir(output_dir):
		os.makedirs(output_dir)
	
	instance = graph_reorganize()
	data_not_loaded = instance.dstruc_loadin()
	if data_not_loaded:
		return
	
	for f in files:
		pathname = os.path.join(dir, f)
		sys.stderr.write("%d/%d:\t%s\n"%(files.index(f)+1,len(files),f))
		inf = open(pathname, 'r')
		outf = open(os.path.join(output_dir,f), 'w')
		instance.transform(inf, outf)
		inf.close()
		outf.close()

def mapping_batch(dir):
	'''
	See comments of graph_reorganize.global_mapping_construct().
	'''
	files = os.listdir(dir)
	sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
	instance = graph_reorganize()	
	for f in files:
		pathname = os.path.join(dir, f)
		sys.stderr.write("%d/%d:\t%s\n"%(files.index(f)+1,len(files),f))
		inf = open(pathname, 'r')
		instance.global_mapping_construct(inf)
		inf.close()

	instance.global_vertex_list = instance.global_vertex_dict.keys()
	instance.global_vertex_list.sort()
	for i in range(len(instance.global_vertex_list)):
		instance.global_vertex_dict[instance.global_vertex_list[i]] = i+1
	global_struc = {'vertex_dict': instance.global_vertex_dict,
			'vertex_list': instance.global_vertex_list,
			'graph_list': instance.global_graph_dict}
	pickle.dump(global_struc, open(instance.pickle_fname, 'w') )
	

if __name__ == '__main__':
	'''
	# this block uses the old class.
	instance = graph_reorganize_old()
	inf = open(sys.argv[1], 'r')
	instance.parse(inf)
	'''
	if sys.argv[1] == 'mapping':
		mapping_batch(sys.argv[2])
		#argv[2] specifies the directory which contains the graph.cc output files
	elif sys.argv[1] == 'transform':
		transform_batch(sys.argv[2], sys.argv[3])
		#argv[2] specifies the directory which contains the graph.cc output files.
		#argv[3] specifies the directory to store the gspan input files.
