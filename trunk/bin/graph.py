#!/usr/bin/env python

import sys
import math, os

class graph_construct:
	def __init__(self, arg1, debug=0):
		import rpy
		import numarray.ma as ma
		import numarray
		self.raw_dataset_source = arg1
		self.dataset_source = os.tempnam('/tmp','annot')
		self.gene_cut_off = 8
		self.jk_cor_cut_off = 7
		self.cor_cut_off =0.6
		self.debug = debug
		self.graph_dict = {}		#a dictionary to store edges
		self.genelabels = []		#a list to store gene labels
		self.cor_vector = []		#a vector to store all the correlations between two vectors.
		
		#self.no_of_genes, #rows
		#self.no_of_cols
		#self.mask_array, the data_struc for the real microarray data
		#self.mask_matrix, used to leave-one-out in correlation calculation
		
		#the variables commented will be initiated later.
		self.preprocess()

	def preprocess(self):
		fd = open(self.raw_dataset_source, 'r')
		fd_new = open(self.dataset_source, 'w')
		line = fd.readline()
		self.no_of_cols = len(line.split())-1
		while line:
			if line[0] == '#' or line =='\n':
				pass
			elif line[0] == '\t':
				self.no_of_cols = len( line.split())
				fd_new.write(line)
			else:
				list = line.split()
				if self.no_of_cols - list.count('NA') <self.gene_cut_off:
					pass
				else:
					fd_new.write(line)
			line=fd.readline()
		fd.close()
		fd_new.close()

		
	def mask_array_construct(self):
		data = rpy.with_mode(0, rpy.r.read_table)(self.dataset_source, row_names=1)
		'''
		!Important!
		if the dataset_source has too few data, conversion from R to python will be a problem.
		The whole data matrix will be converted to a python string matrix.
		R's NA is not converted to nan in python.
		
		The problem has been found. 
		r.as_matrix converts small dataset to character type.
		r.matrix won't rig the class type, but it rigs the structure.
		The only to sovle this is add a colClasses vector to r.read_table.
		such as: colClasses=c('character',rep('double',11))
		But you have to know the no_of_cols in advance.
		
		As our dataset is really big, this problem hasn't appeared.
		
		'''
		#print r.as_matrix(data)
		array = ma.masked_inside(rpy.r.as_matrix(data),  -1.0e20, 1.0e20)
		#all are set to be masked except nan. weird! So have to do a converse.
		self.mask_array = ma.array(array, mask=ma.logical_not(ma.getmask(array)))
		self.genelabels = rpy.r.rownames(data)
		self.no_of_genes = len(self.genelabels)
		self.no_of_cols = len(array[0])
		self.mask_matrix=ma.identity(self.no_of_cols)
		del array ,data
		
	def edge_construct(self):
		for i in range(self.no_of_genes):
			#after preprocessing, theses filters are of no use.
			'''
			if self.no_of_cols - ma.sum(self.mask_array[i].mask()) <self.gene_cut_off:
				if self.debug:
					sys.stderr.write( 'jump_out level 0\t' + self.genelabels[i])
				continue
				#less than 8 valid data spots
			'''
			for j in range(i+1, self.no_of_genes):
				'''
				if self.no_of_cols - ma.sum(self.mask_array[j].mask()) <self.gene_cut_off:
					if self.debug:
						sys.stderr.write(print 'jump_out level 1\t' + self.genelabels[j])
					continue
					#less than 8 valid data spots
				'''
				joint_mask = ma.mask_or(self.mask_array[i].mask(), self.mask_array[j].mask())		#joint mask
				self.cor_vector = []			#initilization
				nn_cor_vector = [] 			#non-negative version of co_vector
				for k in range(self.no_of_cols):
					new_mask = ma.mask_or(joint_mask, self.mask_matrix[k])		#leave k out
					if self.no_of_cols - ma.sum(new_mask) < self.jk_cor_cut_off:
						#if self.debug:
						#	sys.stderr.write( 'jump_out level 2\t%s v.s %s at %d\n'%(self.genelabels[i], self.genelabels[j], k,))
						continue
						#less than jk_cor_cut_off, no correlation
					v1 = ma.array(self.mask_array[i], mask=new_mask).compressed().tolist()
					v2 = ma.array(self.mask_array[j], mask=new_mask).compressed().tolist()
					cor = rpy.r.cor(v1,v2)
					self.cor_vector.append( cor)
					nn_cor_vector.append(math.fabs(cor))
					
					if self.no_of_cols-ma.sum(joint_mask) == self.jk_cor_cut_off:
						break
					#Only jk_cor_cut_off(7) valid quantities shared by two genes. 
					#All the leave-one-out cor's are same. You can only leave NA out.
					
				if len(self.cor_vector) >0:
					min_cor = min(nn_cor_vector)		#minimum in the non-negative version of cor_vector
					if min_cor >= self.cor_cut_off:
						if self.debug:
							sys.stderr.write('cor vector of %s v.s. %s: %s\n'%(self.genelabels[i], self.genelabels[j],self.cor_vector,))
						self.graph_dict[(self.genelabels[i],self.genelabels[j])] = self.cor_vector[nn_cor_vector.index(min_cor)]
						#not simply the smallest, but the smallest absolute value
	
	def cleanup(self):
		os.remove(self.dataset_source)
		del self.mask_array
		
	def output(self, out=sys.stdout):
		out.write( 't\t#\t%s\n'%os.path.basename(self.raw_dataset_source))
		for i in  self.graph_dict:
			out.write( 'e\t%s\t%s\t%f\n'%(i[0], i[1] , self.graph_dict[i],))

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
				
class graph_reorganize:
	'''
	Read the output from graph_construct class(above), convert it to gSpan or Splat input format.
	Store the mapping of vertex_label v.s. number and graph_label v.s. number in file no_label_map in user's home directory.
	Based on this class, optimization of all the graphs is possible.
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
				no_of_graphs = len(self.graph_list_dict)
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
					no_of_vertices = len(self.vertex_list_dict)
					self.vertex_list_dict[vertex] = vertex_attributes(no=no_of_vertices)
					self.vertex_block+='v %d %d\n'%(no_of_vertices,no_of_vertices,)
		self.output()
		
	def output(self):
		import os
		no_label_map_filename = os.path.join(os.path.expanduser('~'),'no_label_map')
		map_fhandler = open(no_label_map_filename, 'w')
		map_fhandler.write('>no_graphlabel_mapping\n')
		for graph_label in self.graph_list_dict:
			outf = open('splan_'+graph_label, 'w')
			graph_attr = self.graph_list_dict[graph_label]
			outf.write('t # %d\n'%graph_attr.no)
			
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
			
if __name__ == '__main__':
	'''
	instance = graph_construct(sys.argv[1])
	instance.mask_array_construct()
	instance.edge_construct()
	instance.cleanup()
	instance.output()
	'''
	instance = graph_reorganize()
	inf = open(sys.argv[1], 'r')
	instance.parse(inf)
