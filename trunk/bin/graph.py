#!/usr/bin/env python

import sys
from rpy import *
import numarray.ma as ma
import numarray
import math, os

class graph_construct:
	def __init__(self, arg1, debug=0):
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
		data =with_mode(0, r.read_table)(self.dataset_source)
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
		array = ma.masked_inside(r.as_matrix(data),  -1.0e20, 1.0e20)
		#all are set to be masked except nan. weird! So have to do a converse.
		self.mask_array = ma.array(array, mask=ma.logical_not(ma.getmask(array)))
		self.genelabels = r.rownames(data)
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
					cor = r.cor(v1,v2)
					self.cor_vector.append( cor)
					nn_cor_vector.append(math.fabs(cor))
					
					if self.no_of_cols-ma.sum(joint_mask) == self.jk_cor_cut_off:
						break
					#Only jk_cor_cut_off(7) valid quantities shared by two genes. 
					#All the leave-one-out cor's are same. You can only leave NA out.
					
				if len(self.cor_vector) >0:
					min_cor = min(nn_cor_vector)		#minimum in the non-negative version of cor_vector
					if min_cor >= self.cor_cut_off:
						#if self.debug:
						#	sys.stderr.write('cor vector of %s v.s. %s: %s\n'%(self.genelabels[i], self.genelabels[j],self.cor_vector,))
						self.graph_dict[(self.genelabels[i],self.genelabels[j])] = self.cor_vector[nn_cor_vector.index(min_cor)]
						#not simply the smallest, but the smallest absolute value
	
	def cleanup(self):
		os.remove(	self.dataset_source)
		del self.mask_array
		
	def output(self, out=sys.stdout):
		out.write( 't\t#\t%s\n'%os.path.basename(self.raw_dataset_source))
		for i in  self.graph_dict:
			out.write( 'e\t%s\t%s\t%f\n'%(i[0], i[1] , self.graph_dict[i],))
		
if __name__ == '__main__':
	instance = graph_construct(sys.argv[1])
	instance.mask_array_construct()
	instance.edge_construct()
	instance.cleanup()
	instance.output()
