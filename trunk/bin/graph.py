#!/usr/bin/env python

import sys
from rpy import *
import numarray.ma as ma
import numarray

class graph_construct:
	def __init__(self, arg1, debug=0):
		self.dataset_source = arg1
		self.debug = debug
		self.graph_dict = {}
		self.genelabels = []
		self.cor_vector = []		#a vector to store all the correlations between two vectors.
		
		#self.no_of_genes, #rows
		#self.no_of_cols
		#self.mask_array, the data_struc for the real microarray data
		#self.mask_matrix, used to leave-one-out in correlation calculation
		
		#the variables commented will be initiated later.

		
	def mask_array_construct(self):
		data =with_mode(0, r.read_table)(self.dataset_source)
		'''
		!Important!
		if the dataset_source has too few data, conversion from R to python will be a problem.
		The whole data matrix will be converted to a python string matrix.
		R's NA is not converted to nan in python.
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
			if self.no_of_cols - ma.sum(self.mask_array[i].mask()) <8:
				if self.debug:
					print 'jump_out level 0\t' + self.genelabels[i]
				continue
				#less than 8 valid data spots
			for j in range(i+1, self.no_of_genes):
				if self.no_of_cols - ma.sum(self.mask_array[j].mask()) <8:
					if self.debug:
						print 'jump_out level 1\t' + self.genelabels[j]
					continue
					#less than 8 valid data spots
				mask_tmp = ma.mask_or(self.mask_array[i].mask(), self.mask_array[j].mask())		#joint mask
				self.cor_vector = []			#initiliation
				for k in range(self.no_of_cols):
					new_mask = ma.mask_or(mask_tmp, self.mask_matrix[k])		#leave k out
					if self.no_of_cols - ma.sum(new_mask) <7:
						if self.debug:
							print 'jump_out level 2\t%s v.s %s at %d'%(self.genelabels[i], self.genelabels[j], k,)
						continue
						#less than 6, no correlation
					v1 = ma.array(self.mask_array[i], mask=new_mask).compressed().tolist()
					v2 = ma.array(self.mask_array[j], mask=new_mask).compressed().tolist()
					self.cor_vector.append( r.cor(v1,v2))
				if len(self.cor_vector) >0:
					min_cor = min(self.cor_vector)
					if min_cor >= 0.6:
						if self.debug:
							print 'cor vector of %s v.s. %s: %s'%(self.genelabels[i], self.genelabels[j],self.cor_vector,)
						self.graph_dict[(self.genelabels[i],self.genelabels[j])] = min_cor
						
					
if __name__ == '__main__':
	instance = graph_construct(sys.argv[1], debug=1)
	instance.mask_array_construct()
	instance.edge_construct()
	print instance.graph_dict
