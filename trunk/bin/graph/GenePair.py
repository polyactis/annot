#!/usr/bin/env python
"""
The backend for the GenePairGui_*.py programs.

"""
import sys, os, getopt, csv
import graph_modeling
from rpy import r
from array import array    # need arrays to pass to ROOT
#from ROOT import gROOT, TCanvas, TPad, TH1F, TImage, TGraph, TMultiGraph

class GenePair:
	'''

	'''
	def __init__(self, input_file):
		self.input_file = input_file
		#the dictionary storing all the expression information
		self.gene_id2expr_array = {}
		
		#the gene_pair calculation dictionary
		self.pairwise_calculate = {'cor': graph_modeling.ind_min_cor}
		#data structure loadin
		self.dstruc_loadin()
		#the file to contain the image
		self.plot_file = '/tmp/genepair.png'

	def dstruc_loadin(self):
		'''
		This method loads in the data structures from database.
		'''
		sys.stderr.write("Loading Data STructure...")
		self.gene_id2expr_array_setup(self.input_file)
		sys.stderr.write("Done\n")

	def gene_id2expr_array_setup(self, f_path):
		reader = csv.reader(file(f_path), delimiter='\t')
		for row in reader:
			gene_id = row[0]
			if gene_id in self.gene_id2expr_array:
				sys.stderr.write("Duplicate Genes: %s\n"%gene_id)
				sys.exit(2)
			else:
				for i in range(1,len(row)):
					if row[i] == 'NA':
						#100000000 is regarded as NAN
						row[i] = float(100000000)
					elif row[i]!='':
						row[i] = float(row[i])
				if row[-1] == '':
					#sometimes the dataset has a last empty column
					new_row = row[1:-1]
				else:
					new_row = row[1:]
				self.gene_id2expr_array[gene_id] = new_row
		del reader
	
	def NA_cleanup(self, vector1, vector2):
		new_vector1 = []
		new_vector2 = []
		for i in range(len(vector1)):
			if vector1[i] == 100000000 or vector2[i] == 100000000:
				continue
			else:
				new_vector1.append(vector1[i])
				new_vector2.append(vector2[i])
		return (new_vector1, new_vector2)

	def xy_list_return(self, vector):
		'''
		return x_list and y_list for plotting, discard those NA values(1e8).
		'''
		x_list = []
		y_list = []
		for i in range(len(vector)):
			if vector[i] != 1e8:
				x_list.append(i+1)
				y_list.append(vector[i])
		return (x_list, y_list)
	
	def get_min_max(self, vector_list):
		'''
		get the minimum and maximum out of vector_list
		'''
		whole_list = []
		for vector in vector_list:
			whole_list += vector
		whole_list.sort()
		min = whole_list[0]
		for i in range(1,len(whole_list)+1):
			if whole_list[-i]>1e8:
				sys.stderr.write("Error, maximum gene expr value > 1e8, NAN\n")
			elif whole_list[-i] == 1e8:
				#it's NAN
				continue
			else:
				#got it
				max = whole_list[-i]
				break

		return (min, max)
	
	def plot(self, vector_list, gene_id_list):
		self.no_of_curves = 0
		x_range = (1, len(vector_list[0]))
		y_range = self.get_min_max(vector_list)
		r.png("%s"%self.plot_file)
		for vector in vector_list:
			(x_list, y_list) = self.xy_list_return(vector)
			self._plot(x_list, y_list, x_range, y_range)
		
		r.legend(x_range[1], y_range[1], gene_id_list, col=range(1, self.no_of_curves+1), lty=1, pch='*', xjust=1)
		r.dev_off()
	
	def _plot(self, x_list, y_list, x_range, y_range):
		self.no_of_curves += 1
		if self.no_of_curves==1:
			r.plot(x_list, y_list, type='o',pch='*',xlab='dataset no.',xlim=x_range,ylim=y_range, \
			ylab='expression value', col=self.no_of_curves)
		else:
			r.lines(x_list, y_list, type='o',pch='*',col=self.no_of_curves)

	
	"""
	#function deprecated
	def root_draw(self, vector1, vector2, gene_id1, gene_id2):
		# create arrays to pass to ROOT later
		array0 = array('d')
		array1 = array('d')
		array2 = array('d')
		#starting from 1
		array0.fromlist(range(1, len(vector1)+1))
		array1.fromlist(vector1)
		array2.fromlist(vector2)
		gROOT.Reset()
		c1 = TCanvas("c1", "Plot", 800,700)
		c1.SetGridx()
		c1.SetGridy()
		gr1 = TGraph(len(vector1), array0, array1)
		gr2 = TGraph(len(vector1), array0, array2)
		mg = TMultiGraph()
		mg.Add(gr1)
		mg.Add(gr2)
		mg.Draw("AL*")
		stop = raw_input("Exit: Y/n")
	"""
	
	def gene_pair_analyze(self, gene_id_list):
		vector_list = []
		#gene_id_list may contain some inexistent genes
		real_gene_id_list = []
		for gene_id in gene_id_list:
			if gene_id in self.gene_id2expr_array:
				real_gene_id_list.append(gene_id)
				vector_list.append(self.gene_id2expr_array[gene_id])	
			else:
				sys.stderr.write("%s doesn't appear in the dataset\n"%(gene_id))
		#(new_vector1, new_vector2) = self.NA_cleanup(vector1, vector2)
		for key in self.pairwise_calculate:
			for i in range(len(vector_list)):
				for j in range(i+1, len(vector_list)):
					edge_data = self.pairwise_calculate[key](vector_list[i], vector_list[j])
					sys.stdout.write("%s %s %s: %s\n"%(real_gene_id_list[i], real_gene_id_list[j], key, edge_data.value))
		if len(real_gene_id_list)>0:
			self.plot(vector_list, real_gene_id_list)
			return self.plot_file
		else:
			return None

	
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	else:
		instance = GenePair(sys.argv[1])
		instance.gene_pair_analyze("YKR006C", "YKR007W")
