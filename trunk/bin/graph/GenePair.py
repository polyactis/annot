#!/usr/bin/env python
"""
The backend for the GenePairGui_*.py programs.

"""
import sys, os, getopt, csv
import graph_modeling
from array import array    # need arrays to pass to ROOT
from ROOT import gROOT, TCanvas, TPad, TH1F, TImage, TGraph, TMultiGraph

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
		
	def gene_pair_analyze(self, gene_id1, gene_id2):
		if gene_id1 in self.gene_id2expr_array and gene_id2 in self.gene_id2expr_array:
			vector1 = self.gene_id2expr_array[gene_id1]
			vector2 = self.gene_id2expr_array[gene_id2]
		else:
			sys.stderr.write("%s and %s don't appear in the dataset\n"%(gene_id1, gene_id2))
			return
		(new_vector1, new_vector2) = self.NA_cleanup(vector1, vector2)
		if len(new_vector1) <3:
			sys.stderr.write("Only %s valid columns.\n"%len(new_vector1))
		else:
			for key in self.pairwise_calculate:
				edge_data = self.pairwise_calculate[key](new_vector1, new_vector2)
				sys.stdout.write("%s: %s\n"%(key, edge_data.value))
			self.root_draw(new_vector1, new_vector2, gene_id1, gene_id2)

	
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	else:
		instance = GenePair(sys.argv[1])
		instance.gene_pair_analyze("YKR006C", "YKR007W")
