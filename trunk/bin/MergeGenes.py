#!/usr/bin/env python
"""
Usage: MergeGenes.py [OPTION] -o OUTPUTDIR FILES

Option:
	FILES are a list of files whose columns are going to be counted.
	-o ... --outputdir=...,	directory containing output files
	-d ... --delimiter=...,	delimiter character used to seperate columns, \t(default)
	-b, --debug	enable debugging
	-h, --help              show this help
	
Examples:
	MergeGenes.py -o datasets/sc_log/ datasets/sc/*

Description:
	Merge the expression values of same genes.
	Inputfiles will be sorted before transforming by calling 'sort'.
"""

import sys, os, re, getopt, csv, math
import MLab
sys.path += [os.path.expanduser('~/script/microarray/bin')]
from microarraydb import microarraydb
from MA import array, average, maximum
from Preprocess import PreprocessEdgeData
from sets import Set

class MergeGenes:
	def __init__(self, file_list, outputdir, delimiter, debug=0):
		""" 
		08-28-05
		"""
		self.files = file_list
		self.files.sort()
		self.outputdir = outputdir
		self.delimiter = delimiter
		self.debug = int(debug)

	def mergeBlock(self, data_ls_2d, mask_ls_2d, gene_id_set, bad_gene_id):
		"""
		08-28-05
			merge the block of same genes
		"""
		gene_id_set.remove(bad_gene_id)
		gene_id = gene_id_set.pop()
		ar = array(data_ls_2d, mask = mask_ls_2d)
		if len(data_ls_2d)==1:	#no need to do average
			ar = array(data_ls_2d[0], mask=mask_ls_2d[0])
			return gene_id, ar
		max_ls = []
		for i in range(ar.shape[0]):
			max_value = maximum(ar[i,:])
			max_ls.append(max_value)
			ar[i,:] = ar[i,:]/max_value
		new_ar = average(ar)*max(max_ls)
		return gene_id, new_ar
		
	def transform_one_file(self, src_pathname, delimiter, outputdir, PreprocessEdgeData_instance,\
		microarraydb_instance):
		"""
		08-28-05
			the src_pathname is sorted in advance
		10-26-05 sort the input file first
		"""
		microarraydb_instance.sort_outputfile(src_pathname)	#10-26-05
		reader = csv.reader(file(src_pathname), delimiter=delimiter)
		filename = os.path.basename(src_pathname)
		output_filename = os.path.join(outputdir, filename)
		writer = csv.writer(open(output_filename, 'w'), delimiter=delimiter)
		gene_id_set = Set()	#check whether the same gene block is passed
		data_ls_2d = []
		mask_ls_2d = []
		for row in reader:
			gene_id = row[0]
			data_ls = []
			mask_ls = []
			for i in range(1, len(row)):
				if row[i] == 'NA':
					data_ls.append(1e20)
					mask_ls.append(1)
				elif row[i] == '':
					#ignore empty entry
					continue
				else:
					value = float(row[i])
					data_ls.append(value)
					mask_ls.append(0)
			if sum(mask_ls) == len(mask_ls):	#all NAs
				continue
			gene_id_set.add(gene_id)
			if self.debug:
				print "gene_id", gene_id
				print "data_ls", data_ls
				print "mask_ls", mask_ls
			if len(gene_id_set)==1:
				data_ls_2d.append(data_ls)
				mask_ls_2d.append(mask_ls)
			elif len(gene_id_set)==2:
				if self.debug:
					print 'data_ls_2d is ',data_ls_2d
					print 'mask_ls_2d is ', mask_ls_2d
				pre_gene_id, ar = self.mergeBlock(data_ls_2d, mask_ls_2d, gene_id_set, gene_id)
				if self.debug:
					print 'final block is ', ar
					raw_input("Continue?(Y/n)")
				if ar.mask():
					writer.writerow([pre_gene_id] + PreprocessEdgeData_instance.ls_NA_fillin(ar))
				else:	#no mask, just list
					writer.writerow([pre_gene_id] + list(ar))
				#clean up data structures
				gene_id_set.clear()
				gene_id_set.add(gene_id)
				data_ls_2d = [data_ls]
				mask_ls_2d = [mask_ls]
			else:
				sys.stderr.write("Error: length of gene_id_set %s is %s.(not 1 or 2)\n"%(repr(gene_id_set), len(gene_id_set)))
				sys.exit(2)
		del reader, writer
	
	def run(self):
		"""
		08-28-05
		"""
		PreprocessEdgeData_instance = PreprocessEdgeData()
		microarraydb_instance = microarraydb()
		if not os.path.isdir(self.outputdir):
			os.makedirs(self.outputdir)
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(self.files))
		for f in self.files:
			sys.stderr.write("%d/%d:\t%s"%(self.files.index(f)+1,len(self.files),f))
			self.transform_one_file(f, self.delimiter, self.outputdir, PreprocessEdgeData_instance, microarraydb_instance)
			sys.stderr.write(".\n")


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "ho:d:b", ["help", "outputdir=", "delimiter=", "debug"])
	except:
		print __doc__
		sys.exit(2)

	delimiter = '\t'
	outputdir = None
	debug = 0
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-o", "--outputdir"):
			outputdir = arg
		elif opt in ("-d", "--delimiter"):
			delimiter = arg
		elif opt in ("-b", "--debug"):
			debug = 1
			
	if len(args)>=1:
		instance = MergeGenes(args, outputdir, delimiter, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
