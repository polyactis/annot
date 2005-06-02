#!/usr/bin/env python
"""
Usage: PreprocessEdgeData.py -k SCHEMA [OPTION] R_FILE

Option:
	R_FILE is the file to store the R code.
	-i ..., --input_file=...	gspan format
	-o ..., --output_file=...	haiyan's matrix format
	-n ..., --no_of_nas=...	minimum number of NAs of one edge, None(default)
	-t ..., --top_percentage=...	0.1(default)
	-b, --debug	just fetch 15000 edges and do a demo
	-h, --help              show this help

Examples:
	PreprocessEdgeData.py -i edge_data_go_282 -o edge_data_go_282.1
	
Description:

	
"""

import sys, os, getopt, csv, re
import MLab
from Numeric import argsort, take, transpose
from MA import array

class PreprocessEdgeData:
	"""
	05-09-05
		a module to preprocess the edge data of a function category, which
		is the correlation vector of an edge across all datasets.
		1. only take the 10% edges with highest std
		2. NA <= 7
		3. transpose the matrix
	"""
	def __init__(self, infname=None, outfname=None, no_of_nas=None, top_percentage=0.1, debug=0):
		self.infname = infname
		self.outfname = outfname
		self.no_of_nas = no_of_nas
		self.top_percentage = top_percentage
		self.debug = int(debug)
		
		self.list_of_mas = []
		self.list_of_top_mas = []
		
		
	def data_read_in(self, infname, no_of_nas):
		"""
		05-09-05
		"""
		sys.stderr.write("Reading data...")
		list_of_mas = []
		reader = csv.reader(open(infname, 'r'),delimiter='\t')
		reader.next()	#ignore the first line
		for row in reader:
			data_ls = []
			mask_ls = []
			for item in row[1:]:	#ignore the first edge id
				if item=='NA':
					data_ls.append(1.1)
					mask_ls.append(1)
				else:
					data_ls.append(float(item))
					mask_ls.append(0)
			if no_of_nas:
				if sum(mask_ls)>=no_of_nas:	#too many NAs
					continue
			list_of_mas.append(array(data_ls, mask=mask_ls))
		del reader
		sys.stderr.write("Done.\n")
		return list_of_mas
	
	def get_top_mas(self, list_of_mas, top_percentage):
		"""
		05-09-05
		"""
		sys.stderr.write("Getting the top %s std edges..."%top_percentage)
		list_of_stds = []
		for ma in list_of_mas:
			std = MLab.std(ma.compressed())	#disregard the NAs
			list_of_stds.append(std)
		
		top_number = int(len(list_of_stds)*top_percentage)	#how many we want
		arg_list  = argsort(list_of_stds)	#sort it, ascending
		arg_list = arg_list.tolist()	#convert from array to list
		arg_list.reverse()	#reverse, descending order
		top_arg_list = arg_list[:top_number]	#get the top_number of arg_list
		if self.debug:
			print "list_of_stds is %s"%repr(list_of_stds)
			print "top_number is %s"%top_number
			print "arg_list is %s"%repr(arg_list)
			print "top_arg_list is %s"%repr(top_arg_list)
		list_of_top_mas = []
		for index in top_arg_list:
			list_of_top_mas.append(list_of_mas[index])
		sys.stderr.write("Done.\n")
		return list_of_top_mas
	
	def transpose_and_output(self, outfname, list_of_top_mas):
		"""
		05-09-05
			--ls_NA_fillin()
		"""
		sys.stderr.write("Outputing the data...")
		ls_2d = []
		for ma in list_of_top_mas:
			ls_2d.append(ma.raw_data())
		matrix = array(ls_2d)
		matrix = transpose(matrix)
		writer = csv.writer(open(outfname, 'w'), delimiter='\t')
		for i in range(matrix.shape[0]):
			ls_with_NA_filled = self.ls_NA_fillin(matrix[i])
			writer.writerow([i+1]+ls_with_NA_filled)
		
		sys.stderr.write("Done.\n")
	
	def ls_NA_fillin(self, ls):
		"""
		05-09-05
			1.1 is NA
		"""
		ls_to_return = []
		for item in ls:
			if item == 1.1:
				ls_to_return.append("NA")
			else:
				ls_to_return.append(item)
		return ls_to_return
	
	def run(self):
		self.list_of_mas = self.data_read_in(self.infname, self.no_of_nas)
		self.list_of_top_mas = self.get_top_mas(self.list_of_mas, self.top_percentage)
		self.transpose_and_output(self.outfname, self.list_of_top_mas)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "input_file=", "output_file=", "no_of_nas=", "top_percentage=", "debug"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "i:o:n:t:b", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	input_file = None
	output_file = None
	no_of_nas = None
	top_percentage = 0.1
	debug = 0

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-i", "--input_file"):
			input_file = arg
		elif opt in ("-o", "--output_file"):
			output_file = arg
		elif opt in ("-n", "--no_of_nas"):
			no_of_nas = int(arg)
		elif opt in ("-t", "--top_percentage"):
			top_percentage = float(arg)
		elif opt in ("-b", "--debug"):
			debug = 1

	if input_file and output_file:
		instance = PreprocessEdgeData(input_file, output_file, no_of_nas, top_percentage, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
