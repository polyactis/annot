#!/usr/bin/env python
"""
Usage: LogAffyDatasets.py [OPTION] -o OUTPUTDIR FILES

Option:
	FILES are a list of files whose columns are going to be counted.
	-o ... --outputdir=...,	directory containing output files
	-d ... --delimiter=...,	delimiter character used to seperate columns, \t(default)
	-s ... --threshold=...,	std/mean threshold, 1.0 (default)
	-t ...,	top percentage	0.90(default)
	-n ... --no_of_valids=...	min number of non-NA, 8(default)
	-l,	take log of the value
	-m,	divide the std with mean
	-u, --debug	enable debugging
	-h, --help              show this help
	
Examples:
	LogAffyDatasets.py -o datasets/sc_log/ datasets/sc/*

Description:
	Program to filter datasets based on std or std/mean.

"""

import sys, os, re, getopt, csv, math
import MLab
from MA import array
from Preprocess import PreprocessEdgeData
from sets import Set

class LogAffyDatasets:
	def __init__(self, file_list, outputdir, delimiter, threshold, top_percentage=0.90, no_of_valids=8, take_log=0, divide_mean=0, debug=0):
		""" 
		07-31-05
		08-09-05
			add type
		12-22-05
			add top_percentage
			add take_log, divide_mean
		"""
		self.files = file_list
		self.files.sort()
		self.outputdir = outputdir
		self.delimiter = delimiter
		self.threshold = float(threshold)
		self.top_percentage = float(top_percentage)
		self.no_of_valids = int(no_of_valids)
		self.take_log = int(take_log)
		self.divide_mean = int(divide_mean)
		self.debug = int(debug)
	
	def get_qualified_counter_set(self, std_counter_ls, top_percentage):
		"""
		12-22-05
			take the top_percentage counter's
		"""
		qualified_counter_set = Set()
		std_counter_ls.sort()
		std_counter_ls.reverse()
		for i in range(int(len(std_counter_ls)*top_percentage)):
			if self.debug:
				print std_counter_ls[i][0], std_counter_ls[i][1]
				#raw_input("Continue?(Y/n)")
			qualified_counter_set.add(std_counter_ls[i][1])
		return qualified_counter_set
	
	def get_ma_array_out_of_list(self, expr_list, take_log):
		"""
		12-22-05
		"""
		new_row = []
		mask_ls = []
		for i in range(len(expr_list)):
			if expr_list[i] == 'NA':
				new_row.append(1e20)
				mask_ls.append(1)
			elif expr_list[i] == '':
				#ignore empty entry
				continue
			else:
				value = float(expr_list[i])
				if take_log:	#12-22-05
					if value<=10:
						value = 10
					value = math.log(value)	#12-22-05
				new_row.append(value)
				mask_ls.append(0)
		ma_array = array(new_row, mask=mask_ls)
		return ma_array
	
	def transform_one_file(self, src_pathname, delimiter, outputdir, b_instance, threshold, no_of_valids, top_percentage, \
		take_log, divide_mean):
		"""
		08-09-05
			add type
		08-29-05
			add no_of_valids to cut genes with too few valid values
		12-22-05
			add top_percentage
			change log(x,2) to log(x) (natural number is base)
		"""
		reader = csv.reader(file(src_pathname), delimiter=delimiter)
		#1st round to read
		counter = 0
		std_counter_ls = []
		for row in reader:
			counter += 1
			gene_id = row[0]
			ma_array = self.get_ma_array_out_of_list(row[1:], take_log)
			"""
			if self.debug:
				print "The data vector is ",ma_array
				print "Its mask is ", ma_array.mask()
			"""
			if len(ma_array.compressed())>=no_of_valids:	#at least two samples, otherwise, correlation can't be calculated
				#08-29-05	no_of_valids controls not too many NA's, which is for graph_modeling
				std = MLab.std(ma_array.compressed())	#disregard the NAs
				if divide_mean:
					ratio = std/MLab.mean(ma_array.compressed())
				else:
					ratio = std
				"""
				if self.debug:
					print "std is ",std
					print "ratio is ", ratio
					raw_input("Continue?(Y/n)")
				"""
				std_counter_ls.append([ratio, counter])
		del reader
		
		qualified_counter_set = self.get_qualified_counter_set(std_counter_ls, top_percentage)
		
		#2nd round to read,  and write out
		reader = csv.reader(file(src_pathname), delimiter=delimiter)
		filename = os.path.basename(src_pathname)
		output_filename = os.path.join(outputdir, filename)
		writer = csv.writer(open(output_filename, 'w'), delimiter=delimiter)
		counter = 0
		for row in reader:
			counter += 1
			if counter in qualified_counter_set:
				gene_id = row[0]
				ma_array = self.get_ma_array_out_of_list(row[1:], take_log)
				writer.writerow([gene_id] + b_instance.ls_NA_fillin(ma_array))
		del reader, writer
	
	def run(self):
		"""
		07-31-05
		"""
		b_instance = PreprocessEdgeData()
		if not os.path.isdir(self.outputdir):
			os.makedirs(self.outputdir)
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(self.files))
		for f in self.files:
			sys.stderr.write("%d/%d:\t%s"%(self.files.index(f)+1,len(self.files),f))
			self.transform_one_file(f, self.delimiter, self.outputdir, b_instance, self.threshold, self.no_of_valids, self.top_percentage, \
				self.take_log, self.divide_mean)
			sys.stderr.write("\n")


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "ho:d:s:t:n:lmu", ["help", "outputdir=", \
			"delimiter=", "threshold=", "type=", "no_of_valids=", "debug"])
	except:
		print __doc__
		sys.exit(2)

	delimiter = '\t'
	outputdir = None
	threshold = 1.0
	top_percentage = 0.90
	no_of_valids = 8
	take_log = 0
	divide_mean = 0
	debug = 0
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-o", "--outputdir"):
			outputdir = arg
		elif opt in ("-d", "--delimiter"):
			delimiter = arg
		elif opt in ("-s", "--threshold"):
			threshold = float(arg)
		elif opt in ("-t"):
			top_percentage = float(arg)
		elif opt in ("-n", "--no_of_valids"):
			no_of_valids = int(arg)
		elif opt in ("-l"):
			take_log = 1
		elif opt in ("-m"):
			divide_mean = 1
		elif opt in ("-u", "--debug"):
			debug = 1
			
	if len(args)>=1:
		instance = LogAffyDatasets(args, outputdir, delimiter, threshold, top_percentage, no_of_valids, take_log, divide_mean, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
