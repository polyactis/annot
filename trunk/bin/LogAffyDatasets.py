#!/usr/bin/env python
"""
Usage: LogAffyDatasets.py [OPTION] -o OUTPUTDIR FILES

Option:
	FILES are a list of files whose columns are going to be counted.
	-o ... --outputdir=...,	directory containing output files
	-d ... --delimiter=...,	delimiter character used to seperate columns, \t(default)
	-s ... --threshold=...,	std/mean threshold, 1.0 (default)
	-y ... --type=...,	geo(1,default), or smd(2)
	-l, --log	apply log transform (IGNORE)
	-u, --debug	enable debugging
	-h, --help              show this help
	
Examples:
	LogAffyDatasets.py -o datasets/sc_log/ datasets/sc/*

Description:
	Program to log transform affymetrix datasets.
	geo: log +std/mean, smd: std
"""

import sys, os, re, getopt, csv, math
import MLab
from MA import array
from Preprocess import PreprocessEdgeData

class LogAffyDatasets:
	def __init__(self, file_list, outputdir, delimiter, threshold, type=1, log=0, debug=0):
		""" 
		07-31-05
		08-09-05
			add type
		"""
		self.files = file_list
		self.files.sort()
		self.outputdir = outputdir
		self.delimiter = delimiter
		self.threshold = float(threshold)
		self.type = int(type)
		self.log = int(log)
		self.debug = int(debug)

	def transform_one_file(self, src_pathname, delimiter, outputdir, b_instance, threshold, type):
		"""
		08-09-05
			add type
		"""
		reader = csv.reader(file(src_pathname), delimiter=delimiter)
		filename = os.path.basename(src_pathname)
		output_filename = os.path.join(outputdir, filename)
		writer = csv.writer(open(output_filename, 'w'), delimiter=delimiter)
		for row in reader:
			gene_id = row[0]
			new_row = []
			mask_ls = []
			for i in range(1, len(row)):
				if row[i] == 'NA':
					new_row.append(1e20)
					mask_ls.append(1)
				elif row[i] == '':
					#ignore empty entry
					continue
				else:
					value = float(row[i])
					if type==1:
						if value<=10:
							value = 10
						value = math.log(value, 2)
					new_row.append(value)
					mask_ls.append(0)
			ma_array = array(new_row, mask=mask_ls)
			if self.debug:
				print "The data vector is ",ma_array
				print "Its mask is ", ma_array.mask()
			if len(ma_array.compressed())>1:	#at least two samples, otherwise, correlation can't be calculated
				std = MLab.std(ma_array.compressed())	#disregard the NAs
				if type==1:
					ratio = std/MLab.mean(ma_array.compressed())
				elif type==2:
					ratio = std
				else:
					sys.stderr.write("type %s not supported\n"%type)
					sys.exit(2)
				if self.debug:
					print "std is ",std
					print "ratio is ", ratio
					raw_input("Continue?(Y/n)")
				if ratio >= threshold:
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
			self.transform_one_file(f, self.delimiter, self.outputdir, b_instance, self.threshold, self.type)
			sys.stderr.write("\n")


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "ho:d:s:y:lu", ["help", "outputdir=", "delimiter=", "threshold=", "type=", "log", "debug"])
	except:
		print __doc__
		sys.exit(2)

	delimiter = '\t'
	outputdir = None
	threshold = 1.0
	type = 1
	log = 0
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
		elif opt in ("-y", "--type"):
			type = int(arg)
		elif opt in ("-l", "--log"):
			log = 1
		elif opt in ("-u", "--debug"):
			debug = 1
			
	if len(args)>=1:
		instance = LogAffyDatasets(args, outputdir, delimiter, threshold, type, log, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
