#!/usr/bin/env python
"""
Usage: DrawHistStd.py [OPTION] -o OUTPUTDIR FILES

Option:
	FILES are a list of files whose columns are going to be counted.
	-o ... --outputdir=...,	directory containing output files
	-d ... --delimiter=...,	delimiter character used to seperate columns, \t(default)
	-y ... --type=...,	geo(1,default, log before std), or smd(2)
	-n ... --no_of_valids=...	min number of non-NA, 8(default)
	-u, --debug	enable debugging
	-h, --help              show this help
	
Examples:
	DrawHistStd.py -o datasets/sc_log/ datasets/sc/*

Description:
	DrawHistStd.py is usd to draw histogram of standard deviations for each microarray dataset.
"""

import sys, os, re, getopt, csv, math
import MLab
from MA import array
from Preprocess import PreprocessEdgeData
from rpy import r

class DrawHistStd:
	def __init__(self, file_list, outputdir, delimiter, threshold, type=1, no_of_valids=8, debug=0):
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
		self.no_of_valids = int(no_of_valids)
		self.debug = int(debug)

	def transform_one_file(self, src_pathname, delimiter, outputdir, b_instance, threshold, type, no_of_valids):
		"""
		08-09-05
			add type
		08-29-05
			add no_of_valids to cut genes with too few valid values
		"""
		reader = csv.reader(file(src_pathname), delimiter=delimiter)
		filename = os.path.basename(src_pathname)
		output_filename = os.path.join(outputdir, filename)
		std_list = []
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
						value = math.log(value)
					new_row.append(value)
					mask_ls.append(0)
			ma_array = array(new_row, mask=mask_ls)
			if self.debug:
				print "The data vector is ",ma_array
				print "Its mask is ", ma_array.mask()
			if len(ma_array.compressed())>=no_of_valids:	#at least two samples, otherwise, correlation can't be calculated
				#08-29-05	no_of_valids controls not too many NA's, which is for graph_modeling
				std = MLab.std(ma_array.compressed())	#disregard the NAs
				if self.debug:
					print "std is ",std
					raw_input("Continue?(Y/n)")
				std_list.append(std)
		del reader
		if len(std_list)>100:
			r.png('%s.png'%output_filename)
			r.hist(std_list, main='histogram',xlab='std',ylab='freq')
			r.dev_off()
	
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
			self.transform_one_file(f, self.delimiter, self.outputdir, b_instance, self.threshold, self.type, self.no_of_valids)
			sys.stderr.write("\n")


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "ho:d:s:y:n:u", ["help", "outputdir=", \
			"delimiter=", "threshold=", "type=", "no_of_valids=", "debug"])
	except:
		print __doc__
		sys.exit(2)

	delimiter = '\t'
	outputdir = None
	threshold = 1.0
	type = 1
	no_of_valids = 8
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
		elif opt in ("-n", "--no_of_valids"):
			no_of_valids = int(arg)
		elif opt in ("-u", "--debug"):
			debug = 1
			
	if len(args)>=1:
		instance = DrawHistStd(args, outputdir, delimiter, threshold, type, no_of_valids, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
