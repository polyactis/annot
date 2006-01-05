#!/usr/bin/env python
"""
Usage: DrawHistCor.py [OPTION] -o OUTPUTDIR FILES

Option:
	FILES are a list of files whose columns are going to be counted.
	-o ... --outputdir=...,	directory containing output files
	-d ... --delimiter=...,	delimiter character used to seperate columns, \t(default)
	-u, --debug	enable debugging
	-h, --help              show this help
	
Examples:
	DrawHistCor.py -o gph_result/sc_log/ gph_result/sc/*

Description:
	draw histogram of pairwise correlations generated from a microarray dataset.
	(output of MpiGrahModeling.py)
"""

import sys, os, re, getopt, csv, math
import MLab
from MA import array
from Preprocess import PreprocessEdgeData
from rpy import r

class DrawHistCor:
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
		01-05-06
			deal with blank files
		"""
		reader = csv.reader(file(src_pathname), delimiter=delimiter)
		filename = os.path.basename(src_pathname)
		output_filename = os.path.join(outputdir, filename)
		cor_list = []
		counter=0	#01-05-06
		for row in reader:
			if counter>0:
				cor = float(row[3])
				cor_list.append(cor)
			counter += 1
		del reader
		if len(cor_list)>100:
			r.png('%s.png'%output_filename)
			r.hist(cor_list, main='histogram',xlab='cor',ylab='freq')
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
		instance = DrawHistCor(args, outputdir, delimiter, threshold, type, no_of_valids, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
