#!/usr/bin/env python
"""
Usage: LogAffyDatasets.py [OPTION] -o OUTPUTDIR FILES

Option:
	FILES are a list of files whose columns are going to be counted.
	-o ... --outputdir=...,	directory containing output files
	-d ... --delimiter=...,	delimiter character used to seperate columns, \t(default)
	-h, --help              show this help
	
Examples:
	LogAffyDatasets.py -o datasets/sc_log/ datasets/sc/*

Description:
	Program to log transform affymetrix datasets.
"""

import sys, os, re, getopt, csv, math

class LogAffyDatasets:
	def __init__(self, file_list, outputdir, delimiter):
		""" 
		07-31-05
		"""
		self.files = file_list
		self.files.sort()
		self.outputdir = outputdir
		self.delimiter = delimiter

	def transform_one_file(self, src_pathname, delimiter, outputdir):
		reader = csv.reader(file(src_pathname), delimiter=delimiter)
		filename = os.path.basename(src_pathname)
		output_filename = os.path.join(outputdir, filename)
		writer = csv.writer(open(output_filename, 'w'), delimiter=delimiter)
		for row in reader:
			new_row = [row[0]]
			for i in range(1, len(row)):
				if row[i] == 'NA' or row[i] == '':
					new_row.append(row[i])
				else:
					value = float(row[i])
					if value<=10:
						value = 10
					new_row.append(math.log(value, 2))
			writer.writerow(new_row)
		del reader, writer
	
	def run(self):
		"""
		07-31-05
		"""
		if not os.path.isdir(self.outputdir):
			os.makedirs(self.outputdir)
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(self.files))
		for f in self.files:
			sys.stderr.write("%d/%d:\t%s"%(self.files.index(f)+1,len(self.files),f))
			self.transform_one_file(f, self.delimiter, self.outputdir)
			sys.stderr.write("\n")


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "ho:d:", ["help", "outputdir=", "delimiter="])
	except:
		print __doc__
		sys.exit(2)

	delimiter = '\t'
	outputdir = None
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-o", "--outputdir"):
			outputdir = arg
		elif opt in ("-d", "--delimiter"):
			delimiter = arg

	if len(args)>=1:
		instance = LogAffyDatasets(args, outputdir, delimiter)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
