#!/usr/bin/env python
"""
Usage: PreprocessEdgeData.py -k SCHEMA [OPTION] R_FILE

Option:
	R_FILE is the file to store the R code.
	-i ..., --input_file=...	gspan format
	-o ..., --output_file=...	haiyan's matrix format
	-n ..., --no_of_nas=...	maximum number of NAs of one edge, None(default)
	-t ..., --top_percentage=...	0.1(default)
	-p, --plain_na	use plain 'NA' instead of random number
	-b, --debug	just fetch 15000 edges and do a demo
	-h, --help		show this help

Examples:
	PreprocessEdgeData.py -i edge_data_go_282 -o edge_data_go_282.1
	
Description:

	
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import sys, os, getopt, csv, re, random
import MLab
from Numeric import argsort, take, transpose
from MA import array, filled

class PreprocessEdgeData:
	"""
	05-09-05
		a module to preprocess the edge data of a function category, which
		is the correlation vector of an edge across all datasets.
		1. only take the 10% edges with highest std
		2. NA <= 7
		3. transpose the matrix
	
	06-15-05
		modify it to verify the filtering data step in iArrayAnalyzer.
	
	"""
	def __init__(self, infname=None, outfname=None, no_of_nas=None, \
		top_percentage=0.1, plain_na=0, debug=0):
		self.infname = infname
		self.outfname = outfname
		self.no_of_nas = no_of_nas
		self.top_percentage = top_percentage
		self.plain_na = int(plain_na)
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
		list_of_gene_ids = []
		for row in reader:
			data_ls = []
			mask_ls = []
			for item in row[1:]:	#ignore the first edge id
				if item=='NA':
					data_ls.append(1e20)
					mask_ls.append(1)
				else:
					data_ls.append(float(item))
					mask_ls.append(0)
			if no_of_nas:
				if sum(mask_ls)>no_of_nas:	#too many NAs
					continue
			
			list_of_gene_ids.append(row[0])
			list_of_mas.append(array(data_ls, mask=mask_ls))
		"""
		#the rest NA replaced with mean
		for i in range(len(list_of_mas)):
			ma = list_of_mas[i]
			list_of_mas[i]  = filled(ma, MLab.mean(ma.compressed()))
		"""
		del reader
		sys.stderr.write("Done.\n")
		return list_of_mas, list_of_gene_ids
	
	def get_top_mas(self, list_of_mas, list_of_gene_ids, threshold, outfname):
		"""
		05-09-05
			start
		06-07-05
			if top_percentage is less than 200, use 200.
		"""
		sys.stderr.write("Getting the std/mean >=%s edges..."%threshold)
		writer = csv.writer(open(outfname, 'w'), delimiter='\t')
		print len(list_of_mas), len(list_of_gene_ids)
		for i in range(len(list_of_mas)):
			ma = list_of_mas[i]
			std = MLab.std(ma.compressed())	#disregard the NAs
			value = std/MLab.mean(ma.compressed())
			if self.debug:
				print list_of_gene_ids[i]
				print ma
				print std
				print value
				raw_input("y/n")
			if value>=threshold:
				ls_with_NA_filled = self.ls_NA_fillin(ma)
				writer.writerow([list_of_gene_ids[i]]+ls_with_NA_filled)
		del writer
		sys.stderr.write("Done.\n")
	

	def ls_NA_fillin(self, ls):
		"""
		06-17-05
		"""
		ls_to_return = []
		mask_ls = ls.mask()
		for i in range(len(ls)):
			if mask_ls[i] == 1:
				ls_to_return.append('NA')
			else:
				ls_to_return.append(ls[i])
		return ls_to_return
	
	def run(self):
		"""
		06-07-05
			--data_read_in()
			--get_top_mas()
				--ls_NA_fillin()
		"""
		list_of_mas, list_of_gene_ids = self.data_read_in(self.infname, self.no_of_nas)
		print len(list_of_mas)
		self.get_top_mas(list_of_mas, list_of_gene_ids, self.top_percentage, self.outfname)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "input_file=", "output_file=", "no_of_nas=", "top_percentage=", "plain_na", "debug"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "i:o:n:t:pb", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	input_file = None
	output_file = None
	no_of_nas = None
	top_percentage = 0.1
	plain_na = 0
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
		elif opt in ("-p", "--plain_na"):
			plain_na = 1
		elif opt in ("-b", "--debug"):
			debug = 1

	if input_file and output_file:
		instance = PreprocessEdgeData(input_file, output_file, no_of_nas, top_percentage, plain_na, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
