#!/usr/bin/env python
"""
Usage: DatasetsIdExchange.py [OPTION] -o OUTPUTDIR -m MAPPING_FILE FILES

Option:
	FILES are a list of files whose columns are going to be counted.
	-o ... --outputdir=...,	directory containing output files
	-d ... --delimiter=...,	delimiter character used to seperate columns, \t(default)
	-g ...,	organism, ('at', default)
	-m ...,	mapping file,(1st column is the target, 2nd column is source)
	-y ...,	type, (1 default)
	-b,	enable debugging
	-h, --help              show this help
	
Examples:
	DatasetsIdExchange.py -o datasets/at_gene/ -m gene2unigene datasets/at_unigene/smd*

Description:
	Program to exchange gene ids in datasets.
	type: 1 means given a list of tg_gene_id_list, just take the sorted-smallest one
		2 means outputs them all(copy)
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path += [os.path.expanduser('~/lib64/python')]
else:   #32bit
	sys.path += [os.path.expanduser('~/lib/python')]
import re, getopt, csv
from codense.common import get_unigene2gene_list, org2tax_id, org_short2long

class DatasetsIdExchange:
	def __init__(self, file_list=[], outputdir=None, delimiter='\t', organism='at', mapping_file=None, type=1, debug=0):
		""" 
		09-30-05
		"""
		self.files = file_list
		self.files.sort()
		self.outputdir = outputdir
		self.delimiter = delimiter
		self.organism = org_short2long(organism)
		self.mapping_file = mapping_file
		self.type = int(type)
		self.debug = int(debug)

	def transform_one_file(self, src_pathname, delimiter, outputdir, mapping_dict, type=1):
		"""
		"""
		reader = csv.reader(file(src_pathname), delimiter=delimiter)
		filename = os.path.basename(src_pathname)
		output_filename = os.path.join(outputdir, filename)
		writer = csv.writer(open(output_filename, 'w'), delimiter=delimiter)
		for row in reader:
			src_gene_id = row[0]
			if src_gene_id in mapping_dict:
				tg_gene_id_list = mapping_dict[src_gene_id]
				if type==1:
					tg_gene_id_list.sort()
					row[0] = tg_gene_id_list[0]	#only get the first one
					writer.writerow(row)
				elif type==2:
					for tg_gene_id in tg_gene_id_list:
						row[0] = tg_gene_id
						writer.writerow(row)
		del reader, writer
	
	def run(self):
		"""
		07-31-05
		"""
		if not os.path.isdir(self.outputdir):
			os.makedirs(self.outputdir)
		tax_id = org2tax_id(self.organism)
		mapping_dict = get_unigene2gene_list(self.mapping_file, tax_id)
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(self.files))
		for f in self.files:
			sys.stderr.write("%d/%d:\t%s"%(self.files.index(f)+1,len(self.files),f))
			self.transform_one_file(f, self.delimiter, self.outputdir, mapping_dict, self.type)
			sys.stderr.write("\n")


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "ho:d:g:m:y:u", ["help", "outputdir=", \
			"delimiter="])
	except:
		print __doc__
		sys.exit(2)

	delimiter = '\t'
	outputdir = None
	organism = 'at'
	mapping_file = None
	type = 1
	debug = 0
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-o", "--outputdir"):
			outputdir = arg
		elif opt in ("-d", "--delimiter"):
			delimiter = arg
		elif opt in ("-g"):
			organism = arg
		elif opt in ("-y"):
			type = int(arg)
		elif opt in ("-m"):
			mapping_file = arg
		elif opt in ("-b"):
			debug = 1
			
	if len(args)>=1 and mapping_file and outputdir:
		instance = DatasetsIdExchange(args, outputdir, delimiter, organism, mapping_file, type, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
