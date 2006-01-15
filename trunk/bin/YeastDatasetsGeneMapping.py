#!/usr/bin/env python
"""
Usage: YeastDatasetsGeneMapping.py [OPTION] -o OUTPUTDIR -m MAP_FILE INPUTFILES

Option:
	FILES are a list of files whose columns are going to be counted.
	-o ... --outputdir=...,	directory containing output files
	-m ..., --map_file=...	the file containing orfname to gene_id mapping
	-b, --debug	enable debugging
	-h, --help              show this help
	
Examples:
	YeastDatasetsGeneMapping.py -o ~/datasets/sc_gene_all -m /tmp/yeast_orf2gene ~/datasets/sc/*

Description:
	Map the orf in yeast datasets into gene_id.
	Source of the map_file: 1. get all orfnames from the datasets
	2. link against the gene.gene_symbol2id table to find the gene id
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import sys, os, re, getopt, csv, math

class YeastDatasetsGeneMapping:
	def __init__(self, file_list, outputdir, map_file=None, debug=0):
		""" 
		08-17-05
		"""
		self.files = file_list
		self.files.sort()
		self.outputdir = outputdir
		self.map_file = map_file
		self.debug = int(debug)
	
	def load_orf2gene(self, map_file):
		"""
		08-17-05
		"""
		sys.stderr.write("Loading orf2gene...")
		orf2gene = {}
		reader = csv.reader(open(map_file, 'r'), delimiter='\t')
		for row in reader:
			orf = row[0]
			gene_id = int(row[1])
			if orf in orf2gene:
				sys.stderr.write("%s appears more than once in the map_file.\n"%orf)
				sys.exit(2)
			orf2gene[orf] = gene_id
		del reader
		sys.stderr.write("Done.\n")
		return orf2gene
	
	def gene_map_one_file(self, orf2gene, src_pathname, outputdir):
		"""
		08-17-05
		"""
		reader = csv.reader(file(src_pathname), delimiter='\t')
		filename = os.path.basename(src_pathname)
		output_filename = os.path.join(outputdir, filename)
		writer = csv.writer(open(output_filename, 'w'), delimiter='\t')
		for row in reader:
			orf = row[0]
			if orf in orf2gene:
				writer.writerow([orf2gene[orf]] + row[1:])
		del reader, writer
	
	def run(self):
		"""
		08-17-05
		"""
		if not os.path.isdir(self.outputdir):
			os.makedirs(self.outputdir)
		orf2gene = self.load_orf2gene(self.map_file)
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(self.files))
		for f in self.files:
			sys.stderr.write("%d/%d:\t%s"%(self.files.index(f)+1,len(self.files),f))
			self.gene_map_one_file(orf2gene, f, self.outputdir)
			sys.stderr.write("\n")


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "ho:m:b", ["help", "outputdir=", "map_file=", "debug"])
	except:
		print __doc__
		sys.exit(2)

	delimiter = '\t'
	outputdir = None
	map_file = None
	debug = 0
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-o", "--outputdir"):
			outputdir = arg
		elif opt in ("-m", "--map_file"):
			map_file = arg
		elif opt in ("-b", "--debug"):
			debug = 1
			
	if outputdir and map_file and len(args)>=1:
		instance = YeastDatasetsGeneMapping(args, outputdir, map_file, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
