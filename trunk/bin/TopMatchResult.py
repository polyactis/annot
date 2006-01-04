#!/usr/bin/env python
"""
Usage: TopMatchResult.py -i -o [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database(transfac)
	-i ...,	inputdir, TRANSFAC match output files
	-o ...,	output file
	-t ...,	top number, 2000(default)
	-b,	debug version.
	-r,	enable report flag
	-h, --help	Display the usage infomation.
	
Examples:
	~/script/annot/bin/TopMatchResult.py -i -o -t 
	
Description:
	The program reads the MATCH output files and only retain the top number of binding sites.
	
"""

import sys, os, getopt, csv, cPickle, fileinput
sys.path += [os.path.expanduser('~/script/annot/bin')]
sys.path += [os.path.expanduser('~/script/transfac/src')]
from codense.common import db_connect
from binding_site2gene_id2mt_no import attr_of_mt_no	#01-03-06	for the top number

class TopMatchResult:
	def __init__(self,hostname='zhoudb', dbname='graphdb', schema=None, inputdir=None,\
		output_file=None, top_number=2000, debug=0, report=0):
		"""
		11-15-05
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.inputdir = inputdir
		self.output_file = output_file
		self.top_number = int(top_number)
		self.debug = int(debug)
		self.report = int(report)
	
	def get_mt_id2binding_sites(self, inputdir, top_number):
		input_files = os.listdir(inputdir)
		no_of_input_files = len(input_files)
		sys.stderr.write("Totally, %s files to handle.\n"%no_of_input_files)
		mt_id2binding_sites = {}
		hit_id = 0	#01-03-06 an id used to track all binding sites
		for i in range(no_of_input_files):
			fname = input_files[i]
			sys.stderr.write("\t %s/%s: %s"%(i+1, no_of_input_files, fname))
			inf = open(os.path.join(inputdir, fname), 'r')
			for line in inf:
				if line[:10]=='Inspecting':
					prom_id = int(line.split()[-1])	#\n is discarded automatically by split()
				elif line[0] == ' ' and line[:6]!=' Total' and line[:6]!=' Frequ':	#the first character is blank, but exclude the end statistic part	
					ls = line[:-1].split('|')
					mt_id = ls[0].strip()	#remove spaces
					bs_disp_start_strand = ls[1].strip()
					bs_disp_start = int(bs_disp_start_strand[:-3])
					strand = bs_disp_start_strand[-2]
					core_similarity_score = float(ls[2])
					matrix_similarity_score = float(ls[3])
					sequence = ls[4].strip()
					#01-03-06
					bs_disp_end = bs_disp_start + len(sequence) - 1
				
					if mt_id not in mt_id2binding_sites:
						mt_id2binding_sites[mt_id] = attr_of_mt_no(top_number, self.debug)
					hit_tuple = (matrix_similarity_score, prom_id, strand, bs_disp_start, bs_disp_end, core_similarity_score, sequence)
					mt_id2binding_sites[mt_id].consume_new_row(hit_tuple, prom_id, hit_id)	#01-03-06 hit_tuple replaces matrix_similarity_score
					hit_id +=1
			del inf
			sys.stderr.write(".\n")
		return mt_id2binding_sites
	
	def output_mt_id2binding_sites(self, output_file, mt_id2binding_sites):
		"""
		01-03-06
		"""
		sys.stderr.write("Outputting mt_id2binding_sites")
		writer = csv.writer(open(output_file, 'w'), delimiter='\t')
		perfect_matrices = []
		for mt_id, mt_id_attr in mt_id2binding_sites.iteritems():
			lowest_score = mt_id_attr.hq_list[0][0][0]
			if lowest_score>=1.0:
				perfect_matrices.append(mt_id)
				continue
			for hit_tuple, prom_id, hit_id in mt_id_attr.hq_list:
				matrix_similarity_score, prom_id, strand, bs_disp_start, bs_disp_end, core_similarity_score, sequence = hit_tuple
				row = [mt_id, prom_id, strand, bs_disp_start, bs_disp_end, core_similarity_score, matrix_similarity_score, sequence]
				writer.writerow(row)
		del writer
		sys.stderr.write(".\n")
		print "perfect_matrices:", perfect_matrices
	
	def run(self):
		"""
		01-03-06
		"""
		mt_id2binding_sites = self.get_mt_id2binding_sites(self.inputdir, self.top_number)
		self.output_mt_id2binding_sites(self.output_file, mt_id2binding_sites)
		
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:o:t:br", ["help", \
			"hostname=", "dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = 'transfac'
	inputdir = None
	output_file = None
	top_number = 2000
	debug = 0
	report = 0
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-z", "--hostname"):
			hostname = arg
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-k", "--schema"):
			schema = arg
		elif opt in ("-i",):
			inputdir = arg
		elif opt in ("-o",):
			output_file = arg
		elif opt in ("-t",):
			top_number = int(arg)
		elif opt in ("-b",):
			debug = 1
		elif opt in ("-r",):
			report = 1
	if inputdir and output_file:
		instance = TopMatchResult(hostname, dbname, schema, inputdir, output_file, top_number, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
