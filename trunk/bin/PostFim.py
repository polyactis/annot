#!/usr/bin/env python
"""
Usage: PostFim.py -i FIM_OUTPUT_FILE -o POSTFIM_OUTPUTFILE -m MIN_LENGTH [OPTIONS]

Option:
	-i ...,	fim_closed or closet+ output file
	-m ...,	minimum length of a frequent itemset
	-o ...,	post fim outputfile
	-y ...,	type, 1(fim_closed), 2(closet+) (IGNORE, only for closet+)
	-b,	debug version
	-r,	enable report flag
	-h, --help              show this help
	
Examples:
	PostFim.py -i mm_fim_114m5x114_o -m 5 -o mm_fim_114m5x114_o.postfim

Description:
	Filter those itemsets that have less than min_length items.
	The input is closet+'s output. Output is in fim_closed format,
		easy for MpiFromDatasetSignatureToPattern.py.
"""

import sys, os, psycopg, getopt, csv
from codense.common import db_connect
sys.path += [os.path.expanduser('~/script/annot/bin')]
if sys.version_info[:2] < (2, 3):       #python2.2 or lower needs some extra
	from python2_3 import *

def PostFim(fim_output_fname, output_fname, min_length, type=1, debug=0, report=0):
	sys.stderr.write("PostFim Filtering %s ...\n"%fim_output_fname)
	reader = csv.reader(open(fim_output_fname, 'r'), delimiter=' ')
	writer = csv.writer(open(output_fname, 'w'), delimiter=' ')
	counter = 0
	for row in reader:
		sig_row = row[:-2]
		frequency = int(row[-1])
		if len(sig_row)>=min_length:
			new_row = sig_row + ['(%s)'%frequency]
			writer.writerow(new_row)
		counter +=1
		if report and counter%5000==0:
			sys.stderr.write('%s%s'%('\x08'*20, counter))
	del reader, writer
	sys.stderr.write("PostFim Done.\n")


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hi:o:m:y:br", ["help"])
	except:
		print __doc__
		sys.exit(2)
	
	fim_output_fname = None
	output_fname = None
	min_length = None
	type = 1
	debug = 0
	report = 0
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-i",):
			fim_output_fname = arg
		elif opt in ("-o",):
			output_fname = arg
		elif opt in ("-m",):
			min_length = int(arg)
		elif opt in ("-y",):
			type = int(arg)
		elif opt in ("-b",):
			debug = 1
		elif opt in ("-r",):
			report = 1
	if fim_output_fname and output_fname and min_length:
		PostFim(fim_output_fname, output_fname, min_length, type, \
			debug, report)
	else:
		print __doc__
		sys.exit(2)
