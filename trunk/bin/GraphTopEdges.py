#!/usr/bin/env python
"""
Usage: GraphTopEdges.py -k SCHEMA -d DATASET_FILE -i INPUT_FILE [OPTIONS]

Option:
	-d ..., --dataset_file=...	DATASET_FILE, is the file containing the expression values
	-i ..., --input_file=...	graph_modeling.cc output file
	-o ..., --output_file=...	haiyan's matrix format
	-t ..., --top_percentage=...	0.01(default)
	-b, --debug	just fetch 15000 edges and do a demo
	-h, --help		show this help

Examples:
	GraphTopEdges.py -d mm_GDS246 -i gph_mm_GDS246 -o gph_mm_GDS246.2
	
Description:
	Get the top percentage edges.
	
"""

import sys, os, getopt, csv, re, random
import MLab
from Numeric import argsort, take, transpose
from MA import array

class GraphTopEdges:
	"""
	06-21-05
		a module to filter the edges of a graph,
		1. only take the top percentage edges
	"""
	def __init__(self, dataset_file=None, infname=None, outfname=None, top_percentage=0.01, debug=0):
		self.dataset_file = dataset_file
		self.infname = infname
		self.outfname = outfname
		self.top_percentage = top_percentage
		self.debug = int(debug)
		
	def data_read_in(self, infname):
		"""
		06-21-05
		"""
		sys.stderr.write("Reading data...")
		reader = csv.reader(open(infname, 'r'),delimiter='\t')
		header_row = reader.next()	#the first line is 't	#	gphName'
		graph_dict = {}
		for row in reader:
			edge_tuple = tuple(row[1:3])	#Note, the 2nd and 3rd are denoted by 1:3
			graph_dict[edge_tuple] = float(row[3])

		del reader
		sys.stderr.write("Done.\n")
		return (header_row, graph_dict)
	
	def get_top_edges_and_output(self, graph_dict, top_number, outfname, header_row):
		"""
		06-21-05
			If the number of edges is less than the top_number, take it directly.
		"""
		sys.stderr.write("Getting the top %s edges..."%top_percentage)
		edge_tuple_list = graph_dict.keys()
		cor_list = graph_dict.values()
		arg_cor_list  = argsort(cor_list)	#sort it, ascending
		arg_cor_list = arg_cor_list.tolist()	#convert from array to list
		arg_cor_list.reverse()	#reverse, descending order
		top_arg_list = arg_cor_list[:top_number]	#get the top_number of arg_list	#06-07-05 if top_number>len(arg_list), it's ok.
		if self.debug:
			print "cor_list is %s"%repr(cor_list)
			print "top_number is %s"%top_number
			print "arg_cor_list is %s"%repr(arg_cor_list)
			print "top_arg_list is %s"%repr(top_arg_list)
		writer = csv.writer(open(outfname, 'w'), delimiter='\t')
		writer.writerow(header_row)
		for index in top_arg_list:
			writer.writerow(['e']+ list(edge_tuple_list[index])+ [cor_list[index]])
		del writer
		sys.stderr.write("Done.\n")
	
	def get_top_number(self, dataset_file, top_percentage):
		"""
		06-21-05
			count the number of lines in dataset_file
		"""
		sys.stderr.write("Getting the top number...")
		inf = open(dataset_file, 'r')
		i = 0
		for line in inf:
			i += 1
		sys.stderr.write("Done.\n")
		return int(i*i*top_percentage)
	
	def run(self):
		"""
		06-07-05
			--data_read_in()
			--get_top_number()
			--get_top_edges_and_output()
		"""
		header_row, graph_dict = self.data_read_in(self.infname)
		top_number = self.get_top_number(self.dataset_file, self.top_percentage)
		self.get_top_edges_and_output(graph_dict, top_number, self.outfname, header_row)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "input_file=", "output_file=", "dataset_file=", "top_percentage=", "debug"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "i:o:d:t:b", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	input_file = None
	output_file = None
	dataset_file = None
	top_percentage = 0.01
	debug = 0

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-i", "--input_file"):
			input_file = arg
		elif opt in ("-o", "--output_file"):
			output_file = arg
		elif opt in ("-d", "--dataset_file"):
			dataset_file = arg
		elif opt in ("-t", "--top_percentage"):
			top_percentage = float(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
	if output_file==None and input_file:
		output_file = input_file+'.2'
	if input_file and output_file and dataset_file:
		instance = GraphTopEdges(dataset_file, input_file, output_file, top_percentage, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
