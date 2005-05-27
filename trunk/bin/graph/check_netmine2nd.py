#!/usr/bin/env python
"""
Usage: check_netmine2nd.py  -i INPUT_FILE -v COR_VECTOR_FILE -o OUTPUTFILE [OPTION]

Option:
	-i ..., --input_file=...	gspan format(by graph_merge.py)
	-w ..., --which_line=...	which line of input_file, 0(default)
	-v ..., --cor_vector_file=...	the edge correlation vector file
	-s ..., --sig_vector_file=...	the edge significance flag vector file.
	-o ..., --output_file=...	the output_file
	-b ..., --label=...	use gene 1(index, default) or 2(no) or 3(id) to label
	-p ..., --p_value_cut_off=...	the p_value_cut_off for an edge to
		be significant, 0.01(default)
	-c ..., --cor_cut_off=...	the cor_cut_off for an edge to 
		be significant, 0.6(default)
		priority gives p_value_cut_off, cor_cut_off is used
		only when p_value_cut_off=0
	-r, --report	report the progress(a number)
	-u, --debug	enable debugging output.
	-h, --help              show this help
	
Examples:
	check_netmine2nd.py -i mm79_5default_from_edge_r4 -o /tmp/gph.2nd -p 0.001 -v mm_79_5.cor_vector

Description:
	(05-26-05)This program is used to check whether the 2nd-order graph constructed by netmine2nd is correct or not.

"""

import sys, os, getopt, csv, re, math
import graph_modeling
from sets import Set

class check_netmine2nd:
	"""
	05-26-05
	"""
	def __init__(self, input_file=None, which_line=0, cor_vector_file=None, sig_vector_file=None,\
		output_file=None, label=None, p_value_cut_off=0.01, cor_cut_off=0.6, report=0, debug=0):
		self.input_file = input_file
		self.which_line = int(which_line)
		self.cor_vector_file = cor_vector_file
		self.sig_vector_file = sig_vector_file
		self.output_file = output_file
		self.label = label
		self.p_value_cut_off = float(p_value_cut_off)
		self.cor_cut_off = float(cor_cut_off)
		self.report = int(report)
		self.debug = int(debug)
	
		self.min_valid_column = 6
		
	def vertex_set_readin(self, input_file, which_line):
		"""
		05-26-05
			read vertices from a specified line and return as a set
		"""
		sys.stderr.write("Getting the vertex_set...")
		reader = csv.reader(open(input_file, 'r'), delimiter='\t')
		i = 0
		for row in reader:
			if i==which_line:
				row = map(int, row)
				vertex_set = Set(row[3:])	#1st is number, 2nd is no of vertices, 3rd is no of edges.
				break
			i+=1
		del reader
		sys.stderr.write("Done.\n")
		if self.debug:
			sys.stderr.write("The vertex_set is %s\n"%(repr(vertex_set)))
		return vertex_set
	
	def cor_vector_readin(self, cor_vector_file, vertex_set):
		"""
		05-26-05
			scan through the cor_vector_file and select the edges whose both vertices
			appear in the vertex_set
		"""
		sys.stderr.write("Getting cor_vector ...")
		reader = csv.reader(open(cor_vector_file, 'r'), delimiter='\t')
		edge_tuple2cor_vector = {}
		for row in reader:
			row = map(int, row)
			if row[0] in vertex_set and row[1] in vertex_set:
				if row[0]<row[1]:
					edge_tuple = (row[0], row[1])
				else:
					edge_tuple = (row[1], row[0])
				if edge_tuple not in edge_tuple2cor_vector:
					edge_tuple2cor_vector[edge_tuple] = []
				for cor in row[2:]:
					if cor==1100:
						edge_tuple2cor_vector[edge_tuple].append(float(100000000))
						#100000000 is NA,
					else:
						edge_tuple2cor_vector[edge_tuple].append(cor/1000.0)
						#the cor_vector_file stores the first three digit after the dot.
						#dividing it by 1000 gives the float correlation
				if self.debug:
					sys.stderr.write("%s cor vector: %s\n"%(repr(edge_tuple), repr(edge_tuple2cor_vector[edge_tuple])))
		sys.stderr.write("Done.\n")
		del reader
		return edge_tuple2cor_vector
	
	def construct_2nd_graph(self, output_file, edge_tuple2cor_vector, p_value_cut_off, min_valid_column, cor_cut_off_vector):
		"""
		05-26-05
			call graph_modeling to do the real calculation, not leave one out.
		"""
		sys.stderr.write("Constructing 2nd-order graph ...")
		writer = csv.writer(open(output_file,'w'), delimiter='\t')
		edge_tuple_list = edge_tuple2cor_vector.keys()
		edge_tuple_list.sort()
		for i in range(len(edge_tuple_list)):
			for j in range(i+1, len(edge_tuple_list)):
				edge_tuple1 = edge_tuple_list[i]
				edge_tuple2 = edge_tuple_list[j]
				edge_data = graph_modeling.ind_cor(edge_tuple2cor_vector[edge_tuple1],\
					edge_tuple2cor_vector[edge_tuple2], -1)	#-1 means no leave one out.
				if self.debug:
					sys.stderr.write("edge: %s-%s; degree: %s; value: %s.\n"%(repr(edge_tuple1),\
						repr(edge_tuple2), edge_data.degree, edge_data.value))
					sys.stderr.write("The cor cur off is %s\n"%(cor_cut_off_vector[edge_data.degree-1]))
				if (edge_data.degree>min_valid_column-2) and \
					edge_data.value>=cor_cut_off_vector[edge_data.degree-1] and\
					edge_data.value<=1.0:	#significant enough
						writer.writerow([repr(edge_tuple1), repr(edge_tuple2), edge_data.value])
		
		sys.stderr.write("Done.\n")
		del writer	
		
	def run(self):
		"""
		05-26-05
			--vertex_set_readin()
			--cor_vector_readin()
			--graph_modeling.cor_cut_off_vector_construct()
			--construct_2nd_graph()
		"""
		vertex_set = self.vertex_set_readin(self.input_file, self.which_line)
		edge_tuple2cor_vector = self.cor_vector_readin(self.cor_vector_file, vertex_set)
		
		cor_cut_off_vector = graph_modeling.cor_cut_off_vector_return(self.p_value_cut_off, self.cor_cut_off)	#construct the cor_cut_off vector
		self.construct_2nd_graph(self.output_file, edge_tuple2cor_vector, self.p_value_cut_off, \
			self.min_valid_column, cor_cut_off_vector)


	
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hi:w:v:s:o:b:p:c:ru", ["help",  "input_file=", \
			"which_line=", "cor_vector_file=", "sig_vector_file=", "output_file=", \
			"label=", "p_value_cut_off=", "cor_cut_off=", "report", "debug"])
	except:
		print __doc__
		sys.exit(2)

	input_file = None
	which_line = 0
	cor_vector_file = None
	sig_vector_file = None
	output_file = None
	label = None
	p_value_cut_off = 0.01
	cor_cut_off = 0.6
	report = 0
	debug = 0
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-i", "--input_file"):
			input_file = arg
		elif opt in ("-w", "--which_line"):
			which_line = int(arg)
		elif opt in ("-v", "--cor_vector_file"):
			cor_vector_file = arg
		elif opt in ("-s", "--sig_vector_file"):
			sig_vector_file = arg
		elif opt in ("-o", "--output_file"):
			output_file = arg
		elif opt in ("-b", "--label"):
			label = int(arg)
		elif opt in ("-p", "--p_value_cut_off"):
			p_value_cut_off = float(arg)
		elif opt in ("-c", "--cor_cut_off"):
			cor_cut_off = float(arg)
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-u", "--debug"):
			debug = 1
			

	if input_file and output_file and cor_vector_file:
		instance =check_netmine2nd(input_file, which_line, cor_vector_file,\
			sig_vector_file, output_file, label, p_value_cut_off, \
			cor_cut_off, report, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
