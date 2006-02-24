#!/usr/bin/env python
"""
Usage: sample_dataset_and_graph.py [OPTION] -i INPUTFILE -o OUTPUTDIR

Option:
	-i ...,	inputfile
	-o ...,	directory containing output files
	-p ...,	percentage of columns to be sampled each time, 0.6(default)
	-m ...,	minimum number of columns for one sampling, 8(default)
	-b,	enable debugging
	-h, --help              show this help
	
Examples:
	

Description:
	sample some columns out from a dataset, compute graphs, check the difference
	between graphs. If OUTPUTDIR doesn't exist, it'll create one.
	
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import getopt, csv, math, random

class sample_dataset_and_graph:
	def __init__(self, inputfile, outputdir, percentage=0.6, min_no_of_columns=8, debug=0):
		""" 
		02-23-06
		"""
		self.inputfile = inputfile
		self.outputdir = outputdir
		self.percentage = float(percentage)
		self.min_no_of_columns = int(min_no_of_columns)
		self.debug = int(debug)
	
	def get_no_of_columns(self, filename, delimiter_char='\t'):
		reader = csv.reader(open(filename, 'r'), delimiter=delimiter_char)
		no_of_columns = 0
		for row in reader:
			no_of_columns = len(row)-1	#1st column, gene id is not included
			break
		del reader
		print "%s has %s columns"%(filename, no_of_columns)
		return no_of_columns
	
	def sample_graph(self, inputfile, no_of_columns, percentage, min_no_of_columns, outputdir, sample_order):
		no_of_sample_columns = int(max(no_of_columns*percentage, min_no_of_columns))
		sample_column_no_list = random.sample(range(1, no_of_columns+1), no_of_sample_columns)
		output_fname = os.path.join(outputdir, '%s.%s.%s'%(os.path.basename(inputfile), no_of_sample_columns, sample_order))
		gph_output_fname = '%s.o'%output_fname
		print "%s sample_column_no_list: %s"%(os.path.basename(output_fname), sample_column_no_list)
		sample_column_no_list = map(str, sample_column_no_list)
		job_fname = 'cut -f 1,%s %s > %s'%(','.join(sample_column_no_list), inputfile, output_fname)
		wl = ['sh', '-c', job_fname]
		os.spawnvp(os.P_WAIT, 'sh', wl)
		os.system('~/script/annot/bin/graph/graph_modeling -o %s -p 0 -c 0 %s'%(gph_output_fname, output_fname))
		
		return gph_output_fname
		
	def compare_graphs(self, seed_gph_fname, other_gph_fname_list):
		sys.stderr.write("Reading seed graph...")
		from sets import Set
		import csv
		reader = csv.reader(open(seed_gph_fname, 'r'), delimiter='\t')
		edge_dict = {}
		vertex_set = Set()
		for row in reader:
			if row[0] == 'e':
				edge_tuple = (int(row[1]), int(row[2]))
				if edge_tuple[0]>edge_tuple[1]:
					edge_tuple = (edge_tuple[1], edge_tuple[0])
				edge_dict[edge_tuple] = 1
				vertex_set.add(edge_tuple[0])
				vertex_set.add(edge_tuple[1])
		del reader
		sys.stderr.write("done.\n")
		
		sys.stderr.write("start to compare with other graphs...\n")
		for input_fname in other_gph_fname_list:
			print '\t' + input_fname
			reader = csv.reader(open(input_fname, 'r'), delimiter='\t')
			no_of_overlapping_edges = 0
			no_of_total_edges = 0
			local_vertex_set = Set()
			for row in reader:
				if row[0] == 'e':
					edge_tuple = (int(row[1]), int(row[2]))
					if edge_tuple[0]>edge_tuple[1]:
						edge_tuple = (edge_tuple[1], edge_tuple[0])
					if edge_tuple in edge_dict:
						edge_dict[edge_tuple] += 1
						no_of_overlapping_edges += 1
					no_of_total_edges += 1
					local_vertex_set.add(edge_tuple[0])
					local_vertex_set.add(edge_tuple[1])
			vertex_overlapping_ratio = len(local_vertex_set&vertex_set)/float(len(local_vertex_set|vertex_set))
			edge_overlapping_ratio = float(no_of_overlapping_edges)/(len(edge_dict) + no_of_total_edges - no_of_overlapping_edges)
			print '\t overlapping ratio: %s(edge), %s(vertex)'%(edge_overlapping_ratio, vertex_overlapping_ratio)
			del reader
		sys.stderr.write("done.\n")
		
		return edge_dict
	
	def run(self):
		no_of_columns = self.get_no_of_columns(self.inputfile)
		gph_output_fname_list = []
		if not os.path.isdir(self.outputdir):
			os.makedirs(self.outputdir)
		for i in range(10):
			gph_output_fname = self.sample_graph(self.inputfile, no_of_columns, self.percentage, self.min_no_of_columns, self.outputdir, i)
			gph_output_fname_list.append(gph_output_fname)
		edge_dict = self.compare_graphs(gph_output_fname_list[0], gph_output_fname_list[1:])

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hi:o:p:m:b", ["help"])
	except:
		print __doc__
		sys.exit(2)
	
	inputfile = None
	outputdir = None
	percentage = 0.6
	min_no_of_columns = 8
	debug = 0
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-i",):
			inputfile = arg
		elif opt in ("-o",):
			outputdir = arg
		elif opt in ("-p",):
			percentage = float(arg)
		elif opt in ("-m",):
			min_no_of_columns = int(arg)
		elif opt in ("-b",):
			debug = 1
			
	if inputfile and outputdir:
		instance = sample_dataset_and_graph(inputfile, outputdir, percentage, min_no_of_columns, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
