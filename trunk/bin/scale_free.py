#!/usr/bin/env python
"""
Usage: scale_free.py [OPTION] INPUTFILE

Option:

Examples:
	scale_free.py gph_result/gph_dataset1
	
Description:
	This program reads a graph file in gspan format and plots
	the degree_distribution log-log curve.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import sys, os, getopt
from graphlib import Graph
from rpy import r

class scale_free:
	def __init__(self, infname):
		self.inf = open(infname, 'r')
		self.graph = Graph.Graph()
		self.degree_dict = {}
	
	def dstruc_loadin(self):
		sys.stderr.write("Loading Data STructure...")
		for line in self.inf:
			if line[0] == 'e':
				#edge here, like 'e 3807 3859 0.804645'
				line_list = line[:-1].split()
				vertex1 = int(line_list[1])
				vertex2 = int(line_list[2])
				edge_data = line_list[3]
				#graphlib is a data structure to store directional graph
				self.graph.add_edge(vertex1, vertex2, edge_data)
				self.graph.add_edge(vertex2, vertex1, edge_data)
		sys.stderr.write("Done\n")
	
	def plot(self):
		vertex_list = self.graph.node_list()
		number_of_nodes = len(vertex_list)
		for vertex in vertex_list:
			degree = self.graph.inc_degree(vertex) + self.graph.out_degree(vertex)
			if degree not in self.degree_dict:
				self.degree_dict[degree] = 1
			else:
				self.degree_dict[degree] += 1
		r.pdf('degree_distribution.pdf')
		x_list = []
		y_list = []
		for degree in self.degree_dict:
			x_list.append(r.log(degree))
			y_list.append(r.log(float(self.degree_dict[degree])/number_of_nodes))
		r.plot(x_list, y_list, type='p', xlab='log(k)', ylab='log(p(k))')
		r.dev_off()
	
	def run(self):
		self.dstruc_loadin()
		self.plot()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "h", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
			
	if len(args)==1:
		instance = scale_free(args[0])
		instance.run()
	else:
		print __doc__
		sys.exit(2)
