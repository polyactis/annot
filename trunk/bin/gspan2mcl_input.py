#!/usr/bin/env python
"""
Usage: gspan2mcl_input.py [OPTION] INPUTFILE OUTPUTFILE

Option:
	INPUTFILE is in gspan format.
	OUTPUTFILE will be in mcl_input format.
	-w, --weight	add edge weight to the mcl matrix
	-h, --help	show this help
	
Examples:
	gspan2mcl_input.py gph_result/sc/gph_dataset1 gph_result/sc_mcl/mcl_gph_dataset1

Description:
	This program converts gspan formatted file generated by graph_reorganize.py
	to file in mcl_input format.
	Use batch.qsub.py to iterate over files in a directory.
"""


import sys, os, re, getopt
from graphlib import Graph
class gspan2mcl_input:
	'''
	'''
	def __init__(self, infname, ofname, weight):
		self.infname = infname
		self.inf = open(self.infname, 'r')
		self.of = open(ofname, 'w')
		self.weight = int(weight)
		self.graph = Graph.Graph()
		#a number extracted from filename used as a pseudo splat_id
		self.no = 0
	
	def dstruc_loadin(self):
		sys.stderr.write("Loading Data STructure...")
		p_no = re.compile(r'\d+$')
		try:
			self.no = int(p_no.search(self.infname).group())	#integer conversion
		except AttributeError, error:
			sys.stderr.write('%s\n'%error)
			sys.stderr.write("can't find the number of this dataset from %s\n"%self.infname)
			sys.exit(2)
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
	
	def output(self):
		vertex_list = self.graph.node_list()
		dim = len(vertex_list)
		out_block = '(splat_id %s )\n'%self.no	#here it is '=' not '+='
		out_block += '(mclheader\n'
		out_block += 'mcltype matrix\n'
		out_block += 'dimensions %dx%d\n)\n'%(dim,dim)
		out_block += '(mcldoms\n'
		vertex_list.sort()
		for vertex in vertex_list:
			out_block += '%s '%vertex
		out_block += '$\n)\n'
		out_block += '(mclmatrix\nbegin\n'
		self.of.write(out_block)
		for vertex in vertex_list:
			self.of.write('%s '%vertex)
			out_neighbors = self.graph.out_nbrs(vertex)
			for neighbor in out_neighbors:
				if self.weight:
					edge_id = self.graph.edge_by_node(vertex,neighbor)
					edge_data = self.graph.edge_data(edge_id)
					self.of.write('%s:%s '%(neighbor, edge_data))
				else:
					self.of.write('%s '%neighbor)
			self.of.write('$\n')
		self.of.write(')\n\n')
	
	def run(self):
		self.dstruc_loadin()
		self.output()
	
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "wh", ["weight", "help"])
	except:
		print __doc__
		sys.exit(2)
	
	weight = 0
	for opt, arg in opts:
		if opt in ("-w", "--weight"):
			weight = 1
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)

			
	if len(args) == 2:
		instance = gspan2mcl_input(args[0], args[1], weight)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
