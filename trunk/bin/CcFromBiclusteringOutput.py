#!/usr/bin/env python
"""
Usage: CcFromBiclusteringOutput.py [OPTION] INPUTFILE OUTPUTFILE

Option:
	
	
Examples:
	CcFromBiclusteringOutput.py Fsc54_5G1E6D40Q40S200C50H4J40W40Z0010
		Fsc54_5G1E6D40Q40S200C50H4J40W40Z0010E
	
Description:
	For each biclustering output, crack the edge_list into several connected components.
	And regard them as individual cluster.

"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import getopt, csv, re
from codense.codense2db import cluster_dstructure
from graph.cc_from_edge_list import cc_from_edge_list
from sets import Set

class CcFromBiclusteringOutput:
	"""
	04-12-05
		most of it copied from EdgeClusterFromCopathOutput.py
	"""
	def __init__(self, infname=None, outfname=None):
		self.infname = infname
		self.outfname = outfname
		
		#the main data repositary
		self.cooccurrent_cluster_id2cluster = {}
		#used to find the cooccurrent_cluster_id, like '6.6.9', cooccurrent_cluster_id is '6.6'
		self.p_cooccurrent_cluster_id = re.compile(r'\d+\.\d+')

	def copath_parser(self, row, writer, argument=None, argument2=None):
		"""
		04-12-05
			copied from codense2db.py, changed a lot
		"""
	
		cooccurrent_cluster_id = self.p_cooccurrent_cluster_id.match(row[0]).group()
		vertex_set = row[2][1:-2].split(';')
		vertex_set = map(int, vertex_set)
		edge_list = row[3][2:-4].split(' );(')
		edge_set = []
		for edge in edge_list:
			edge = edge.split(',')
			edge = map(int, edge)
			#in ascending order
			edge.sort()
			edge_set.append(edge)
		#04-29-05	cc module come into play to get the connected components
		instance = cc_from_edge_list()
		instance.run(edge_set)
		cc_list = instance.cc_list
		for cc_edge_list in cc_list:
			cluster = cluster_dstructure()
			cluster.cooccurrent_cluster_id = cooccurrent_cluster_id	#it's not used in the output()
			#initialize two sets
			cluster.vertex_set = self.vertex_set_from_cc_edge_list(cc_edge_list)
			cluster.edge_set = cc_edge_list
			self.output(writer, cluster)
		

	def vertex_set_from_cc_edge_list(self, cc_edge_list):
		"""
		04-29-05
		"""
		vertex_set = Set()
		for edge in cc_edge_list:
			vertex_set.add(edge[0])
			vertex_set.add(edge[1])
		return list(vertex_set)
	
	def output(self, writer, cluster):
		"""
		04-12-05
			output the cooccurrent_cluster_id2cluster in the same format as the INPUTFILE.
			
			--return_vertex_set_string
			--return_edge_set_string
		"""

		#ascending order
		cluster.vertex_set.sort()
		cluster.edge_set.sort()
		
		no_of_nodes = len(cluster.vertex_set)
		cluster.splat_connectivity = 2.0*len(cluster.edge_set)/(no_of_nodes*(no_of_nodes-1))	#2.0 keeps the whole result float
		row = [cluster.cooccurrent_cluster_id, cluster.splat_connectivity, self.return_vertex_set_string(cluster.vertex_set), \
			self.return_edge_set_string(cluster.edge_set)]
		writer.writerow(row)

	
	def return_vertex_set_string(self, vertex_set):
		"""
		04-12-05
			convert the list(vertex_set) into a string form as in copath's result
		"""
		vertex_set = map(repr, vertex_set)
		vertex_set_string = ';'.join(vertex_set)
		vertex_set_string = '{'+vertex_set_string+';}'
		return vertex_set_string
	
	def return_edge_set_string(self, edge_set):
		"""
		04-12-05
			convert the 2d list (edge_set) into a string form as in copath's result
		"""
		edge_set_string = ''
		for edge in edge_set:
			edge_set_string += '(%s,%s );'%(edge[0],edge[1])
		edge_set_string = '{' + edge_set_string + '}'
		return edge_set_string
	
	def run(self):
		"""
		04-12-05
			
			(loop over the INPUTFILE)
				--copath_parser()
			
			--output()
				--return_vertex_set_string()
				--return_edge_set_string()
		"""

		inf = csv.reader(open(self.infname, 'r'), delimiter='\t')
		outf = open(self.outfname, 'w')
		writer = csv.writer(outf, delimiter='\t')
		for row in inf:
			self.copath_parser(row, writer)
		del inf
		del writer
		outf.close()
		
		
if __name__ == '__main__':
	if len(sys.argv) != 3:
		print __doc__
		sys.exit(2)
	else:
		instance = CcFromBiclusteringOutput(sys.argv[1], sys.argv[2])
		instance.run()
