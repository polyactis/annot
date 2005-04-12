#!/usr/bin/env python
"""
Usage: EdgeClusterFromCopathOutput.py [OPTION] INPUTFILE OUTPUTFILE

Option:
	
	
Examples:
	EdgeClusterFromCopathOutput.py Fsc54_5G1E6D40Q40S200C50H4J40W40Z0010
		Fsc54_5G1E6D40Q40S200C50H4J40W40Z0010E
	
Description:
	This program restores the edge clusters(2nd-order cluster) from copath's 1st-order cluster
	output.
	The result is got by running copath with min_graph_size = 1. So no connected component
	is lost. The program doesn't restore the graph between the edges.

"""

import sys, os, getopt, csv, re
from codense.codense2db import cluster_dstructure

class EdgeClusterFromCopathOutput:
	"""
	04-12-05
	"""
	def __init__(self, infname, outfname):
		self.infname = infname
		self.outfname = outfname
		
		#the main data repositary
		self.cooccurrent_cluster_id2cluster = {}
		#used to find the cooccurrent_cluster_id, like '6.6.9', cooccurrent_cluster_id is '6.6'
		self.p_cooccurrent_cluster_id = re.compile(r'\d+\.\d+')

	def copath_parser(self, row, argument=None, argument2=None):
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
		
		if cooccurrent_cluster_id not in self.cooccurrent_cluster_id2cluster:
			cluster = cluster_dstructure()
			cluster.cooccurrent_cluster_id = cooccurrent_cluster_id	#it's not used in the output()
			#initialize two sets
			cluster.vertex_set = []
			cluster.edge_set = []

			self.cooccurrent_cluster_id2cluster[cooccurrent_cluster_id] = cluster
		
		#pass it to ease programming
		cluster = self.cooccurrent_cluster_id2cluster[cooccurrent_cluster_id]
		cluster.vertex_set += vertex_set
		cluster.edge_set += edge_set

	def output(self, outfname, cooccurrent_cluster_id2cluster):
		"""
		04-12-05
			output the cooccurrent_cluster_id2cluster in the same format as the INPUTFILE.
			
			--return_vertex_set_string
			--return_edge_set_string
		"""
		outf = open(outfname, 'w')
		writer = csv.writer(outf, delimiter='\t')
		sys.stderr.write("Outputting the edge clusters...")
		for cooccurrent_cluster_id,cluster in cooccurrent_cluster_id2cluster.iteritems():
			#ascending order
			cluster.vertex_set.sort()
			cluster.edge_set.sort()
			
			no_of_nodes = len(cluster.vertex_set)
			cluster.splat_connectivity = 2.0*len(cluster.edge_set)/(no_of_nodes*(no_of_nodes-1))	#2.0 keeps the whole result float
			row = [cooccurrent_cluster_id, cluster.splat_connectivity, self.return_vertex_set_string(cluster.vertex_set), \
				self.return_edge_set_string(cluster.edge_set)]
			writer.writerow(row)
		del writer
		outf.close()
		sys.stderr.write("Done.\n")
	
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
		for row in inf:
			self.copath_parser(row)
		del inf
		
		self.output(self.outfname, self.cooccurrent_cluster_id2cluster)
		
		
if __name__ == '__main__':
	if len(sys.argv) != 3:
		print __doc__
		sys.exit(2)
	else:
		instance = EdgeClusterFromCopathOutput(sys.argv[1], sys.argv[2])
		instance.run()
