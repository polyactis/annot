#!/usr/bin/env python
"""
Usage: GO_common_ancestor.py -l LIST_FILE [OPTION] Outputfile

Option:
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database, go(default)
	-g ..., --organism=...	two letter organism abbreviation, sc(default)
	-t ..., --type=...	which branch of GO, 0, 1(default) or 2
	-l ..., --list_file=...	the file contains a list of go id's
	-h, --help              show this help

Examples:
	GO_common_ancestor.py -d GO -l hs_informative.stat -g sc go.xml
	GO_common_ancestor.py -l hs_informative.stat -g hs go.xml
	GO_common_ancestor.py -d GO -k go -l hs_yh60_42.stat -t 0 molecular_function.xml

Description:
	This class builds a big GO graph from GO database(two tables, term and term2term).
	And given a list of GO terms, it outputs a common ancestor DAG
	in GraphXML format, which is readable by hypergraph to visualize.
	Organism information is used to show the associated genes.
	type illustration:
	0:	molecular_function
	1:	biological_process
	2:	cellular_component
"""

import sys, os, psycopg, getopt, csv
from GO_graphxml import GO_graphxml
from graphlib import Graph
from sets import Set
from kjbuckets import *

class GO_common_ancestor(GO_graphxml):
	'''
	'''
	def __init__(self, dbname, schema, type, output, orgn, go_set_fname):
		GO_graphxml.__init__(self, dbname, schema, type, output, orgn)
		self.go_set_file = csv.reader(file(go_set_fname), delimiter='\t')
		self.go_id_set = Set()

	def dstruc_loadin(self):
		GO_graphxml.dstruc_loadin(self)
		for row in self.go_set_file:
			self.go_id_set.add(row[0])
		
	def run(self):
		self.dstruc_loadin()
		kjgraph = kjGraph()
		for go_id in self.go_id_set:
			if go_id in self.acc2id_dict:
				go_node = self.acc2id_dict[go_id]
			#term doesn't exist
			else:
				sys.stderr.write('%s inexists\n'%go_id)
				continue
			#this termid is obsolete
			if self.termid_dict[go_node].is_obsolete == 1:
				sys.stderr.write("%s is obsolete\n"%go_id)
				continue
			subgraph = self.go_graph.back_bfs_subgraph(go_node)
			kjgraph += self.graph2kjgraph(subgraph)
		graph = self.kjgraph2graph(kjgraph)
		self.subgraph_to_graphxml(graph)

	def graph2kjgraph(self, graph):
		kjgraph = kjGraph()
		for edge in graph.edges.values():
			kjgraph.add((edge[0], edge[1]))
		return kjgraph

	def kjgraph2graph(self, kjgraph):
		graph = Graph.Graph()
		for item in kjgraph.items():
			graph.add_edge(item[0], item[1])
		return graph

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)

	long_options_list = ["help", "dbname=", "schema=", "type=", "organism=", "list_file="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hd:k:t:g:l:", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	dbname = 'graphdb'
	schema = 'go'
	organism = 'sc'
	type = 1
	list_file = ''
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-k", "--schema"):
			schema = arg
		elif opt in ("-t", "--type"):
			type = int(arg)
		elif opt in ("-g", "--organism"):
			organism = arg
		elif opt in ("-l", "--list_file"):
			list_file = arg

	if len(args) == 1 and list_file:
		instance = GO_common_ancestor(dbname, schema, type, args[0], organism, list_file)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
