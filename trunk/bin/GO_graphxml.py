#!/usr/bin/env python
"""
Usage: GO_graphxml.py [OPTION] Outputfile

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database, go(default)
	-g ..., --organism=...	two letter organism abbreviation, sc(default)
	-t ..., --type=...	which branch of GO, 0, 1(default) or 2
	-h, --help              show this help

Examples:
	GO_graphxml.py -d GO -g sc go.xml
	GO_graphxml.py -d graphdb -g hs go.xml
	GO_graphxml.py -d GO -k go -t 0 molecular_function.xml

Description:
	This class builds a big GO graph from GO database(two tables, term and term2term).
	And given a GO term, it can output either its parent subgraph or child subgraph
	in GraphXML format, which is readable by hypergraph to visualize.
	Organism information is used to show the associated genes.
	type illustration:
	0:	molecular_function
	1:	biological_process
	2:	cellular_component
"""

import sys, os, psycopg, getopt
from xml.dom.minidom import getDOMImplementation
from graphlib import Graph
from sets import Set

class termid_attr:
	def __init__(self, acc, name, is_obsolete):
		self.acc = acc
		self.name = name
		self.is_obsolete = is_obsolete
		self.individual = Set()
		self.family = Set()
		
class GO_graphxml:
	'''
	This class builds a big GO graph from GO database(two tables, term and term2term).
	And given a GO term, it can output either its parent subgraph or child subgraph
	in GraphXML format, which is readable by hypergraph to visualize.
	'''
	#orgn ='sc' is for interface compliance. as this class is inherited by Go_no_parent_filter.
	def __init__(self, hostname, dbname, schema, type, output, orgn='sc'):
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		
		self.type_dict = {0:'molecular_function', 1:'biological_process', 2:'cellular_component'}
		self.type = self.type_dict[int(type)]
		self.ofname = output
		self.org_short2long = {'at':'Arabidopsis thaliana',
			'ce':'Caenorhabditis elegans',
			'dm':'Drosophila melanogaster',
			'hs':'Homo sapiens',
			'mm':'Mus musculus',
			'sc':'Saccharomyces cerevisiae',
			'Arabidopsis thaliana':'Arabidopsis thaliana',
			'Caenorhabditis elegans':'Caenorhabditis elegans',
			'Drosophila melanogaster':'Drosophila melanogaster',
			'Homo sapiens':'Homo sapiens',
			'Mus musculus':'Mus musculus',
			'Gorilla gorilla Pan paniscus Homo sapiens':'Homo sapiens',
			'Saccharomyces cerevisiae':'Saccharomyces cerevisiae'}
		self.organism = self.org_short2long[orgn]
		self.log_file = open('/tmp/GO_graphxml.log', 'w')		
		self.termid_dict = {}
		self.acc2id_dict = {}
		self.go_graph = Graph.Graph()
		
	def dstruc_loadin(self):
		'''
		Construcs three data structures from database,
		self.termid_dict
		self.acc2id_dict
		self.go_graph
		'''
		sys.stderr.write("Loading Data STructure...")
		self.curs.execute("select id, acc, name, is_obsolete from term")
		rows = self.curs.fetchall()
		for row in rows:
			self.termid_dict[row[0]] = termid_attr(row[1], row[2], row[3])
			self.acc2id_dict[row[1]] = row[0]
		
		#construct the graph
		self.curs.execute("select term1_id, term2_id, relationship_type_id from term2term")
		rows = self.curs.fetchall()
		for row in rows:
			#we don't want edge to contain obsolete term
			#if self.termid_dict[row[0]].is_obsolete == 0 and self.termid_dict[row[1]].is_obsolete ==0:
			self.go_graph.add_edge(row[0], row[1])
		
		#for each termid, loads in its associated genes
		self.curs.execute("select t.id, a.gene_id from graph.association a, graph.gene_id_to_no g, go.term t\
			where a.gene_id = g.gene_id and a.go_id=t.acc and g.organism='%s'"%self.organism)
		rows = self.curs.fetchall()
		for row in rows:
			self.termid_dict[row[0]].individual.add(row[1])
		
		#collect the associated genes from descendent terms		
		for termid in self.termid_dict:
			#assign its own associated genes
			self.termid_dict[termid].family |= self.termid_dict[termid].individual
			#from descendents
			if termid not in self.go_graph.nodes:
				continue
			descendents = self.go_graph.forw_bfs(termid)
			for descendent in descendents:
				#sets.Set is different from kjSet in the union operation.
				self.termid_dict[termid].family |= self.termid_dict[descendent].individual
		
		sys.stderr.write("Done\n")
		
	def subgraph_to_graphxml(self, subgraph):
		'''
		outputs subgraph(a graphlib.Graph instance) in GraphXML format.
		'''
		of = open(self.ofname, 'w')
		impl = getDOMImplementation('')
		#Create a document type node using the doctype name "GraphXML"
		#system ID is GraphXML.dtd, and blank public ID
		doctype = impl.createDocumentType(u"GraphXML", '', 'GraphXML.dtd' )
		doc = impl.createDocument(None, 'GraphXML', doctype)
		#'GraphXML' is the top_element
		top_element = doc.documentElement
		graph = doc.createElement('graph')
		top_element.appendChild(graph)
		
		for go_node in subgraph.node_list():
			#node is child of graph
			#label and dataref are childs of node
			node = doc.createElement('node')
			node.setAttribute('name', '%d'%go_node)
			
			label = doc.createElement('label')
			text = doc.createTextNode("%s %d(%d)"%(self.termid_dict[go_node].acc,\
				len(self.termid_dict[go_node].individual), len(self.termid_dict[go_node].family)) )
			label.appendChild(text)
			
			dataref = doc.createElement('dataref')
			ref = doc.createElement('ref')
			ref.setAttribute("xlink:href", self.termid_dict[go_node].name)
			dataref.appendChild(ref)
			
			node.appendChild(label)
			node.appendChild(dataref)
			graph.appendChild(node)
		
		for edge_id in subgraph.edge_list():
			#edge is child of graph
			go_edge = subgraph.edge_by_id(edge_id)
			edge = doc.createElement('edge')
			edge.setAttribute('source', '%d'%go_edge[0])
			edge.setAttribute('target', '%d'%go_edge[1])
			graph.appendChild(edge)
		
		doc.writexml(of)
		of.close()
	
	def cgi_run(self, go_acc, backward):
		'''
		For cgi purpose, given a go accession number and direction information
		'''
		self.dstruc_loadin()
		
		if go_acc in self.acc2id_dict:
			go_node = self.acc2id_dict[go_acc]
		#term doesn't exist
		else:
			print 'ERROR: %s inexists\n'%go_acc
			return
		#this termid is obsolete
		if self.termid_dict[go_node].is_obsolete == 1:
			print "ERROR: %s is obsolete\n"%go_acc
			return
		if backward:
			subgraph = self.go_graph.back_bfs_subgraph(go_node)
			self.subgraph_to_graphxml(subgraph)
		else:
			subgraph = self.go_graph.forw_bfs_subgraph(go_node)
			self.subgraph_to_graphxml(subgraph)

	def run(self):
		'''
		given a GO term, find its parent subgraph(back_bfs_subgraph)
		or child subgraph(forw_bfs_subgraph)
		'''
		self.dstruc_loadin()
		while 1:
			go_id = raw_input("Input a GO id:\t")
			#program exits.
			if go_id == 'q':
				break
			elif go_id in self.acc2id_dict:
				go_node = self.acc2id_dict[go_id]
			#term doesn't exist
			else:
				sys.stderr.write('%s inexists\n'%go_id)
				continue
			#this termid is obsolete
			if self.termid_dict[go_node].is_obsolete == 1:
				sys.stderr.write("%s is obsolete\n"%go_id)
				continue
			#output the genes associated with this go_id.
			for gene_id in self.termid_dict[go_node].family:
				self.log_file.write('%s\t%s\n'%(gene_id, go_id))

			direction = raw_input("backward(b) or forward(f):\t")
			if direction == 'b':
				subgraph = self.go_graph.back_bfs_subgraph(go_node)
				self.subgraph_to_graphxml(subgraph)
			elif direction == 'f':
				subgraph = self.go_graph.forw_bfs_subgraph(go_node)
				self.subgraph_to_graphxml(subgraph)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)

	long_options_list = ["help", "hostname=", "dbname=", "schema=", "type=", "organism="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:g:", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = 'go'
	organism = 'sc'
	type = 1
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
		elif opt in ("-t", "--type"):
			type = int(arg)
		elif opt in ("-g", "--organism"):
			organism = arg

	if len(args) == 1:
		instance = GO_graphxml(hostname, dbname, schema, type, args[0], organism)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
