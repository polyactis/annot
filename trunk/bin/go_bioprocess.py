#!/usr/bin/env python
"""
Usage: go_bioprocess.py -k SCHEMA -p PARSER -u UNKNOWNFILE [OPTION] go_function_file

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-c --commit	commits the database transaction
	-u ..., --unknown=...	the file containing the unknown gene_ids(seperated by ';')
	-p ..., --parser=...	which parser to use
	-g ..., --organism=...	IGNORE it, interface relics
	-h, --help              show this help
	
Examples:
	go_bioprocess.py -k shu -p min -u yeast_unknown yestprocess2.txt

Description:
	This program extracts go functional class information from raw file,
	combines the stable unknown gene set and sets up the schema.go
	table in the database.
	It depends on schema.gene.
"""

import pickle,sys,os, psycopg, getopt, csv
from graphlib import Graph, GraphAlgo

														
class go_term:
	'''
	a structure for holding information related to a GeneOntology term.
	'''
	def __init__(self):
		no = None
		name = None
		no_of_genes = None
		whole_gene_array = None
		gene_array = None

class ming_parser:
	def __init__(self):
		self.go_dict = {}
	
	def parse(self, inf, vertex_dict):
		reader = csv.reader(inf, delimiter='\t')
		for row in reader:
			go_id = row[0]
			info = go_term()
			info.name = row[1]
			info.whole_gene_array = []
			for gene in row[2:]:
				if gene != '':
					info.whole_gene_array.append(gene)
			info.no_of_genes = len(info.whole_gene_array)
			info.gene_array = []
			for gene in info.whole_gene_array:
				if gene in vertex_dict:
					info.gene_array.append(vertex_dict[gene])
			self.go_dict[go_id] = info
			
		key_list = self.go_dict.keys()
		key_list.sort()
		for i in range(len(key_list)):
			id = key_list[i]
			self.go_dict[id].no = i+1
		
		return self.go_dict

class min_parser:
	def __init__(self):
		self.go_dict = {}
	
	def parse(self, inf, vertex_dict):
		reader = csv.reader(inf, delimiter='\t')
		for row in reader:
			go_id = row[0]
			if go_id not in self.go_dict:
				info = go_term()
				info.name = row[1]
				info.whole_gene_array = []
				self.go_dict[go_id] = info
			ls = row[2].split('|')
			for item in ls:
				self.go_dict[go_id].whole_gene_array.append(item)
		
		for go_id in self.go_dict:
			entry = self.go_dict[go_id]
			entry.no_of_genes = len(entry.whole_gene_array)
			entry.gene_array = []
			for gene in entry.whole_gene_array:
				if gene in vertex_dict:
					entry.gene_array.append(vertex_dict[gene])
			
		key_list = self.go_dict.keys()
		key_list.sort()
		for i in range(len(key_list)):
			id = key_list[i]
			self.go_dict[id].no = i+1
		
		return self.go_dict


parser_map = {"ming":ming_parser(),\
			"min":min_parser()}

class go_table_setup:
	def __init__(self, fname, hostname, dbname, schema, parser, u_fname, orgn, needcommit=0):
		self.go_inf = open(fname, 'r')
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.parser = parser_map[parser]
		self.reader = csv.reader(open(u_fname, 'r'), delimiter=';')
		self.needcommit = int(needcommit)
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
		#mapping between gene_id and gene_no
		self.vertex_dict = {}
		#stores all the unknown genes in schema.gene
		self.unknown_gene_list = []
		#GO DAG
		self.go_graph = Graph.Graph()

	def dstruc_loadin(self):
		sys.stderr.write("Loading Data STructure...")

		#sets up the self.vertex_dict
		self.curs.execute("select gene_id, gene_no from gene")
		rows = self.curs.fetchall()
		for row in rows:
			self.vertex_dict[row[0]] = row[1]
		#sets up the self.unknown_gene_list
		self._unknown_gene_list = self.reader.next()
		for gene in self._unknown_gene_list:
			if gene in self.vertex_dict:
				self.unknown_gene_list.append(self.vertex_dict[gene])
			
		#get the non-obsolete biological_process GO DAG
		self.curs.execute("select t2t.term1_id, t2t.term2_id, t1.acc, t2.acc from \
			go.term2term t2t, go.term t1, go.term t2 where t2t.term1_id=t1.id and \
			t2t.term2_id=t2.id and t1.is_obsolete=0 and t2.is_obsolete=0 and \
			t1.term_type='biological_process' and t2.term_type='biological_process' ")
		rows = self.curs.fetchall()
		for row in rows:
		#setup the go_graph structure
			self.go_graph.add_edge(row[2], row[3])

		sys.stderr.write("Done\n")

	def submit(self):
		sys.stderr.write("Database transacting...")
		#following string operations are because of format restrictions of database array input
		string__unknown_gene_list = repr(self._unknown_gene_list)
		string__unknown_gene_list = string__unknown_gene_list.replace("'", '')
		string__unknown_gene_list = '{' + string__unknown_gene_list[1:-1] + '}'
		string_unknown_gene_list = repr(self.unknown_gene_list)
		string_unknown_gene_list = '{' + string_unknown_gene_list[1:-1] + '}'
		
		self.curs.execute("insert into go(go_id, go_no, no_of_genes, name, whole_gene_array, gene_array, depth) \
			values('%s', %d, %d, '%s', '%s', '%s', 2)"%('GO:0000004', 0, len(self._unknown_gene_list), \
			'biological_process unknown', string__unknown_gene_list, string_unknown_gene_list))
		go_dict = self.parser.parse(self.go_inf, self.vertex_dict)
		for term in go_dict:
			string_whole_gene_array = repr(go_dict[term].whole_gene_array)
			string_gene_array = repr(go_dict[term].gene_array)
			string_whole_gene_array = string_whole_gene_array.replace("'", '')
			string_whole_gene_array = '{' + string_whole_gene_array[1:-1] + '}'
			string_gene_array = '{' + string_gene_array[1:-1] + '}'
			go_dict[term].name = go_dict[term].name.replace("'",'')
			depth = len(GraphAlgo.shortest_path(self.go_graph, 'GO:0008150', term))
			self.curs.execute("insert into go(go_id, go_no, no_of_genes, name, whole_gene_array, gene_array, depth) \
				values('%s', %d, %d, '%s', '%s', '%s', %d)"%(term, go_dict[term].no, go_dict[term].no_of_genes,\
				go_dict[term].name, string_whole_gene_array, string_gene_array, depth))
		if self.needcommit:
			self.conn.commit()
		sys.stderr.write("done.\n")
			
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:cu:p:g:", ["help", "hostname", "dbname=", "schema=", "commit", "unknown=","parser=","organism="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	commit = 0
	unknown = ''
	parser = ''
	organism = 'sc'
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
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-u", "-unknown"):
			unknown = arg
		elif opt in ("-p", "-parser"):
			parser = arg
		elif opt in ("-g", "-organism"):
			organism = arg
			
	if schema and unknown and parser and len(args)>0:
		instance = go_table_setup(args[0], hostname, dbname, schema, parser, unknown, organism, commit)
		instance.dstruc_loadin()
		instance.submit()
	else:
		print __doc__
		sys.exit(2)
