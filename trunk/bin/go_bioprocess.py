#!/usr/bin/env python
"""
Usage: go_bioprocess.py -d DATABASENAME -k SCHEMA -p PARSER -u UNKNOWNFILE -g ORGANISM [OPTION] go_function_file

Option:
	-d ..., --dbname=...	the database name
	-k ..., --schema=...	which schema in the database
	-c ..., --commit=...	1 or 0(default) specifies commit or not
	-u ..., --unknown=...	the file containing the unknown gene_ids(seperated by ';')
	-p ..., --parser=...	which parser to use
	-g ..., --organism=...	two letter organism abbreviation
	-h, --help              show this help
	
Examples:
	go_bioprocess.py -d mdb -k shu -p shu  -g sc -u yeast_unknown yestprocess2.txt

Description:
	This program extracts go functional class information from raw file,
	combines the stable unknown gene set and sets up the schema.go
	table in the database.
"""

import pickle,sys,os, psycopg, getopt, csv
import Martel
from xml.sax import saxutils
from Martel import LAX

termtype = Martel.Str("TERMTYPE:[") + Martel.ToSep("termtype", "]") + Martel.AnyEol()
id = Martel.Str("id:") + Martel.ToEol("id")
name = Martel.Str("name:") + Martel.ToEol("name")
alt_id = Martel.Str("alt_id:") + Martel.ToEol("alt_id")
namespace = Martel.Str("namespace:") + Martel.ToEol("namespace")
definition = Martel.Str("def:") + Martel.ToEol("def")
comment = Martel.Str("comment:") + Martel.ToEol("comment")
is_a = Martel.Str("is_a:") + Martel.ToEol("is_a")
no_of_genes = Martel.Digits("no_of_genes")
genes = Martel.ToEol("genes")
dataline = Martel.Group("dataline", no_of_genes+\
												Martel.Str(";")+\
												genes)
												
term = Martel.Group("term",  Martel.Rep(definition)+\
										Martel.Rep(is_a)+\
										termtype+\
										id+\
										name+\
										Martel.Rep(alt_id)+\
										namespace+\
										Martel.Rep(definition)+\
										Martel.Rep(comment)+\
										Martel.Rep(is_a)+\
										dataline+\
										Martel.Rep(Martel.AnyEol()) )
go_bioprocess = Martel.Group("go_bioprocess",  Martel.Rep(term), {"format":"go"} )
														
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


class shu_parser:
	'''
	parse the yeast known_gene(GO bioprocess) file,
	construct a known-gene dictionary with gene name as key and GO id as value.
	'''
	def __init__(self):
		self.go_dict = {}
		
	def parse(self, inf, vertex_dict):
		self.term_iterator = go_bioprocess.make_iterator("term")
		self.lax=LAX.LAX()
		for record in self.term_iterator.iterateFile(inf, self.lax):
			gene_array = []
			whole_gene_array = record['genes'][0].split(';')
			for gene in whole_gene_array:
				if gene in vertex_dict:
					gene_array.append(vertex_dict[gene])
			gene_array.sort()
			info = go_term()
			info.name = record['name'][0]
			info.no_of_genes = int(record['no_of_genes'][0])
			info.whole_gene_array = whole_gene_array
			info.gene_array = gene_array
			self.go_dict[record['id'][0]] = info
			
		key_list = self.go_dict.keys()
		key_list.sort()
		for i in range(len(key_list)):
			id = key_list[i]
			self.go_dict[id].no = i+1
		
		return self.go_dict
		#parser = go_bioprocess.make_parser()
		#parser.setContentHandler(saxutils.XMLGenerator())
		#parser.parseFile(self.inf)

parser_map = {"shu":shu_parser(),\
			"ming":ming_parser(),\
			"min":min_parser()}

class go_table_setup:
	def __init__(self, fname, dbname, schema, parser, u_fname, orgn, needcommit=0):
		self.go_inf = open(fname, 'r')
		self.conn = psycopg.connect('dbname=%s'%dbname)
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
		self.vertex_dict = {}
		self.unknown_gene_list = []
		
	def dstruc_loadin(self):
		self.curs.execute("select gene_id, gene_no from graph.gene_id_to_no where organism='%s' "%self.organism)
		rows = self.curs.fetchall()
		for row in rows:
			self.vertex_dict[row[0]] = row[1]
		
		self._unknown_gene_list = self.reader.next()
		for gene in self._unknown_gene_list:
			if gene in self.vertex_dict:
				self.unknown_gene_list.append(self.vertex_dict[gene])
	
	def submit(self):
		string__unknown_gene_list = repr(self._unknown_gene_list)
		string__unknown_gene_list = string__unknown_gene_list.replace("'", '')
		string__unknown_gene_list = '{' + string__unknown_gene_list[1:-1] + '}'
		string_unknown_gene_list = repr(self.unknown_gene_list)
		string_unknown_gene_list = '{' + string_unknown_gene_list[1:-1] + '}'
		
		self.curs.execute("insert into go values('%s', %d, %d, '%s', '%s', '%s')"%\
			('GO:0000004', 0, len(self._unknown_gene_list), 'biological_process unknown', string__unknown_gene_list, string_unknown_gene_list))
		go_dict = self.parser.parse(self.go_inf, self.vertex_dict)
		for term in go_dict:
			string_whole_gene_array = repr(go_dict[term].whole_gene_array)
			string_gene_array = repr(go_dict[term].gene_array)
			string_whole_gene_array = string_whole_gene_array.replace("'", '')
			string_whole_gene_array = '{' + string_whole_gene_array[1:-1] + '}'
			string_gene_array = '{' + string_gene_array[1:-1] + '}'
			go_dict[term].name = go_dict[term].name.replace("'",'')
			self.curs.execute("insert into go values('%s', %d, %d, '%s', '%s', '%s')"%\
				(term, go_dict[term].no, go_dict[term].no_of_genes, go_dict[term].name,\
				string_whole_gene_array, string_gene_array))
		if self.needcommit:
			self.conn.commit()
			
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hd:k:c:u:p:g:", ["help", "dbname=", "schema=", "commit=", "unknown=","parser=","organism="])
	except:
		print __doc__
		sys.exit(2)
	
	dbname = ''
	schema = ''
	commit = 0
	unknown = ''
	parser = ''
	organism = ''
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-k", "--schema"):
			schema = arg
		elif opt in ("-c", "--commit"):
			commit = int(arg)
		elif opt in ("-u", "-unknown"):
			unknown = arg
		elif opt in ("-p", "-parser"):
			parser = arg
		elif opt in ("-g", "-organism"):
			organism = arg
			
	if dbname and schema and unknown and parser and organism and len(args)>0:
		instance = go_table_setup(args[0], dbname, schema, parser, unknown, organism, commit)
		instance.dstruc_loadin()
		instance.submit()
	else:
		print __doc__
		sys.exit(2)
