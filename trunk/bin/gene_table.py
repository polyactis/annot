#!/usr/bin/env python
"""
Usage: gene_table.py -d DATABASENAME -k SCHEMA -g ORGANISM [OPTION]

Option:
	-d ..., --dbname=...	the database name
	-k ..., --schema=...	which schema in the database
	-c ..., --commit=...	1 or 0(default) specifies commit or not
	-g ..., --organism=...	two letter organism abbreviation
	-h, --help              show this help
	
Examples:
	gene_table.py -d mdb -k shu -g sc

Description:
	Depends on two tables,
	graph.gene_id_to_no and schema.go.
	Sets up table schema.gene.
"""

import pickle, sys, os, psycopg, getopt


class gene_table_setup:
	'''
	Initialize the gene_table of yeast from two pickled data structures,
		known_genes_dict
		global_struc['vertex_list']
	'''
	def __init__(self, dbname, schema, orgn, needcommit=0):
		self.conn = psycopg.connect('dbname=%s'%dbname)
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
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

	def dstruc_loadin(self):
		self.geneno_go_dict = {}
		self.curs.execute("select go_no, gene_array from go")
		rows = self.curs.fetchall()
		for row in rows:
			go_no = row[0]
			gene_no_list = row[1][1:-1].split(',')
			gene_no_list = map(int, gene_no_list)
			for gene_no in gene_no_list:
				if gene_no in self.geneno_go_dict:
					self.geneno_go_dict[gene_no].append(go_no)
				else:
					self.geneno_go_dict[gene_no] = [go_no]
		
		self.geneno_to_geneid = {}
		self.curs.execute("select gene_id, gene_no from graph.gene_id_to_no where organism='%s'"%self.organism)
		rows = self.curs.fetchall()
		for row in rows:
			self.geneno_to_geneid[row[1]] = row[0]
	
	def run(self):
		for gene_no in self.geneno_go_dict:
			gene_id = self.geneno_to_geneid[gene_no]
			known = '1'
			if len(self.geneno_go_dict[gene_no])==1 and self.geneno_go_dict[gene_no][0]==0 :
				known = '0'
			self.curs.execute("insert into gene(gene_id, gene_no, known, \
				go_functions) values ('%s', %d, '%s', ARRAY%s)"%\
				(gene_id, gene_no, known, repr(self.geneno_go_dict[gene_no]) ))
		if self.needcommit:
			self.conn.commit()


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hd:k:c:g:", ["help", "dbname=", "schema=", "commit=", "organism="])
	except:
		print __doc__
		sys.exit(2)
	
	dbname = ''
	schema = ''
	commit = 0
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
		elif opt in ("-g", "-organism"):
			organism = arg
			
	if dbname and schema and organism:
		instance = gene_table_setup(dbname, schema, organism, commit)
		instance.dstruc_loadin()
		instance.run()
	else:
		print __doc__
		sys.exit(2)
