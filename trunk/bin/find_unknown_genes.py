#!/usr/bin/env python
"""
Usage: find_unknown_genes.py -g ORGANISM FILE

Option:
	FILE is the file to store the biological_process unknown genes.
	-d ..., --dbname=...	the database name, graphdb(default)
	-g ..., --organism=...	two letter organism abbreviation
	-h, --help              show this help
	
Examples:
	find_unknown_genes.py -d graphdb -g sc sc_unknown

Description:
	This program finds all the biological_process unknown genes in the gene_pool
	from all the datasets.
	It depends on two tables, gene_id_to_no and raw_association.
"""

import sys, os, psycopg, getopt, csv
from sets import Set

class find_unknown_genes:
	'''
	'''
	def __init__(self, dbname, orgn, unknown_file):
		self.conn = psycopg.connect('dbname=%s'%dbname)
		self.curs = self.conn.cursor()
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
		self.unknown_file = open(unknown_file, 'w')
		#data structure to store biological_process known genes
		self.known_genes_set = Set()
		#data structure to store biological_process unknown genes
		self.unknown_genes_set = Set()
	
	def dstruc_loadin(self):
		#select only those biological_process known genes form table raw_association
		#also the associated go term is not obsolete.
		self.curs.execute("select a.gene_id from graph.raw_association a, go.term t \
			where a.go_id = t.acc and t.term_type='biological_process' and t.is_obsolete=0 and a.go_id!='GO:0000004'")
		rows = self.curs.fetchall()
		for row in rows:
			self.known_genes_set.add(row[0])
		
		#select those biological_process unknown genes from table raw_association.
		self.curs.execute("select gene_id from graph.raw_association where go_id='GO:0000004'")
		rows = self.curs.fetchall()
		for row in rows:
			self.unknown_genes_set.add(row[0])
		
	def run(self):
		self.dstruc_loadin()
		'''
		some genes in table gene_id_to_no are not assigned to either known or unknown.
		Here I regard them as unknown.
		'''
		self.curs.execute("select gene_id from graph.gene_id_to_no where organism='%s'"%
			self.organism)
		rows = self.curs.fetchall()
		for row in rows:
			if row[0] not in self.known_genes_set:
				self.unknown_genes_set.add(row[0])
		sys.stderr.write('Total %d biological_process known genes.\n'%len(self.known_genes_set))
		sys.stderr.write('Total %d biological_process unknown genes.\n'%len(self.unknown_genes_set))
		self.unknown_file.write('%s\n'%';'.join(self.unknown_genes_set))
		
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hd:g:", ["help", "dbname=", "organism="])
	except:
		print __doc__
		sys.exit(2)
	
	dbname = 'graphdb'
	organism = ''
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-g", "--organism"):
			organism = arg
			
	if organism and len(args) == 1:
		instance = find_unknown_genes(dbname, organism, args[0])
		instance.run()
	else:
		print __doc__
		sys.exit(2)
