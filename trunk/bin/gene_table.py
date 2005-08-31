#!/usr/bin/env python
"""
Usage: gene_table.py -k SCHEMA -g ORGANISM [OPTION] DATADIR

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	DATADIR is the directory containing all the datasets
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-g ..., --organism=...	two letter organism abbreviation
	-u, --union	takes the union of all genes in the datasets, default is intersection
	-c, --commit	commits the database transaction
	-h, --help              show this help
	
Examples:
	gene_table.py -k hs_yh60 -g hs datasets/hs_wanted

Description:
	This program sets up schema.gene from the datasets, which
	are probably a subset from that organism's total datasets.
	It depends on table graph.gene_id_to_no.
"""

import pickle, sys, os, psycopg, csv, getopt
from sets import Set

class gene_table:
	'''
	Initialize the local gene_id:gene_no mapping in table schema.gene
	'''
	def __init__(self, dir, hostname, dbname, schema, orgn, union=0, needcommit=0):
		"""
		08-30-05
			add rn to org_short2long
		"""
		self.dir = dir
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.union = int(union)
		self.needcommit = int(needcommit)
		self.org_short2long = {'at':'Arabidopsis thaliana',
			'ce':'Caenorhabditis elegans',
			'dm':'Drosophila melanogaster',
			'hs':'Homo sapiens',
			'mm':'Mus musculus',
			'sc':'Saccharomyces cerevisiae',
			'rn':'Rattus norvegicus',
			'Rattus norvegicus':'Rattus norvegicus',
			'Arabidopsis thaliana':'Arabidopsis thaliana',
			'Caenorhabditis elegans':'Caenorhabditis elegans',
			'Drosophila melanogaster':'Drosophila melanogaster',
			'Homo sapiens':'Homo sapiens',
			'Mus musculus':'Mus musculus',
			'Gorilla gorilla Pan paniscus Homo sapiens':'Homo sapiens',
			'Saccharomyces cerevisiae':'Saccharomyces cerevisiae'}
		self.organism = self.org_short2long[orgn]
		#mapping between gene_id and gene_no
		self.gene_id2gene_no = {}
		#mapping between gene_id and its occurence
		self.gene_id2freq = {}
		#unique gene collection, for database submission
		self.gene_set = Set()
		
	def dstruc_loadin(self):
		#setup self.gene_id2gene_no
		self.curs.execute("select gene_id, gene_no from graph.gene_id_to_no where organism='%s'"%self.organism)
		rows = self.curs.fetchall()
		for row in rows:
			self.gene_id2gene_no[row[0]] = row[1]
	
	def run(self):
		#load in the data structures first.
		self.dstruc_loadin()
		#iterate over all the datasets, find all the genes
		files = os.listdir(self.dir)
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
		new_yeast_gene_list = []
		for f in files:
			sys.stderr.write("%d/%d:\t%s\n"%(files.index(f)+1,len(files),f))
			f_path = os.path.join(self.dir, f)
			reader = csv.reader(file(f_path), delimiter='\t')
			for row in reader:
				if row[0] in self.gene_id2gene_no:
				#03-29-05, ignore genes not in graph.gene_id_to_no.
					if row[0] in self.gene_id2freq:
					#not first encountering
						self.gene_id2freq[row[0]] += 1
					else:
					#first encountering
						self.gene_id2freq[row[0]] = 1
			del reader
		
		if self.union:
		#take the union set
			self.gene_set = Set(self.gene_id2freq.keys())
		else:
		#take the intersection set
			for (gene_id, freq) in self.gene_id2freq.iteritems():
				if freq == len(files):
				#occur in all datasets
					self.gene_set.add(gene_id)
		sys.stderr.write("%d genes to be submitted\n"%len(self.gene_set))
		#database submission
		self.submit()
		
	def submit(self):
		sys.stderr.write("Database transacting...")
		for gene_id in self.gene_set:
			self.curs.execute("insert into gene(gene_id, gene_no) values ('%s', %d)"%\
				(gene_id, self.gene_id2gene_no[gene_id] ))
		if self.needcommit:
			self.conn.commit()
		sys.stderr.write("done.\n")

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:g:uc", ["help", "hostname=", "dbname=", "schema=", "organism=", "union", "commit"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	organism = ''
	union = 0
	commit = 0
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
		elif opt in ("-g", "--organism"):
			organism = arg
		elif opt in ("-u", "--union"):
			union = 1
		elif opt in ("-c", "--commit"):
			commit = 1
			
	if schema and organism and len(args) == 1:
		instance = gene_table(args[0], hostname, dbname, schema, organism, union, commit)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
