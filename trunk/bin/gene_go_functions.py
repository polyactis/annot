#!/usr/bin/env python
"""
Usage: gene_go_functions.py -k SCHEMA [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-c, --commit	commits the database transaction
	-g ..., --organism=...	IGNORE it, interface relics
	-h, --help              show this help
	
Examples:
	gene_go_functions.py -k shu
	gene_go_functions.py -k hs_yh60 -c

Description:
	Depends on two tables, schema.gene and schema.go.
	Fills up column go_functions in schema.gene.
"""

import pickle, sys, os, psycopg, getopt

class gene_go_functions:
	'''
	Initialize the go_functions column in table schema.gene
	'''
	def __init__(self, hostname, dbname, schema, orgn, needcommit=0):
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
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
		#mapping between gene_no and GO list
		self.geneno_go_dict = {}
		#mapping between gene_no and gene_id
		self.geneno_to_geneid = {}

	def dstruc_loadin(self):
		#setup self.geneno_go_dict
		self.curs.execute("select go_no, gene_array from go")
		rows = self.curs.fetchall()
		for row in rows:
			go_no = row[0]
			if row[1] == '{}':
				continue
			gene_no_list = row[1][1:-1].split(',')
			gene_no_list = map(int, gene_no_list)
			for gene_no in gene_no_list:
				if gene_no in self.geneno_go_dict:
					self.geneno_go_dict[gene_no].append(go_no)
				else:
					self.geneno_go_dict[gene_no] = [go_no]
		#setup self.geneno_to_geneid
		self.curs.execute("select gene_id, gene_no from gene")
		rows = self.curs.fetchall()
		for row in rows:
			self.geneno_to_geneid[row[1]] = row[0]
	
	def run(self):
		"""
		03-06-05
			change the way dealing with left-out known genes. the old way is deleting them.
			Now. keep them, mark them known, but assign go_no=0. So they are fake unknown genes.
		"""
		#load in the data structures first.
		self.dstruc_loadin()
		excluded_known_genes = 0
		sys.stderr.write("Database transacting...")		
		for (gene_no,gene_id) in self.geneno_to_geneid.iteritems():
			if gene_no not in self.geneno_go_dict:
				'''
				this gene is shared among all datasets. It is known but not present in
				informative nodes. Regard them as fake unknown genes.
				'''
				excluded_known_genes += 1
				self.curs.execute("update gene set known='1', go_functions='{0}' where gene_no=%d"%(gene_no))
				continue
			known = '1'
			if self.geneno_go_dict[gene_no]==[0] :
				known = '0'
			string_geneno_go_dict = repr(self.geneno_go_dict[gene_no])
			string_geneno_go_dict = '{' + string_geneno_go_dict[1:-1] + '}'
			self.curs.execute("update gene set known='%s', go_functions='%s' where gene_no=%d"%\
				(known, string_geneno_go_dict, gene_no))
		if self.needcommit:
			self.conn.commit()
		sys.stderr.write("done.\n")
		sys.stderr.write("%d genes are known genes but not in the informative nodes.\n\
			Fake them as unknown genes(known=TRUE, go_functions={0})\n"%excluded_known_genes)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:cg:", ["help", "hostname=", "dbname=", "schema=", "commit", "organism="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	commit = 0
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
		elif opt in ("-g", "--organism"):
			organism = arg
			
	if schema:
		instance = gene_go_functions(hostname, dbname, schema, organism, commit)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
