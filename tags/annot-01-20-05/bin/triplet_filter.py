#!/usr/bin/env python
"""
Usage: triplet_filter.py -k SCHEMA [OPTION] DATADIR NEWDIR

Option:
	DATADIR is the directory of the graph construction results.
	NEWDIR is the directory to store the filtered results.
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --type=...	0(default, only known), 1(known+unknown) or 2(all)
			if 2 is selected, specify the organism
	-g ..., --organism=...	two letter organism abbreviation
	-h, --help              show this help
	
Examples:
	triplet_filter.py -k shu normal/ known/
	triplet_filter.py -k shu -t 2 -g sc normal/ known/	take all yeast genes

Description:
	This program depends on schema.gene.
	It does the filtration for the triplets found by triplet_counting.py.
	Three types of filtration:
	0: only known genes in that schema are preserved.
	1: only known genes and unknown genes in that schema are preserved.
	2: all the genes in table gene_id_to_no are preserved.

Note:
	There're some known genes that are not in schema.
	So 1 & 2 are different.
"""


import sys, pickle, os, psycopg, getopt

			
class triplet_filter:
	'''
	This class filters triplets containing schema-specific unknown genes out.
	'''
	def __init__(self, hostname, dbname, schema, type, orgn):
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.type = int(type)
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
		if self.type == 2:
			self.organism = self.org_short2long[orgn]
		self.global_vertex_dict = {}

	def dstruc_loadin(self):
		'''
		This method loads in the data structures from database.
		'''
		if self.type == 0:
			self.curs.execute("select gene_id, gene_no from gene where known=TRUE order by gene_no")
		if self.type == 1:
			self.curs.execute("select gene_id, gene_no from gene order by gene_no")
		if self.type == 2:
			self.curs.execute("select gene_id, gene_no from graph.gene_id_to_no \
				where organism='%s' order by gene_no"%self.organism)
			
		rows = self.curs.fetchall()
		for row in rows:
			self.global_vertex_dict[row[1]] = row[0]
		return 0
		
	def transform(self, inf, outf):
		'''
		'''
		line = inf.readline()
		while line:
			list = line[:-1].split(',')
			if int(list[0]) in self.global_vertex_dict and int(list[1]) in self.global_vertex_dict and int(list[2]) in self.global_vertex_dict:
				outf.write(line)
			line = inf.readline()

def transform_batch(dbname, schema, type, orgn, dir, output_dir):
	files = os.listdir(dir)
	sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
	if not os.path.isdir(output_dir):
		os.makedirs(output_dir)
	
	instance = triplet_filter(dbname, schema, type, orgn)
	data_not_loaded = instance.dstruc_loadin()
	if data_not_loaded:
		return
	
	for f in files:
		pathname = os.path.join(dir, f)
		sys.stderr.write("%d/%d:\t%s\n"%(files.index(f)+1,len(files),f))
		inf = open(pathname, 'r')
		outf = open(os.path.join(output_dir,f), 'w')
		instance.transform(inf, outf)
		inf.close()
		outf.close()
		

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:g:", ["help", "hostname=", "dbname=", "schema=", "type=", "organism="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	type = 0
	organism = ''
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
		elif opt in ("-g", "-organism"):
			organism = arg
			
	if schema and len(args) == 2:
		if type == 2 and organism == '':
			print __doc__
			sys.exit(2)
		transform_batch(hostname, dbname, schema, type, organism, args[0], args[1])
	else:
		print __doc__
		sys.exit(2)
