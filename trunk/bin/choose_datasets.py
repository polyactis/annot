#!/usr/bin/env python
"""
Usage: choose_datasets.py  [OPTION] DATADIR

Option:
	DATADIR is the directory of the datasets to be choosen
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database, graph(default)
	-g ..., --organism=...	two letter organism abbreviation, hs(default)
	-h, --help              show this help
	
Examples:
	choose_datasets.py -k graph datasets/hs
	choose_datasets.py -k graph  -g sc normal/

Description:

"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import sys, pickle, os, psycopg, getopt, csv
from sets import Set

class choose_datasets:
	'''
	This class's purpose is to look into gene abundance of individual
	dataset based on the standard gene set got from database.

	'''
	def __init__(self, hostname, dbname, schema, orgn, dir):
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
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
		self.dir = dir
		self.gene_set = Set()
		
		
	def dstruc_loadin(self):
		'''
		This method loads in the data structures from database.
		'''
		self.curs.execute("select gene_id from gene_id_to_no where organism='%s'"%self.organism)
		
		rows = self.curs.fetchall()
		for row in rows:
			self.gene_set.add(row[0])
		return 0
	
	def run(self):
		self.dstruc_loadin()
		
		files = os.listdir(self.dir)
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
		for f in files:
			local_gene_set = Set()
			sys.stderr.write("%d/%d:\t%s"%(files.index(f)+1,len(files),f))
			f_path = os.path.join(self.dir, f)
			reader = csv.reader(file(f_path), delimiter='\t')
			for row in reader:
				local_gene_set.add(row[0])
			#This dataset has some genes that the standard set doesn't have.
			genes_more = local_gene_set - self.gene_set
			len_more = len(genes_more)
			per_more = float(len_more)/len(self.gene_set)
			#The standard set has some genes that this dataset doesn't have.
			genes_less = self.gene_set - local_gene_set
			len_less = len(genes_less)
			per_less = float(len_less)/len(self.gene_set)
			#output in a readable format
			sys.stderr.write("\t%f\t%d\t%f\t%d\n"%(per_more, len_more, per_less, len_less))
			del reader
			
	
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:g:", ["help", "hostname=", "dbname=", "schema=", "organism="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = 'graph'
	organism = 'hs'
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
		elif opt in ("-g", "-organism"):
			organism = arg
			
	if len(args) == 1:
		instance = choose_datasets(hostname, dbname, schema, organism, args[0])
		instance.run()
	else:
		print __doc__
		sys.exit(2)
