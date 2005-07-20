#!/usr/bin/env python
"""
Usage: gene_id_to_no.py -g ORGANISM [OPTION] DATADIR

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-c, --commit	commits the database transaction
	-g ..., --organism=...	two letter organism abbreviation
	-l, --log	record down some stuff in the logfile(gene_id_to_no.log)
	-n, --gene	The gene_id is NCBI Gene ID. (gene_id = gene_no)
	-h, --help              show this help
	
Examples:
	gene_id_to_no.py -g sc datasets/yeast_data/normal
	gene_id_to_no.py -d mdb -g hs -c datasets/hs
	gene_id_to_no.py -g hs -c -n datasets/hs_gene

Description:
	This program extracts all gene_id's from datasets,
	numbers them and insert gene_id:number pair into
	database table gene_id_to_no.
"""

import sys, os, csv, psycopg, getopt
from sets import Set

class gene_id_to_no:
	def __init__(self, dir, hostname, dbname, orgn, log, gene=0, needcommit=0):
		self.dir = dir
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to graph")
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
		self.log = int(log)
		self.gene = int(gene)
		if self.log:
			self.logfile = open('/tmp/gene_id_to_no.log', 'w')
		self.needcommit = int(needcommit)
		#records down the maximum gene_no ever assigned.
		self.max_gene_no = 0
		#mapping between gene_id and gene_no
		self.vertex_dict = {}
		self.vertex_dict_extension = {}
		#a unique collection of all genes in the datasets of the directory
		self.gene_set = Set()
		
	def dstruc_loadin(self):
		'''
		look up the database to find all the assigned genes of this organism.
		this is to maintain backward compatibility.
		Expand the gene pool on the basis of currently assigned genes.
		'''
		self.curs.execute("select gene_id, gene_no from gene_id_to_no where organism='%s'"%self.organism)
		rows = self.curs.fetchall()
		for row in rows:
			if row[1] > self.max_gene_no:
			#this gene_no is bigger than max_gene_no
				self.max_gene_no = row[1]
			self.vertex_dict[row[0]] = row[1]

	def submit(self):
		sys.stderr.write("%d\tgenes to be inserted\n"%(len(self.vertex_dict_extension)) )
		sys.stderr.write("Database transacting...")
		for (gene_id,gene_no) in self.vertex_dict_extension.iteritems():
			self.curs.execute("insert into gene_id_to_no(gene_id, gene_no, organism) values('%s', %d, '%s')"%\
				(gene_id, gene_no, self.organism))
		if self.needcommit:
			self.conn.commit()
		sys.stderr.write("done.\n")
		
	def run(self):
		"""
		07-19-05
			deal with self.gene, gene_id=gene_no
		"""
		#load in the data structure first.
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
				self.gene_set.add(row[0])
			del reader
		if self.gene:
			sys.stderr.write("NCBI Gene ID. gene_id=gene_no.\n")
		else:
			#expand current gene pool
			sys.stderr.write('Old max_gene_no: %d\n'%self.max_gene_no)
		for gene_id in self.gene_set:
			if gene_id not in self.vertex_dict:
				if self.gene:
					self.vertex_dict_extension[gene_id] = int(gene_id)
				else:
					self.max_gene_no += 1
					self.vertex_dict_extension[gene_id] = self.max_gene_no
				if self.log:
					self.logfile.write('%s = %d\n'%(gene_id, self.vertex_dict_extension[gene_id]))
		if self.gene==0:
			sys.stderr.write('New max_gene_no: %d\n'%self.max_gene_no)
		self.submit()


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:g:lnc", \
			["help", "hostname=", "dbname=", "organism=", "log", "gene", "commit"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	commit = 0
	organism = ''
	log = 0
	gene = 0
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-z", "--hostname"):
			hostname = arg
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-g", "--organism"):
			organism = arg
		elif opt in ("-l", "--log"):
			log = 1
		elif opt in ("-n", "--gene"):
			gene = 1
		elif opt in ("-c", "--commit"):
			commit = 1
			
	if dbname and organism and len(args)>0:
		instance = gene_id_to_no(args[0], hostname, dbname, organism, log, gene, commit)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
