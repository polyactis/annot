#!/usr/bin/env python
"""
Usage: gene_id_to_no.py -d DATABASENAME -g ORGANISM [OPTION] DATADIR

Option:
	-d ..., --dbname=...	the database name
	-c ..., --commit=...	1 or 0(default) specifies commit or not
	-g ..., --organism=...	two letter organism abbreviation
	-h, --help              show this help
	
Examples:
	gene_id_to_no.py -d mdb -g sc datasets/yeast_data/normal
	gene_id_to_no.py -d mdb -g hs -c 1 datasets/hs

Description:
	This program extracts all gene_id's from datasets,
	numbers them and insert gene_id:number pair into
	database table gene_id_to_no.
"""

import sys, os, csv, psycopg, getopt

class gene_id_to_no:
	def __init__(self, dir, dbname, orgn, needcommit=0):
		self.dir = dir
		self.conn = psycopg.connect('dbname=%s'%dbname)
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
		self.needcommit = int(needcommit)
		
		if self.organism == 'Saccharomyces cerevisiae':
			import pickle
			pickle_fname = os.path.join(os.path.expanduser('~'),'pickle/yeast_global_struc')
			global_struc = pickle.load(open(pickle_fname,'r'))
			self.vertex_dict = global_struc['vertex_dict']
			self.no_of_yeast_genes = len(self.vertex_dict)
		else:
			self.vertex_dict = {}
			
	def run(self):
		files = os.listdir(self.dir)
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
		new_yeast_gene_list = []
		for f in files:
			sys.stderr.write("%d/%d:\t%s\n"%(files.index(f)+1,len(files),f))
			f_path = os.path.join(self.dir, f)
			reader = csv.reader(file(f_path), delimiter='\t')
			for row in reader:
				if row[0] not in self.vertex_dict:
					if self.organism == 'Saccharomyces cerevisiae':
						self.vertex_dict[row[0]] = 0
						new_yeast_gene_list.append(row[0])
					else:
						self.vertex_dict[row[0]] = 1
			del reader
		
		if self.organism == 'Saccharomyces cerevisiae':
			new_yeast_gene_list.sort()
			for gene in new_yeast_gene_list:
				self.vertex_dict[gene] = self.no_of_yeast_genes + new_yeast_gene_list.index(gene) +1
		else:
			key_list = self.vertex_dict.keys()
			key_list.sort()
			for i in xrange(len(key_list)):
				self.vertex_dict[key_list[i]] = i+1
				
		if self.needcommit:
			self.submit()
	
	def submit(self):
		for item in self.vertex_dict:
			self.curs.execute("insert into gene_id_to_no(gene_id, gene_no, organism) values('%s', %d, '%s')"%\
				(item, self.vertex_dict[item], self.organism))
		self.conn.commit()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hd:g:c:", ["help", "dbname=", "organism=", "commit="])
	except:
		print __doc__
		sys.exit(2)
	
	dbname = ''
	commit = 0
	organism = ''
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-g", "--organism"):
			organism = arg
		elif opt in ("-c", "--commit"):
			commit = int(arg)
			
	if dbname and organism and len(args)>0:
		instance = gene_id_to_no(args[0], dbname, organism, commit)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
