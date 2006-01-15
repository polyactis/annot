#!/usr/bin/env python
"""
Usage: gene_table.py -k SCHEMA -g ORGANISM [OPTION] DATADIR

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	DATADIR is the directory containing all the datasets
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-g ..., --organism=...	two letter organism abbreviation
	-o ...,	output_table(gene, default)
	-m ...,	minimum frequency(0, default, take all genes)
	-c, --commit	commits the database transaction
	-h, --help              show this help
	
Examples:
	gene_table.py -k hs_yh60 -g hs datasets/hs_wanted

Description:
	This program sets up schema.gene from the datasets, which
	are probably a subset from that organism's total datasets.
	It depends on table graph.gene_id_to_no.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import pickle, sys, os, psycopg, csv, getopt
from sets import Set
from codense.common import db_connect, get_gene_id2mt_no_list, get_global_gene_id2gene_no, \
	org2tax_id, org_short2long, get_gene_no2tf_set
	
class gene_table:
	'''
	Initialize the local gene_id:gene_no mapping in table schema.gene
	'''
	def __init__(self, dir=None, hostname='zhoudb', dbname='graphdb', schema=None, \
		output_table='gene', orgn='', min_frequency=0, needcommit=0):
		"""
		08-30-05
			add rn to org_short2long
		"""
		self.dir = dir
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.output_table = output_table
		self.organism = org_short2long(orgn)
		self.min_frequency = int(min_frequency)
		self.needcommit = int(needcommit)		
	
	def return_gene_id_set(self, dir, gene_id2gene_no, min_frequency):
		"""
		09-19-05
			rewrite and split from run()
		09-27-05
			use min_frequency to select genes
		"""
		#iterate over all the datasets, find all the genes
		files = os.listdir(dir)
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
		new_yeast_gene_list = []
		gene_id2freq = {}
		gene_id_set = Set()
		for f in files:
			sys.stderr.write("%d/%d:\t%s\n"%(files.index(f)+1,len(files),f))
			f_path = os.path.join(dir, f)
			reader = csv.reader(file(f_path), delimiter='\t')
			for row in reader:
				if row[0] in gene_id2gene_no:
				#03-29-05, ignore genes not in graph.gene_id_to_no.
					if row[0] in gene_id2freq:
					#not first encountering
						gene_id2freq[row[0]] += 1
					else:
					#first encountering
						gene_id2freq[row[0]] = 1
			del reader
		
		for (gene_id, freq) in gene_id2freq.iteritems():
			if freq >= min_frequency:
				gene_id_set.add(gene_id)
		sys.stderr.write("%d genes to be submitted\n"%len(gene_id_set))
		return gene_id_set
	
	def submit(self, curs, output_table, gene_id_set, gene_id2gene_no, gene_id2mt_no_list):
		sys.stderr.write("Database transacting...")
		for gene_id in gene_id_set:
			mt_no_list = gene_id2mt_no_list.get(gene_id)
			gene_no = gene_id2gene_no[gene_id]
			if mt_no_list:
				mt_no_list.sort()
				curs.execute("insert into %s(gene_id, gene_no, trans_fac)  values\
					('%s', %s, '{%s}')"%(output_table, gene_id, gene_no, repr(mt_no_list)[1:-1]))
			else:
				curs.execute("insert into %s(gene_id, gene_no) values ('%s', %d)"%\
					(output_table, gene_id, gene_no ))
		sys.stderr.write("done.\n")
	
	def run(self):
		"""
		09-19-05
			rewrite
			
			--db_connect()
			--get_global_gene_id2gene_no()
			--org2tax_id()
			--get_gene_id2mt_no_list()
			--return_gene_id_set()
			--submit()
		"""
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		gene_id2gene_no = get_global_gene_id2gene_no(curs, self.organism)
		tax_id = org2tax_id(self.organism)
		"""
		#01-14-06 comment it out for future 
		gene_no2tf_set = get_gene_no2tf_set(curs)	#12-15-05 just yeast.
		#12-15-05 convert gene_no(integer) into gene_id(string)
		gene_id2mt_no_list = {}
		for gene_no, tf_set in gene_no2tf_set.iteritems():
			gene_id2mt_no_list[repr(gene_no)] = list(tf_set)
		"""
		gene_id2mt_no_list = get_gene_id2mt_no_list(tax_id)
		gene_id_set = self.return_gene_id_set(self.dir, gene_id2gene_no, self.min_frequency)
		self.submit(curs, output_table, gene_id_set, gene_id2gene_no, gene_id2mt_no_list)
		if self.needcommit:
			conn.commit()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:o:g:m:c", ["help", "hostname=", "dbname=", "schema=", "organism=", "commit"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	output_table = 'gene'
	organism = ''
	min_frequency = 0
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
		elif opt in ("-o",):
			output_table = arg
		elif opt in ("-g", "--organism"):
			organism = arg
		elif opt in ("-m",):
			min_frequency = int(arg)
		elif opt in ("-c", "--commit"):
			commit = 1
			
	if schema and organism and len(args) == 1 and output_table:
		instance = gene_table(args[0], hostname, dbname, schema, output_table, organism, min_frequency, commit)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
