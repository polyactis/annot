#!/usr/bin/env python
"""
Usage: prepare_gene_id2no.py -k SCHEMA [OPTION] OUTPUTFILE

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-h, --help              show this help
	
Examples:
	prepare_gene_id2no.py -k shu ~/bin/hhu_clustering/shu_gene_id2no

Description:
	Output a gene_id2no file in gene_no order. It's used to map
	gene_no to netmine(copath)'s index.
"""

import sys, os, psycopg, getopt, csv
from codense.common import db_connect

class prepare_gene_id2no:
	'''
	05-20-05
	'''
	def __init__(self, hostname, dbname, schema, outputfile):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.outputfile = outputfile
	
	def run(self):
		"""
		"""
		(conn, curs) = db_connect(self.hostname, self.dbname, self.schema)
		writer = csv.writer(open(self.outputfile, 'w'), delimiter='\t')
		
		sys.stderr.write("Preparing gene_id2no for %s..."%self.schema)
		curs.execute("select gene_id, gene_no from gene order by gene_no")
		rows = curs.fetchall()
		for row in rows:
			writer.writerow(row)
		del writer
		sys.stderr.write("done.\n")

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:", ["help", "hostname=", "dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
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
	outputfile = args[0]
	
	if schema and outputfile:
		instance = prepare_gene_id2no(hostname, dbname, schema, outputfile)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
