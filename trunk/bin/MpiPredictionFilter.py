#!/usr/bin/env python
"""
Usage: PredictionFilterByClusterSize.py -k SCHEMA -s SPLAT_TABLE -m MCL_TABLE
	-p P_GENE_TABLE -t SPLAT_VIEW -n MCL_VIEW -q P_GENE_VIEW [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-s ...,	SPLAT_TABLE
	-m ...,	MCL_TABLE
	-p ...,	P_GENE_TABLE
	-t ...,	SPLAT_VIEW
	-n ...,	MCL_VIEW
	-q ...,	P_GENE_VIEW
	-e ...,	max cluster size, 40 (default)
	-c,	commit the database transaction
	-b,	debug version.
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	PredictionFilterByClusterSize.py -k mm_fim_97 
	
Description:
	Program to filter clusters by max_size. Database transaction
	is done by creating views.

"""

import sys, os, getopt
sys.path += [os.path.expanduser('~/script/annot/bin')]
from codense.common import db_connect
from gene_stat import gene_stat

class PredictionFilterByClusterSize:
	"""
	09-26-05
		upgrade the function-form to class, 
		p_gene_view is not a view anymore, a real table derived from p_gene_table
			because runtime select cluster_size_cut_off<=max_size blows the memory
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, splat_table=None,\
		mcl_table=None, p_gene_table=None, splat_view=None, mcl_view=None, p_gene_view=None, max_size=40, \
		commit=0, debug=0, report=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.splat_table = splat_table
		self.mcl_table = mcl_table
		self.p_gene_table = p_gene_table
		self.splat_view = splat_view
		self.mcl_view = mcl_view
		self.p_gene_view = p_gene_view
		self.max_size = int(max_size)
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
	
	def view_from_table(self, curs, table, view):
		curs.execute("CREATE OR REPLACE VIEW %s AS SELECT * FROM %s"%(view, table))
		
	def new_p_gene_table(self, curs, p_gene_table, p_gene_view, max_size):
		gene_stat_instance = gene_stat()
		gene_stat_instance.createGeneTable(curs, p_gene_view)
		curs.execute("INSERT INTO %s SELECT * from %s where cluster_size_cut_off<=%s"%(p_gene_view, p_gene_table, max_size))
		
	def run(self):
		(conn, curs) =  db_connect(hostname, dbname, schema)
		self.view_from_table(curs, self.splat_table, self.splat_view)
		self.view_from_table(curs, self.mcl_table, self.mcl_view)
		self.new_p_gene_table(curs, self.p_gene_table, self.p_gene_view, self.max_size)
		if commit:
			curs.execute("End")

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:s:m:p:t:n:q:e:cbr", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	splat_table = None
	mcl_table = None
	p_gene_table = None
	splat_view = None
	mcl_view = None
	p_gene_view = None
	max_size = 40
	commit = 0
	debug = 0
	report = 0
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
		elif opt in ("-s"):
			splat_table = arg
		elif opt in ("-m"):
			mcl_table = arg
		elif opt in ("-p"):
			p_gene_table = arg
		elif opt in ("-t"):
			splat_view = arg
		elif opt in ("-n"):
			mcl_view = arg
		elif opt in ("-q"):
			p_gene_view = arg
		elif opt in ("-e"):
			max_size = int(arg)
		elif opt in ("-c"):
			commit = 1
		elif opt in ("-b"):
			debug = 1
		elif opt in ("-r"):
			report = 1
	if schema and splat_table and p_gene_table and mcl_table and splat_view and mcl_view and p_gene_view:
		instance = PredictionFilterByClusterSize(hostname, dbname, schema, splat_table, mcl_table, p_gene_table,\
			splat_view, mcl_view, p_gene_view, max_size, commit, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
