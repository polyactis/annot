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

def PredictionFilterByClusterSize(hostname='zhoudb', dbname='graphdb', schema=None, splat_table=None,\
	mcl_table=None, p_gene_table=None, splat_view=None, mcl_view=None, p_gene_view=None, max_size=40, \
	commit=0, debug=0, report=0):
	(conn, curs) =  db_connect(hostname, dbname, schema)
	curs.execute("CREATE OR REPLACE VIEW %s AS SELECT * FROM %s"%(splat_view, splat_table))
	curs.execute("CREATE OR REPLACE VIEW %s AS SELECT * FROM %s"%(mcl_view, mcl_table))
	curs.execute("CREATE OR REPLACE VIEW %s AS SELECT * FROM %s \
		where cluster_size_cut_off<=%s"%(p_gene_view, p_gene_table, max_size))
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
		PredictionFilterByClusterSize(hostname, dbname, schema, splat_table, mcl_table, p_gene_table,\
			splat_view, mcl_view, p_gene_view, max_size, commit, debug, report)
	else:
		print __doc__
		sys.exit(2)
