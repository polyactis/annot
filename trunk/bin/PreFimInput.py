#!/usr/bin/env python
"""
Usage: PreFimInput.py -k SCHEMA [OPTION] OUTPUTFILE

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-m ..., --min_sup=...	minimum support of an edge in the database, 0(default)
	-x ..., --max_sup=...	maximum support of an edge, 200(default)
	-b, --debug	debug version
	-h, --help              show this help
	
Examples:
	PreFimInput.py -k shu ~/bin/hhu_clustering/shu.fim.input

Description:
	Prepare input for fim_closed.
"""

import sys, os, psycopg, getopt, csv
from codense.common import db_connect
sys.path += [os.path.expanduser('~/script/annot/bin')]
if sys.version_info[:2] < (2, 3):       #python2.2 or lower needs some extra
        from python2_3 import *

def output_edge_data_fim_edge_oriented(outputfile, min_support=3, max_support=200, hostname='zhoudb', dbname='graphdb', \
	schema='sc_new_38', debug=0, edge_table='edge_cor_vector'):
	"""
	08-07-05
		copied from misc.py
	"""
	conn, curs = db_connect(hostname, dbname, schema)
	
	fname = outputfile
	writer = csv.writer(open(fname,'w'),delimiter='\t')
	
	sys.stderr.write("Getting edge matrix for all functions...\n")
	curs.execute("DECLARE crs CURSOR FOR select edge_id,sig_vector \
		from %s"%(edge_table))
	curs.execute("fetch 5000 from crs")
	rows = curs.fetchall()
	counter = 0
	while rows:
		for row in rows:
			edge_id = row[0]
			sig_vector = row[1][1:-1].split(',')
			sig_vector = map(int, sig_vector)
			
			if sum(sig_vector)>=min_support and sum(sig_vector)<=max_support:
				new_row = []
				for i in range(len(sig_vector)):
					if sig_vector[i]==1:
						new_row.append(i+1)
				writer.writerow(new_row)
			counter +=1
		sys.stderr.write('%s%s'%('\x08'*20, counter))
		if debug:
			break
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
	
	sys.stderr.write("Done\n")

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:m:x:b", ["help", "hostname=", "dbname=", "schema=", "min_sup", "max_sup=", \
			"debug"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	min_sup = 0
	max_sup = 200
	debug = 0
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
		elif opt in ("-m", "--min_sup"):
			min_sup = int(arg)
		elif opt in ("-x", "--max_sup"):
			max_sup = int(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
	outputfile = args[0]
	
	if schema and outputfile:
		output_edge_data_fim_edge_oriented(outputfile, min_sup, max_sup, hostname, dbname, schema, debug)
	else:
		print __doc__
		sys.exit(2)
