#!/usr/bin/env python
"""
Usage: PreFimInput.py -s SIG_VECTOR_FILE [OPTION] OUTPUTFILE

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default, IGNORE)
	-d ..., --dbname=...	the database name, graphdb(default, IGNORE)
	-k ..., --schema=...	which schema in the database, IGNORE
	-s ...,	sig_vector file
	-m ..., --min_sup=...	minimum support of an edge in the database, 0(default)
	-x ..., --max_sup=...	maximum support of an edge, 200(default)
	-t ...,	delimiter for the outputfile, 'tab'(default)
	-y ...,	type, 1(fim_closed, Jianfei Zhu's default), 2(closet+)
	-b, --debug	debug version
	-r,	enable report flag
	-h, --help              show this help
	
Examples:
	PreFimInput.py -s mm_fim_114_4.gene_no.sig_vector -m 5 -x 114 mm_fim_114m5x114_i

Description:
	Prepare input for fim_closed.
"""

import sys, os, psycopg, getopt, csv
from codense.common import db_connect
sys.path += [os.path.expanduser('~/script/annot/bin')]
if sys.version_info[:2] < (2, 3):       #python2.2 or lower needs some extra
        from python2_3 import *

def output_edge_data_fim_edge_oriented(outputfile, min_support=3, \
	max_support=200, hostname='zhoudb', dbname='graphdb', \
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

def output_edge_data_fim_edge_oriented_from_sig_vector_file(outputfile, min_support=3, \
	max_support=200, sig_vector_fname='', delimiter_char='\t', type=1, debug=0, report=0):
	"""
	12-30-05
		read edge data from sig_vector_file
	01-01-06
		add delimiter_char and type
		
		For type==2 (closet+) 's spec file
		1. file name of the data set;
		2. total number of distinct items in the dataset;
		3. Maximum length of one transaction
			(i.e., maximum number of distinct items in each transaction);
		4. total number of transactions
	"""	
	reader = csv.reader(open(sig_vector_fname, 'r'), delimiter='\t')
	writer = csv.writer(open(outputfile, 'w'),delimiter=delimiter_char)
	
	if type==2:	#closet+ needs a spec file
		spec_fname = '%s.spec'%outputfile
		spec_f = open(spec_fname, 'w')
	sys.stderr.write("Getting edge matrix for all functions...\n")
	counter = 0
	real_counter = 0
	for row in reader:
		sig_vector = row[2:]
		sig_vector = map(int, sig_vector)
		if sum(sig_vector)>=min_support and sum(sig_vector)<=max_support:
			new_row = []
			for i in range(len(sig_vector)):
				if sig_vector[i]==1:
					new_row.append(i+1)
			if type==2:	#closet+ needs no of items for each transaction
				new_row.insert(0, len(new_row))
			writer.writerow(new_row)
			real_counter += 1
		counter +=1
		if report and counter%5000==0:
			sys.stderr.write('%s%s'%('\x08'*20, counter))
		if debug:
			break
	if type==2:
		spec_f.write('%s\n'%os.path.basename(outputfile))	#1. file name of the data set;
		spec_f.write('%s\n'%(len(sig_vector)+1))	#2. total number of distinct items in the dataset;
			#01-02-06 +1 because for closet+, item no starts from 0.
		spec_f.write('%s\n'%(len(sig_vector)+1))	#3. Maximum length of one transaction
		spec_f.write('%s\n'%real_counter)	#4. total number of transactions
		del spec_f
	del reader, writer
	sys.stderr.write("Done\n")
	
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:s:m:x:t:y:br", ["help", "hostname=", \
			"dbname=", "schema=", "min_sup", "max_sup=", "debug"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	sig_vector_fname = ''
	min_sup = 0
	max_sup = 200
	delimiter_char = '\t'
	type = 1
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
		elif opt in ("-s",):
			sig_vector_fname = arg
		elif opt in ("-m", "--min_sup"):
			min_sup = int(arg)
		elif opt in ("-x", "--max_sup"):
			max_sup = int(arg)
		elif opt in ("-t",):
			delimiter_char = arg
		elif opt in ("-y",):
			type = int(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r",):
			report = 1
	outputfile = args[0]
	#12-30-05 discard the schema form, take input from sig_vector_fname
	if sig_vector_fname and outputfile:
		output_edge_data_fim_edge_oriented_from_sig_vector_file(outputfile, min_sup, max_sup, \
			sig_vector_fname, delimiter_char, type, debug, report)
	else:
		print __doc__
		sys.exit(2)
