#!/usr/bin/env python
"""
Usage: haiyan_cor_vector2db.py -k SCHEMA -p MAPPING_FILE INPUT_FILE OUTPUT_FILE

Option:
	INPUT_FILE is the file containing all the clusters(Haiyan's format).
	OUTPUT_FILE is the file to contain all the clusters in Kangyu's format.
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	the splat_result table (edge_set)
	-p ..., --mapping_file=...	the file to get the mapping between haiyan's index and my gene_no
	-c, --commit	commit this database transaction
	
Examples:
	./script/annot/bin/haiyan_cor_vector2db.py -c -k sc_54  -p sc_54_gene_id2no /tmp/sv_f.txt /tmp/cor_vector

Description:
	Submit Haiyan's cor_vector file into a database table.
	
"""


import sys, os, getopt, psycopg, csv
from common import *

	
class haiyan_cor_vector2db:
	def __init__(self, hostname, dbname, schema, table, mapping_file, needcommit, infname, outfname):
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.table = table
		self.mapping_file = mapping_file
		self.needcommit = int(needcommit)
		self.inf = csv.reader(open(infname, 'r'), delimiter='\t')
		self.outf = csv.writer(open(outfname, 'w'), delimiter = '\t')
	
	def run(self):
		#setup the self.haiyan_no2gene_no	
		self.haiyan_no2gene_no = get_haiyan_no2gene_no(self.mapping_file)
		for row in self.inf:
			edge = row[:2]
			edge = map(int, edge)
			edge = dict_map(self.haiyan_no2gene_no, edge)
			edge.sort()
			edge = '{'+ repr(edge)[1:-1] + '}'
			cor_vector = row[2:-1]
			cor_vector = map(divided_by_1000, cor_vector)
			cor_vector = '{' + ','.join(cor_vector) + '}'
			new_row = [edge, cor_vector]
			#self.outf.writerow(new_row)
			try:
				#inserting into the splat_table
				self.curs.execute("insert into %s(edge_name, cor_vector) \
							values ('%s', '%s')"%\
							(self.table, edge, cor_vector))
			except:
				sys.stderr.write('Error occurred when inserting pattern. Aborted.\n')
				self.conn.rollback()
				sys.exit(1)
		if self.needcommit:
			self.conn.commit()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema=", "table=", "mapping_file=", "commit"]
	
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:p:c", long_options_list)
	except:
		print __doc__
		sys.exit(2)
		
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	commit = 0
	table = 'edge_cor_vector'
	mapping_file = None
	
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
		elif opt in ("-t", "--table"):
			table = arg
		elif opt in ("-p", "--mapping_file"):
			mapping_file = arg
		elif opt in ("-c", "--commit"):
			commit = 1

	if schema and mapping_file and len(args) == 2:
		instance = haiyan_cor_vector2db(hostname, dbname, schema, table, mapping_file, commit, args[0], args[1])
		instance.run()
	else:
		print __doc__
		sys.exit(2)
