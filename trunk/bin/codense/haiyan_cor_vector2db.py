#!/usr/bin/env python
"""
Usage: haiyan_cor_vector2db.py -k SCHEMA -p MAPPING_FILE -i COR_VECTOR_FILE -s SIGNIFICANCE_FILE

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	the splat_result table (edge_set)
	-p ..., --mapping_file=...	the file to get the mapping between haiyan's index and my gene_no
	-i ..., --input_file=...	edge correlation vector file, (one input of codense and copath)
	-s ..., --significance_file=...	the edge significance flag vector file,
	-c, --commit	commit this database transaction
	
Examples:
	./script/annot/bin/haiyan_cor_vector2db.py -c -k sc_54  -p sc_54_gene_id2no 
		-i /tmp/sv_f.txt -s /tmp/sig_vector

Description:
	Submit Haiyan's cor_vector file into a database table.
	
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path += [os.path.expanduser('~/lib64/python')]
else:   #32bit
	sys.path += [os.path.expanduser('~/lib/python')]
import getopt, psycopg, csv
from common import *

	
class haiyan_cor_vector2db:
	"""
	03-02-05
		submit cor_vector + sig_vector
	"""
	def __init__(self, hostname, dbname, schema, table, mapping_file, needcommit, infname, significance_file):
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.table = table
		self.mapping_file = mapping_file
		self.needcommit = int(needcommit)
		self.inf = csv.reader(open(infname, 'r'), delimiter='\t')
		self.sig_f = csv.reader(open(significance_file, 'r'), delimiter='\t')
	
	def run(self):
		#setup the self.haiyan_no2gene_no	
		self.haiyan_no2gene_no = get_haiyan_no2gene_no(self.mapping_file)
		for row in self.inf:
			#get a row from the significance_file
			sig_row = self.sig_f.next()
			edge = row[:2]
			edge = map(int, edge)
			edge = dict_map(self.haiyan_no2gene_no, edge)
			edge.sort()
			edge = '{'+ repr(edge)[1:-1] + '}'
			cor_vector = row[2:]
			cor_vector = map(divided_by_1000, cor_vector)
			cor_vector = '{' + ','.join(cor_vector) + '}'
			sig_vector = sig_row[2:]
			sig_vector = '{' + ','.join(sig_vector) + '}'
			

			#inserting into the splat_table
			self.curs.execute("insert into %s(edge_name, cor_vector, sig_vector) \
				values ('%s', '%s', '%s')"%\
				(self.table, edge, cor_vector, sig_vector))
		if self.needcommit:
			self.curs.execute("create index %s_edge_name_idx on %s(edge_name)"%(self.table, self.table))
			self.conn.commit()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema=", "table=", "mapping_file=", \
		"input_file=", "significance_file=", "commit"]
	
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:p:i:s:c", long_options_list)
	except:
		print __doc__
		sys.exit(2)
		
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	commit = 0
	table = 'edge_cor_vector'
	mapping_file = None
	input_file = None
	significance_file = None
	
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
		elif opt in ("-i", "--input_file"):
			input_file = arg
		elif opt in ("-s", "--significance_file"):
			significance_file = arg
		elif opt in ("-c", "--commit"):
			commit = 1

	if schema and mapping_file and input_file and significance_file:
		instance = haiyan_cor_vector2db(hostname, dbname, schema, table, mapping_file, commit, input_file, significance_file)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
