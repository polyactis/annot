#!/usr/bin/env python
"""
Usage: fill_dataset_no2id.py -k -i -o [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ...,	table ('dataset_no2id', default)
	-m ...,	mapping file,('~/mapping/$schema\_datasets_mapping', default)
	-c,	commit the database transaction
	-h, --help              show this help

Examples:

Description:
	Fill in table schema.dataset_no2id  based on the mapping_file.
"""

import os, sys, csv, getopt, re
sys.path += [os.path.expanduser('~/script/annot/bin')]
from codense.common import db_connect

class fill_dataset_no2id:
	"""
	10-31-05
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, table='dataset_no2id',\
		mapping_file=None, commit=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.table = table
		self.mapping_file = mapping_file
		
		self.commit = int(commit)
	
	def get_dataset_no2id(self, mapping_file):
		sys.stderr.write("Getting dataset_no2id...\n")
		#to get the number out of the dataset name
		p_no = re.compile(r'\d+$')
		reader = csv.reader(open(mapping_file), delimiter = '\t')
		dataset_no2id = {}
		for row in reader:
			dataset_id, dataset_no = row
			dataset_no = int(p_no.search(dataset_no).group())
			dataset_no2id[dataset_no] = dataset_id
		sys.stderr.write("End getting dataset_no2id.\n")
		return dataset_no2id
	
	def submit_to_table(self, curs, table, dataset_no2id):
		sys.stderr.write("Submitting dataset_no2id to table...\n")
		for dataset_no, dataset_id in dataset_no2id.iteritems():
			curs.execute("insert into %s(dataset_no, dataset_id) values(%s, '%s')"%(table, dataset_no, dataset_id))
		sys.stderr.write("End submitting.\n")
	
	def run(self):
		conn, curs = db_connect(self.hostname, self.dbname, self.schema)
		if self.mapping_file==None:
			self.mapping_file = os.path.expanduser('~/mapping/%s_datasets_mapping'%self.schema)
		dataset_no2id = self.get_dataset_no2id(self.mapping_file)
		self.submit_to_table(curs, self.table, dataset_no2id)
		if self.commit:
			curs.execute("end")

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema="]
	opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:m:c", long_options_list)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = 'dataset_no2id'
	mapping_file = None
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
		elif opt in ("-t"):
			table = arg
		elif opt in ("-m"):
			mapping_file = arg
		elif opt in ("-c"):
			commit = 1
		
	if schema:
		instance = fill_dataset_no2id(hostname, dbname, schema, table, mapping_file, commit)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
