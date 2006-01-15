#!/usr/bin/env python
"""
Usage: go_index.py -k SCHEMA [OPTION] go_index_file

Option:
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-c, --commit	commit the database transaction
	-h, --help      show this help
	
Examples:
	go_index.py -k shu process.pathindex.annotation.txt

Description:
	this program fills in go_index in the table schema.go.
	the go_index_file is outputted by Ming-chi's program.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import pickle,sys,os, psycopg, getopt, csv

class go_table_setup:
	def __init__(self, fname, dbname, schema, needcommit=0):
		self.reader = csv.reader(open(fname, 'r'), delimiter='\t')
		self.conn = psycopg.connect('dbname=%s'%dbname)
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.needcommit = int(needcommit)
		self.go_id_dict = {}

	def dstruc_loadin(self):
		self.curs.execute("select go_id from go")
		rows = self.curs.fetchall()
		for row in rows:
			self.go_id_dict[row[0]] = []
	
	def run(self):
		for row in self.reader:
			go_id = 'GO:'+row[0]
			if go_id in self.go_id_dict:
				self.go_id_dict[go_id].append(row[1])
		for item in self.go_id_dict:
			self.curs.execute("update go set go_index=ARRAY%s where go_id='%s'"%(repr(self.go_id_dict[item]), item))
		if self.needcommit:
			self.conn.commit()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hd:k:c", ["help", "dbname=", "schema=", "commit"])
	except:
		print __doc__
		sys.exit(2)
	
	dbname = 'graphdb'
	schema = ''
	commit = 0
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-k", "--schema"):
			schema = arg
		elif opt in ("-c", "--commit"):
			commit = 1
			
	if schema and len(args) == 1:
		instance = go_table_setup(args[0], dbname, schema, commit)
		instance.dstruc_loadin()
		instance.run()
	else:
		print __doc__
		sys.exit(2)
