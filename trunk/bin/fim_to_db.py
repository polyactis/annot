#!/usr/bin/env python
"""
Usage: fim_to_db.py -k SCHEMA [OPTION] FIM_OUTPUT_FILE

Option:
	FIM_OUTPUT_FILE is a file containing results of fim.
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-c, --commit	commit the database transaction
	-r, --report	report the progress(a number)
	-h, --help              show this help
	
Examples:
	fim_to_db.py -k sc_yh60_fp out.dataset
	
Description:
	Parse the fim ouput files.
	import results into schema.fim_result.

"""

import sys, os, psycopg, csv, getopt


class fim_to_db:
	def __init__(self, infname, dbname, schema, report, needcommit=0):
		self.reader = csv.reader(file(infname), delimiter=' ')
		self.conn = psycopg.connect('dbname=%s'%dbname)
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.curs.execute("truncate fim_result")
		self.curs.execute("alter sequence fim_result_mcl_id_seq restart with 1")
		self.report = int(report)
		self.needcommit = int(needcommit)
		
	def submit(self):
		no = 0
		for row in self.reader:
			#the last item one the line is support, like (10).
			support = int(row.pop()[1:-1])
			if len(row) >= 3:
				row = map(int, row)
				string_row = '{' + repr(row)[1:-1] + '}'
				self.curs.execute("insert into fim_result(vertex_set, support)\
					values('%s', %d)"%(string_row, support))
				no+=1
				if self.report:
					sys.stderr.write('%s%s'%('\x08'*20, no))
		if self.needcommit:
			self.conn.commit()
		sys.stderr.write('\n\tTotal clusters: %d\n'%no)


				
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hrcd:k:", ["help", "report", "dbname=", "schema=", "commit"])
	except:
		print __doc__
		sys.exit(2)
	
	dbname = 'graphdb'
	schema = ''
	commit = 0
	report = 0
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
		elif opt in ("-r", "--report"):
			report = 1

	if schema and len(args)==1:
		instance = fim_to_db(args[0], dbname, schema, report, commit)
		instance.submit()
	else:
		print __doc__
		sys.exit(2)
