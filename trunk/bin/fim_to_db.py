#!/usr/bin/env python
"""
Usage: fim_to_db.py -k SCHEMA [OPTION] FIM_OUTPUT_FILE

Option:
	FIM_OUTPUT_FILE is a file containing results of fim.
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-n, --set_connectivity	this will set all fim_result's connectivity to 1.
	-c, --commit	commit the database transaction
	-r, --report	report the progress(a number)
	-h, --help              show this help
	
Examples:
	fim_to_db.py -k sc_yh60_fp out.dataset
	
Description:
	Parse the fim ouput files.
	import results into schema.fim_result.
	set_connectivity is to make downstream stats convenient.
	
"""

import sys, os, psycopg, csv, getopt


class fim_to_db:
	def __init__(self, infname, hostname, dbname, schema, set_connectivity, report, needcommit=0):
		self.reader = csv.reader(file(infname), delimiter=' ')
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.curs.execute("truncate fim_result")
		self.curs.execute("alter sequence fim_result_mcl_id_seq restart with 1")
		self.report = int(report)
		self.set_connectivity = int(set_connectivity)
		self.needcommit = int(needcommit)
		
	def submit(self):
		no = 0
		for row in self.reader:
			#the last item one the line is support, like (10).
			support = int(row.pop()[1:-1])
			if len(row) >= 3:
				row = map(int, row)
				string_row = '{' + repr(row)[1:-1] + '}'
				if self.set_connectivity:
					self.curs.execute("insert into fim_result(vertex_set, connectivity, support)\
						values('%s', 1, %d)"%(string_row, support))
				else:
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
		opts, args = getopt.getopt(sys.argv[1:], "hrz:d:k:nc", ["help", "report", "hostname=", "dbname=", "schema=", "set_connectivity", "commit"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	set_connectivity = 0
	report = 0
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
		elif opt in ("-n", "--set_connectivity"):
			set_connectivity = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-c", "--commit"):
			commit = 1


	if schema and len(args)==1:
		instance = fim_to_db(args[0], hostname, dbname, schema, set_connectivity, report, commit)
		instance.submit()
	else:
		print __doc__
		sys.exit(2)
