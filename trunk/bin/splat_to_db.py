#!/usr/bin/env python
"""
Usage: splat_to_db.py -k SCHEMA [OPTION] DATAFILE

Option:
	DATAFILE usually is patterns-splat, output of splat.
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	the table to store the splat results, splat_result(default)
	-c, --commit	commit this database transaction
	-r, --report	report the progress(a number)
	-h, --help              show this help
	
Examples:
	splat_to_db.py -k shu patterns-splat
	splat_to_db.py -k shu -c patterns-splat
	splat_to_db.py -k sc_yh60_splat_5 -c -t splat_result_4 patterns-splat
	
Description:
	Parse the splat results and import into schema.splat_result.

"""


import sys,os,cStringIO,psycopg,getopt
from codense.common import db_connect

class splat_result_iterator:
	'''looping over a splat result file, generate a single pattern result'''
	def __init__(self, inf):
		self.inf = inf
		self.pattern = ''
	def __iter__(self):
		return self		
	def next(self):
		self.read()
		return cStringIO.StringIO(self.pattern)
	def read(self):
		self.pattern = ''
		line = self.inf.readline()
		while line != '\n':
			if line == '':
				raise StopIteration
				break
			self.pattern += line
			line = self.inf.readline()
		self.pattern += line
	
class splat_to_db:
	'''
	03-19-05
		modify toward module-reuse direction
	'''
	def __init__(self, infname=None, hostname='zhoudb', dbname='graphdb', schema=None,\
		table='splat_result', prefix=None, report=0, needcommit=0):
		self.no_of_edges = ''
		self.recurrence_pattern = ''
		self.recurrence_array = []
		self.edge_set = []
		
		self.infname = infname
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.report = int(report)
		self.needcommit = int(needcommit)
		self.table = table
		self.prefix = prefix
					
	def parse(self, pattern):
		self.recurrence_array = []
		line = pattern.readline()
		no_in_string = line[:-1]
		self.no_of_edges = int(no_in_string)
		line = pattern.readline()
		self.recurrence_pattern = line[:-1]
		for i in range(len(self.recurrence_pattern)):
			if self.recurrence_pattern[i] == '1':
				self.recurrence_array.append(i+1)
		
		line = pattern.readline()
		self.edge_set = []	#initialize the edge_set structure
		while line != '\n':
			list = line.split(')')
			list.pop()	#throw away the last '\n'
			for item in list:
				l_list = item[1:].split()
				edge = [int(l_list[0]), int(l_list[1])]
				self.edge_set.append(edge)
			line = pattern.readline()
			
	def run(self):
		"""
		03-19-05
		"""
		self.inf = open(self.infname, 'r')
		(conn, curs) = self.db_connect(self.hostname, self.dbname, self.schema)

		#if table is not splat_result, create it like splat_result
		if self.table != 'splat_result':
			try:
				curs.execute("create table %s(\
					splat_id		serial primary key,\
					no_of_edges	integer,\
					recurrence_pattern	bit varying(200),\
					recurrence_array	integer[],\
					edge_set	integer[][],\
					connectivity	float)"%self.table)
			except:
				sys.stderr.write("Error occurred when creating table %s\n"%self.table)
		iter = splat_result_iterator(self.inf)
		no = 0
		for pattern in iter:
			self.parse(pattern)
			string_edge_set = repr(self.edge_set)
			string_edge_set = string_edge_set.replace('[','{')
			string_edge_set = string_edge_set.replace(']','}')
			string_recurrence_array = repr(self.recurrence_array)
			string_recurrence_array = '{'+string_recurrence_array[1:-1]+'}'
			try:
				curs.execute("insert into %s(no_of_edges, \
							recurrence_pattern, recurrence_array, edge_set) values (%d,B'%s','%s','%s')"%\
							(self.table, self.no_of_edges, self.recurrence_pattern,\
							string_recurrence_array, string_edge_set ))
			except:
				sys.stderr.write('Error occurred when inserting pattern. Aborted.\n')
				conn.rollback()
				sys.exit(1)
			
			no+=1
			if self.report and no%1000==0:
				sys.stderr.write('%s%d'%('\x08'*20, no))
		if self.report:
			sys.stderr.write('%s%d'%('\x08'*20, no))
		if self.needcommit:
			conn.commit()
		sys.stderr.write('\n\tTotal patterns: %d\n'%no)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hrz:d:k:t:p:c", ["help", "report", "hostname=", "dbname=", "schema=", "table=", "prefix=", "commit"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	commit = 0
	table = 'splat_result'
	prefix = ''
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
		elif opt in ("-t", "--table"):
			table = arg
		elif opt in ("-p", "--prefix"):
			prefix = arg
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-r", "--report"):
			report = 1
	if schema and len(args)==1:
		instance = splat_to_db(args[0], hostname, dbname, schema, table, prefix, report, commit)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
