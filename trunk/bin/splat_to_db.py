#!/usr/bin/env python
"""
Usage: splat_to_db.py -k SCHEMA -p PREFIX [OPTION] DATAFILE

Option:
	DATAFILE usually is patterns-splat, output of splat.
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-p ..., --prefix=...	specify the splat_id's prefix
	-c ..., --commit=...	0(default) or 1 specifies commit or not
	-h, --help              show this help
	
Examples:
	splat_to_db.py -k shu -p sc_shu patterns-splat
	splat_to_db.py -k shu -p sc_shu -c 1 patterns-splat
	
Description:
	Parse the splat results and import into schema.splat_result.

"""


import sys,os,cStringIO,psycopg,getopt

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
	'''
	def __init__(self, infname, dbname, schema, prefix, needcommit=0):
		self.splat_id = ''
		self.no_of_edges = ''
		self.recurrence_pattern = ''
		self.recurrence_array = []
		self.edge_set = []
		
		self.inf = open(infname, 'r')
		self.needcommit = int(needcommit)
		self.conn = psycopg.connect('dbname=%s'%dbname)
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
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
		iter = splat_result_iterator(self.inf)
		no = 0
		for pattern in iter:
			self.parse(pattern)
			self.splat_id = '%s_%d'%(self.prefix, (no+1))
			string_edge_set = repr(self.edge_set)
			string_edge_set = string_edge_set.replace('[','{')
			string_edge_set = string_edge_set.replace(']','}')
			string_recurrence_array = repr(self.recurrence_array)
			string_recurrence_array = '{'+string_recurrence_array[1:-1]+'}'
			try:
				self.curs.execute("insert into splat_result(splat_id, no_of_edges, \
							recurrence_pattern, recurrence_array, edge_set) values ('%s',%d,B'%s','%s','%s')"%\
							(self.splat_id, self.no_of_edges, self.recurrence_pattern,\
							string_recurrence_array, string_edge_set ))
			except:
				sys.stderr.write('Error occured when inserting pattern. Aborted.\n')
				self.conn.rollback()
				sys.exit(1)
			no+=1
			sys.stderr.write('%s%d'%('\x08'*20, no))
		if self.needcommit:
			self.conn.commit()
		sys.stderr.write('\n\tTotal patterns: %d\n'%no)
		sys.stderr.write('\tLast pattern: %s\n'%self.splat_id)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hd:k:p:c:", ["help", "dbname=", "schema=", "prefix=", "commit="])
	except:
		print __doc__
		sys.exit(2)
	
	dbname = 'graphdb'
	schema = ''
	commit = 0
	prefix = ''
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-k", "--schema"):
			schema = arg
		elif opt in ("-c", "--commit"):
			commit = int(arg)
		elif opt in ("-p", "--prefix"):
			prefix = arg

	if schema and prefix and len(args)==1:
		instance = splat_to_db(args[0], dbname, schema, prefix, commit)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
