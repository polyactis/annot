#!/usr/bin/env python
"""
Usage: mcl_to_db.py -k SCHEMA -s SIZE [OPTION] MCLINPUTDIR

Option:
	MCLINPUTDIR is directory containing the mcl output files.
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-c, --commit	commit the database transaction
	-r, --report	report the progress(a number)
	-h, --help              show this help
	
Examples:
	mcl_to_db.py -k shu gph_result/mcl_input
	
Description:
	Parse the mcl output results, which are prefixed by 'out' and 
	import results into schema.mcl_result.

"""

import sys,os,psycopg,cStringIO,getopt

class mcl_result_iterator:
	'''looping over a mcl result file, generate a block of clusters related to one splat pattern'''
	def __init__(self, inf):
		self.inf = inf
		self.mcl_cluster_block = ''
	def __iter__(self):
		return self		
	def next(self):
		self.read()
		return cStringIO.StringIO(self.mcl_cluster_block)
	def read(self):
		self.mcl_cluster_block = ''
		mcl_begin = 0
		line = self.inf.readline()
		while line != '\n':
			'''
			\n
			(mclruninfo
			...
			is discarded.
			'''
			if line == '':
				raise StopIteration
				break
			if mcl_begin == 1:
				self.mcl_cluster_block += line
			if line.find('(splat_id') == 0:
				mcl_begin = 1
				self.mcl_cluster_block += line
			line = self.inf.readline()
		self.mcl_cluster_block += line

class mcl_to_db:
	def __init__(self, dir, dbname, schema, report, needcommit=0):
		self.dir = os.path.abspath(dir)
		self.conn = psycopg.connect('dbname=%s'%dbname)
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.report = int(report)
		self.needcommit = int(needcommit)
		self.splat_id = ''
		self.cluster_set = []
		self.parameter = ''
				
	def submit(self):
		f_list = os.listdir(self.dir)
		no = 0
		for fname in f_list:
			if fname.find('out') == 0:
				path = os.path.join(self.dir, fname)
				inf = open(path, 'r')
				iter = mcl_result_iterator(inf)
				for mcl_cluster_block in iter:
					self.parse(mcl_cluster_block)
					for i in range(len(self.cluster_set)):

						vertex_set = self.cluster_set[i]
						self.mcl_id = '%s_%s_%d'%(self.splat_id,self.parameter,(i+1))
						string_vertex_set = repr(vertex_set)
						string_vertex_set = string_vertex_set.replace('[','{')
						string_vertex_set = string_vertex_set.replace(']','}')
						try:
							self.curs.execute("insert into mcl_result(splat_id, vertex_set, parameter) \
							values ('%s','%s','%s')"%(self.splat_id, string_vertex_set, self.parameter))
						
						except:
							sys.stderr.write('Error occured when inserting pattern. Aborted.\n')
							self.conn.rollback()
							sys.exit(1)
						
						no+=1
					if self.report:
						sys.stderr.write('%s%s'%('\x08'*20, no))
		if self.needcommit:
			self.conn.commit()
		sys.stderr.write('\n\tTotal clusters: %d\n'%no)

	def parse(self,inf):
		line = inf.readline()
		self.cluster_set = []
		cluster_begin = 0
		vertex_list = []
		vertex_line = ''
		while line:
			if line.find('(splat_id') == 0:
				self.splat_id = line[10:-3]
			if line.find('(parameter') == 0:
				self.parameter = line[11:-3]
			if cluster_begin and (line == ')\n'):
				break;
			elif cluster_begin:
				if line[-2] == '$':
					vertex_line += line
					vertex_list = vertex_line.split()
					vertex_list.pop(0)	#throw away the cluster no.
					vertex_list.pop()		#throw away '$'
					if len(vertex_list) >2:
						for i in range(len(vertex_list)):
							vertex_list[i] = int(vertex_list[i])
						self.cluster_set.append(vertex_list)
					vertex_line = ''
				else:
					vertex_line += line
			if line == 'begin\n':
				cluster_begin = 1
			line = inf.readline()
				
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
		instance = mcl_to_db(args[0], dbname, schema, report, commit)
		instance.submit()

	else:
		print __doc__
		sys.exit(2)
