#!/usr/bin/env python
"""
Usage: mcl_to_db.py -k SCHEMA -s SIZE [OPTION] MCLINPUTDIR

Option:
	MCLINPUTDIR is directory containing the mcl output files.
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --type=...	result type, 'mcl'(default) or 'mcode'
	-c, --commit	commit the database transaction
	-r, --report	report the progress(a number)
	-h, --help              show this help
	
Examples:
	mcl_to_db.py -k shu gph_result/mcl_input
	mcl_to_db.py -k shu -t mcode gph_result/pattern-splat.out

Description:
	Parse the mcl or mcode output results, which are prefixed by 'out' and 
	import results into schema.mcl_result.

"""

import sys, os, psycopg, csv, cStringIO, getopt

class mcl_iterator:
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

class mcode_iterator:
	'''looping over a mcode result file, generate a single pattern result'''
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

class mcl_to_db:
	def __init__(self, dir, hostname, dbname, schema, type, report, needcommit=0):
		self.dir = os.path.abspath(dir)
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		try:
			self.curs.execute("drop index splat_id_idx")
		except psycopg.ProgrammingError, error:
			self.conn.rollback()
			self.curs.execute("set search_path to %s"%schema)
		self.curs.execute("truncate mcl_result")
		self.curs.execute("alter sequence mcl_result_mcl_id_seq restart with 1")
		self.type = type
		self.report = int(report)
		self.needcommit = int(needcommit)
		self.splat_id = ''
		self.cluster_set = []
		self.parameter = ''
		self.connectivity = 0
		self.mcl_to_be_inserted = []
		
		#iterator mapping
		self.iterator_dict = {'mcl': mcl_iterator,
			'mcode': mcode_iterator}
		#parser mapping
		self.parser_dict = {'mcl': self.mcl_parser,
			'mcode': self.mcode_parser}
		
	def submit(self):
		f_list = os.listdir(self.dir)
		no = 0
		for fname in f_list:
			if fname.find('out') == 0:
				path = os.path.join(self.dir, fname)
				inf = open(path, 'r')
				iter = self.iterator_dict[self.type](inf)
				for mcl_cluster_block in iter:
					self.parser_dict[self.type](mcl_cluster_block)
					for i in range(len(self.cluster_set)):

						vertex_set = self.cluster_set[i]
						string_vertex_set = repr(vertex_set)
						string_vertex_set = string_vertex_set.replace('[','{')
						string_vertex_set = string_vertex_set.replace(']','}')
						self.mcl_to_be_inserted.append([self.splat_id, string_vertex_set, self.parameter, self.connectivity])
						no+=1
					if self.report:
						sys.stderr.write('%s%s'%('\x08'*20, no))
		if self.needcommit:
			for entry in self.mcl_to_be_inserted:
				try:
					self.curs.execute("insert into mcl_result(splat_id, vertex_set, parameter, connectivity) \
					values (%s,'%s','%s',%f)"%(entry[0], entry[1], entry[2], entry[3]))
				except:
					sys.stderr.write('Error occured when inserting pattern. Aborted.\n')
					self.conn.rollback()
					sys.exit(1)
			self.curs.execute("create index splat_id_idx on mcl_result(splat_id)")
			self.conn.commit()
		sys.stderr.write('\n\tTotal clusters: %d\n'%no)

	def mcode_parser(self, inf):
		reader = csv.reader(inf, delimiter='\t')
		#first row is the splat_id
		row = reader.next()
		if row == []:
			#the last empty block returned by mcode_iterator
			return
		self.splat_id = row[0]
		#second row is the vertices delimited by '\t'
		row = reader.next()
		#throw away the last empty element caused by the last '\t'
		row.pop()
		self.cluster_set = [map(int, row)]
		#third row is the number of nodes and the number of edges
		row = reader.next()
		no_of_nodes = float(row[0])
		no_of_edges = float(row[1])
		self.connectivity = 2*no_of_edges/(no_of_nodes*(no_of_nodes-1))
		

	def mcl_parser(self,inf):
		line = inf.readline()
		self.cluster_set = []
		cluster_begin = 0
		vertex_list = []
		vertex_line = ''
		while line:
			if line.find('(splat_id') == 0:
				self.splat_id = int(line[10:-3])
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
		opts, args = getopt.getopt(sys.argv[1:], "hrcz:d:k:t:", ["help", "report", "hostname=", "dbname=", "schema=", "type=", "commit"])
	except:
		print __doc__
		sys.exit(2)

	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	type = 'mcl'
	commit = 0
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
		elif opt in ("-t", "--type"):
			type = arg
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-r", "--report"):
			report = 1

	if schema and len(args)==1:
		instance = mcl_to_db(args[0], hostname, dbname, schema, type, report, commit)
		instance.submit()

	else:
		print __doc__
		sys.exit(2)
