#!/usr/bin/env python
"""
Usage: mcl_to_db.py -k SCHEMA -s SIZE [OPTION] MCLINPUTDIR

Option:
	MCLINPUTDIR is directory containing the mcl output files.
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-f ..., --form=...	result type, 'mcl'(default), 'mcode' or 'kMax'
	-t ..., --table=...	the table to store the mcl results, mcl_result(default)
	-c, --commit	commit the database transaction
	-r, --report	report the progress(a number)
	-h, --help              show this help
	
Examples:
	mcl_to_db.py -k shu gph_result/mcl_input
	mcl_to_db.py -k shu -f mcode gph_result/mcode
	mcl_to_db.py -k sc_yh60_splat -f kMax -t mcl_result_kMax gph_result/mcl_intput_kMax

Description:
	Parse the mcl or mcode output results, which are prefixed by 'out' and 
	import results into schema.mcl_result.

"""

import sys, os, psycopg, csv, cStringIO, getopt
from splat_to_db import splat_result_iterator
from sets import Set

class kMax_iterator:
	'''looping over a kMax_batch_run.py output file'''
	def __init__(self, inf):
		self.inf = inf
		self.pattern = ''
	def __iter__(self):
		return self		
	def next(self):
		self.read()
		return cStringIO.StringIO(self.pattern)
	def read(self):
		kMax_begin = 0
		self.pattern = ''
		line = self.inf.readline()
		while line != '>\n':
			if line == '':
				raise StopIteration
				break
				
			if kMax_begin:
				#this must be ahead of the next if block
				self.pattern += line
			if line[0] == 't':
				kMax_begin = 1
				self.pattern += line
			line = self.inf.readline()

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

class cluster:
	def __init__(self, splat_id=0, vertex_set=None, parameter='', connectivity=0):
		self.splat_id = splat_id
		self.vertex_set = vertex_set
		self.parameter = parameter
		self.connectivity = connectivity

class mcl_to_db:
	def __init__(self, dir, hostname, dbname, schema, type, table, report, needcommit=0):
		self.dir = os.path.abspath(dir)
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.type = type
		self.table = table
		self.report = int(report)
		self.needcommit = int(needcommit)
		try:
			self.curs.execute("drop index %s_splat_id_idx"%self.table)
		except psycopg.ProgrammingError, error:
			sys.stderr.write("Warning: drop index failed.\n")
			self.conn.rollback()
			self.curs.execute("set search_path to %s"%schema)
		self.curs.execute("truncate mcl_result")
		self.curs.execute("alter sequence mcl_result_mcl_id_seq restart with 1")
		
		#iterator mapping
		self.iterator_dict = {'mcl': mcl_iterator,
			'mcode': mcode_iterator,
			'kMax':kMax_iterator}
		#parser mapping
		self.parser_dict = {'mcl': self.mcl_parser,
			'mcode': self.mcode_parser,
			'kMax': self.kMax_parser}
		
	def submit(self):
		#if table is not splat_result, create it like splat_result
		if self.table != 'mcl_result':
			try:
				self.curs.execute("create table %s(\
					mcl_id	serial primary key,\
					splat_id	integer,\
					vertex_set	integer[],\
					parameter	varchar,\
					connectivity	float,\
					p_value_min	float,\
					go_no_vector	integer[],\
					unknown_gene_ratio	float,\
					recurrence_array	integer[])"%self.table)
			except:
				sys.stderr.write("Error occurred when creating table %s\n"%self.table)
		f_list = os.listdir(self.dir)
		no = 0
		for fname in f_list:
			if fname.find('out') == 0:
				sys.stderr.write("%d/%d:\t%s\n"%(f_list.index(fname)+1,len(f_list),fname))
				path = os.path.join(self.dir, fname)
				inf = open(path, 'r')
				iter = self.iterator_dict[self.type](inf)
				for mcl_cluster_block in iter:
					cluster_set = self.parser_dict[self.type](mcl_cluster_block)
					for cluster in cluster_set:

						vertex_set = cluster.vertex_set
						string_vertex_set = '{' + repr(vertex_set)[1:-1] + '}'
						entry = [cluster.splat_id, string_vertex_set, cluster.parameter, cluster.connectivity]
						#try:
						self.curs.execute("insert into %s(splat_id, vertex_set, parameter, connectivity) \
							values (%d,'%s','%s',%f)"%(self.table, entry[0], entry[1], entry[2], entry[3]))
						'''except:
							sys.stderr.write('Error occured when inserting pattern. Aborted.\n')
							self.conn.rollback()
							sys.exit(1)'''
						no+=1
						if self.report and no%1000==0:
							sys.stderr.write('%s%s'%('\x08'*20, no))
		if self.report:
			sys.stderr.write('%s%s'%('\x08'*20, no))
		if self.needcommit:
			self.curs.execute("create index %s_splat_id_idx on %s(splat_id)"%(self.table, self.table))
			self.conn.commit()
		sys.stderr.write('\n\tTotal clusters: %d\n'%no)

	def mcode_parser(self, inf):
		#initlialize the cluster_set data structure
		cluster_set = []
		reader = csv.reader(inf, delimiter='\t')
		#first row is the splat_id
		row = reader.next()
		'''This is bad, a bug in min's mcode wrapper program
		if row == []:
			#the last empty block returned by mcode_iterator
			return
		'''
		splat_id = int(row[0])
		#second row is the vertices delimited by '\t'
		row = reader.next()
		#throw away the last empty element caused by the last '\t'
		row.pop()
		vertex_set = map(int, row)
		#third row is the number of nodes and the number of edges
		row = reader.next()
		no_of_nodes = float(row[0])
		no_of_edges = float(row[1])
		connectivity = 2*no_of_edges/(no_of_nodes*(no_of_nodes-1))
		#append the unit cluster into the cluster_set
		unit = cluster(splat_id=splat_id, vertex_set=vertex_set, connectivity=connectivity)
		cluster_set.append(unit)
		return cluster_set

	def mcl_parser(self, inf):
		#initlialize the cluster_set data structure
		cluster_set = []
		line = inf.readline()
		cluster_begin = 0
		vertex_list = []
		vertex_line = ''
		while line:
			if line.find('(splat_id') == 0:
				splat_id = int(line[10:-3])
			if line.find('(parameter') == 0:
				parameter = line[11:-3]
			if cluster_begin and (line == ')\n'):
				break;
			elif cluster_begin:
				if line[-2] == '$':
					vertex_line += line
					vertex_list = vertex_line.split()
					vertex_list.pop(0)	#throw away the cluster no.
					vertex_list.pop()		#throw away '$'
					if len(vertex_list) >2:
						vertex_list = map(int, vertex_list)
						unit = cluster(splat_id, vertex_list, parameter)
						cluster_set.append(unit)
					vertex_line = ''
				else:
					vertex_line += line
			if line == 'begin\n':
				cluster_begin = 1
			line = inf.readline()
		return cluster_set

	def kMax_parser(self, inf):
		#initlialize the cluster_set data structure
		cluster_set = []
		#for dumb purpose
		parameter = ''
		#first line contains the splat_id
		line = inf.readline()
		splat_id = int(line[2:-1])
		#second line contains the parameter, old kMax_batch_run.py doesn't output the parameter
		line = inf.readline()
		parameter = line[2:-1]
		inf2 = cStringIO.StringIO(inf.read())
		iter = splat_result_iterator(inf2)
		for pattern in iter:
			#initlialize the set first
			vertex_set = Set()
			#first line is no_of_edges
			line = pattern.readline()
			no_of_edges = int(line[:-1])
			#second line is recurrence pattern, ignore.
			line = pattern.readline()
			#from the third line, it's the edge_set
			line = pattern.readline()
			while line !='\n':
				edge_set = line[1:-2].split(')(')
				for edge in edge_set:
					vertex_list = edge.split()
					vertex_set.add(int(vertex_list[0]))
					vertex_set.add(int(vertex_list[1]))
				line = pattern.readline()
			#convert the set to list
			vertex_list = list(vertex_set)
			vertex_list.sort()
			no_of_nodes = float(len(vertex_list))
			connectivity = 2*no_of_edges/(no_of_nodes*(no_of_nodes-1))
			unit = cluster(splat_id, vertex_list, parameter, connectivity)
			cluster_set.append(unit)
		return cluster_set

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hrz:d:k:f:t:c", ["help", "report", "hostname=", "dbname=", "schema=", "form=", "table=", "commit"])
	except:
		print __doc__
		sys.exit(2)

	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	type = 'mcl'
	table = 'mcl_result'
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
		elif opt in ("-t", "--table"):
			table = arg
		elif opt in ("-f", "--form"):
			type = arg
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-r", "--report"):
			report = 1

	if schema and len(args)==1:
		instance = mcl_to_db(args[0], hostname, dbname, schema, type, table, report, commit)
		instance.submit()

	else:
		print __doc__
		sys.exit(2)
