#!/usr/bin/env python
"""
Usage:	connectivity_no_splat.py -k SCHEMA -t TABLE -g GRAPHDIR [OPTION]
	connectivity_no_splat.py -k SCHEMA -t fim_result -a CHOICE [OPTION]

Option:
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	'fim_result' or 'mcl_result', compute which table's connectivity
	-g ..., --graphdir=...	the directory contains all graphs in gspan format.
	-a ..., --avg=...	CHOICE is 'p'(avg the top #support connectivities) or
		'a'(avg all the connectivities)
	-c, --commit	commit the database transaction
	-r, --report	report the progress(a number)
	-h, --help	show this help
	
Examples:
	connectivity_no_splat.py -k shu -c -t mcl_result -g gph_result/sc_yh60_fp_type1
		:compute the connectivity of mcl_result
	connectivity_no_splat.py -k ming1 -r -t fim_result -g gph_result/sc_yh60_fp_type1
		:compute the connectivity_array of fim_result
	connectivity_no_splat.py -k sc_yh60_fp -c -t fim_result -a p
		:compute the connectivity of fim_result, not avg all.

Description:
	This program computes the connectivity of mcl_result and fim_result.
	5000 records are done in one time. Reduce the memory usage.
	This is done through postgresql's CURSOR mechanism.
	Two kinds of usage:
	First is to compute connectivity of mcl_result and connectivity_array of fim_result.
	Second is to compute the connectivity of fim_result.

"""

import sys, os, psycopg, re, getopt
from kjbuckets import *

class connectivity_no_splat:
	'''
	This program computes the connectivity of mcl_result and fim_result.
	5000 records are done in one time. Reduce the memory usage.
	This is done through postgresql's CURSOR mechanism.
	'''
	def __init__(self, dbname, schema, table, graphdir, avg, report, needcommit=0):
		self.conn = psycopg.connect('dbname=%s'%dbname)
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.table = table
		#mapping between table and update functions
		update_dict = {'fim_result': self.fim_atom_update,
			'mcl_result': self.mcl_atom_update}
		if self.table not in update_dict:
			sys.stderr.write("%s doesn't exist.\n"%self.table)
			sys.exit(2)
		self._atom_update = update_dict[self.table]
		self.graphdir = graphdir
		self.avg = avg
		self.report = int(report)
		self.needcommit = int(needcommit)
		#a data structure to store graphs corresponding to datasets
		self.dataset2graph = {}
		#a data structure to store vertex_sets corresponding to datasets
		self.dataset2vertex_set = {}
		#list of dataset no. in ascending order
		self.dataset_no_list = []
		#a kjSet to store vertices in all datasets
		self.vertex_pool = kjSet()
		self.no_of_fim_records = 0
		self.no_of_mcl_records = 0
		#mapping between table and records counter
		records_dict = {'fim_result': self.no_of_fim_records,
			'mcl_result': self.no_of_mcl_records}
		self.records = records_dict[self.table]
		self.p_no = re.compile(r'\d+$')
		self.log_file = open('/tmp/connectivity_no_splat.log', 'w')
		#mapping between table and connectivity indices
		idx_dict = {'fim_result': 'connectivity_fim_result_idx',
			'mcl_result': 'connectivity_mcl_result_idx'}
		self.connectivity_idx = idx_dict[self.table]
		try:
			self.curs.execute("drop index %s"%self.connectivity_idx)
		except psycopg.ProgrammingError, error:
			self.conn.rollback()
			self.curs.execute("set search_path to %s"%schema)
			
	def dstruc_loadin(self):
		sys.stderr.write("Loading Data STructure...\n")
		
		f_list = os.listdir(self.graphdir)
		for fname in f_list:
			try:
				no = int(self.p_no.search(fname).group())	#integer conversion
			except AttributeError, error:
				sys.stderr.write('%s\n'%error)
				sys.stderr.write("can't get dataset no. from the filename %s\n"%fname)
				sys.exit(2)
			if no in self.dataset2graph:
				sys.stderr.write('dataset %d already exists. ???\n'%no)
				sys.exit(2)
			sys.stderr.write("%d/%d:\t%s\n"%(f_list.index(fname)+1,len(f_list),fname))
			self.dataset2graph[no] = kjGraph()
			path_fname = os.path.join(self.graphdir, fname)
			inf = open(path_fname, 'r')
			for line in inf:
				if line[0] == 'e':
					#edge here, like 'e 3807 3859 0.804645'
					line_list = line[:-1].split()
					vertex1 = int(line_list[1])
					vertex2 = int(line_list[2])
					#store the edge in ascending order
					if vertex1 < vertex2:
						self.dataset2graph[no].add(vertex1, vertex2)
					else:
						self.dataset2graph[no].add(vertex2, vertex1)
			#collect all vertices into a kjSet
			self.dataset2vertex_set[no] = kjSet(self.dataset2graph[no].keys() + self.dataset2graph[no].values())
			inf.close()
		#initialize the dataset no list
		self.dataset_no_list = self.dataset2graph.keys()
		self.dataset_no_list.sort()
		sys.stderr.write("Done\n")	
	
	
	def run(self):
		if self.avg:
			self.second_run()
		else:
			self.first_run()
		
	def second_run(self):
		'''
		minor connectivity computing function
		only for fim_result, compute its connectivity from connectivity_array
		'''
		if self.table != 'fim_result':
			sys.stderr.write('--avg only works on fim_result table\n')
			sys.exit(2)
		sys.stderr.write("Compute the connectivity of %s, avg type: %s.\n"%(self.table, self.avg))
		self.curs.execute("begin")
		self.curs.execute("DECLARE crs CURSOR FOR select mcl_id, connectivity_array, support from %s"%self.table)
		self.curs.execute("fetch 5000 from crs")
		rows = self.curs.fetchall()
		while rows:
			for row in rows:
				self._second_run(row)
			if self.report:
				sys.stderr.write("%s%s"%("\x08"*20,self.records))
			self.curs.execute("fetch 5000 from crs")
			rows = self.curs.fetchall()
	
	def _second_run(self, row):
		mcl_id = row[0]
		connectivity_array = row[1][1:-1].split(',')
		connectivity_array = map(float, connectivity_array)
		support = row[2]
		if avg == 'p':
			#average the #support biggest connectivities
			connectivity_array.sort()
			connectivity = sum(connectivity_array[-support:])/support
		elif avg == 'a':
			#average all
			connectivity = sum(connectivity_array)/len(connectivity_array)
		else:
			sys.stderr.write("avg type  %s, not available\n"%avg)
			sys.exit(2)
		#log it
		self.log_file.write("%d\t%f\n"%(mcl_id, connectivity))
		try:
			self.curs.execute("update %s set connectivity=%f where mcl_id=%d"% \
			(self.table, connectivity, mcl_id))
		except:
			sys.stderr.write('Error occurred while setting connectivity\n')
			sys.exit(1)
		self.records += 1

	def first_run(self):
		'''
		major connectivity computing function, mcl_result's connectivity and fim_result's connectivity_array and recurrence_array
		'''
		self.dstruc_loadin()
		self.curs.execute("begin")
		self.curs.execute("DECLARE crs CURSOR FOR select mcl_id, splat_id, vertex_set from %s"%self.table)
		self.curs.execute("fetch 5000 from crs")
		rows = self.curs.fetchall()
		while rows:
			for row in rows:
				self._atom_update(row)

			if self.report:
				sys.stderr.write("%s%s"%("\x08"*20,self.records))
			self.curs.execute("fetch 5000 from crs")
			rows = self.curs.fetchall()

		sys.stderr.write("\tNumber of distinct genes in all patterns: %d\n"%len(self.vertex_pool) )

		if self.needcommit:
			self.curs.execute("create index %s on %s(connectivity)"%(self.connectivity_idx, self.table))
			self.curs.execute("end")
			sys.stderr.write('\n\tTotal %d records updated\n'%self.records)
		else:
			self.conn.rollback()
			sys.stderr.write('\n\tNo real updates\n')
		
	def mcl_atom_update(self, row):
		'''
		Compute the connectivity of mcl_result based on the graph it's from.
		'''
		mcl_id = row[0]
		splat_id = row[1]
		vertex_list = row[2][1:-1].split(',')
		vertex_list = map(int, vertex_list)
		#vertices in ascending order, cause edges are in ascending order
		vertex_list.sort()
		vertex_set = kjSet(vertex_list)	
		#throw it in the vertex_pool
		self.vertex_pool += vertex_set
		no_of_vertices = len(vertex_list)
		no_of_edges = 0
		edge_list = []
		for i in xrange(no_of_vertices):
			for j in xrange(i+1, no_of_vertices):
				if self.dataset2graph[splat_id].member(vertex_list[i],vertex_list[j]):
					edge_list.append([vertex_list[i], vertex_list[j]])
					no_of_edges += 1
		connectivity = 2.0*no_of_edges/((no_of_vertices-1)*no_of_vertices)
		#log it
		self.log_file.write("%d\t%f\n"%(mcl_id, connectivity))	
		try:
			self.curs.execute("update %s set connectivity=%f where mcl_id=%d"% \
			(self.table, connectivity, mcl_id))
			#insert the edge_list of a mcl cluster into splat_result as a pseudo splat result
			#leave the possibility of computing the connectivity of fim_result.
			self.curs.execute("insert into splat_result(splat_id, no_of_edges, edge_set) \
				values(%d, %d, ARRAY%s)"%(mcl_id, no_of_edges, repr(edge_list)))
		except:
			sys.stderr.write('Error occurred while setting connectivity\n')
			sys.exit(1)
		self.records += 1

	def fim_atom_update(self, row):
		'''
		first find which graph contains this vertex_set
		second compute the corresponding connectivity
		finally, get a connectivity_array, recurrence_array is its by_product
		'''
		mcl_id = row[0]
		splat_id = row[1]
		vertex_list = row[2][1:-1].split(',')
		vertex_list = map(int, vertex_list)
		vertex_list.sort()
		vertex_set = kjSet(vertex_list)
		#throw it in the vertex_pool
		self.vertex_pool += vertex_set
		no_of_vertices = len(vertex_list)
		recurrence_array = []
		connectivity_array = []
		for no in self.dataset_no_list:
			if vertex_set.subset(self.dataset2vertex_set[no]):
				no_of_edges = 0
				recurrence_array.append(no)
				for i in xrange(no_of_vertices):
					for j in xrange(i+1, no_of_vertices):
						if self.dataset2graph[no].member(vertex_list[i],vertex_list[j]):
							no_of_edges += 1
				connectivity_array.append( 2.0*no_of_edges/((no_of_vertices-1)*no_of_vertices) )
		#log it
		self.log_file.write("%d\t%s\t%s\n"%(mcl_id, repr(recurrence_array), repr(connectivity_array)))
		try:
			self.curs.execute("update %s set recurrence_array=ARRAY%s, connectivity_array=ARRAY%s where mcl_id=%d"% \
			(self.table, repr(recurrence_array), repr(connectivity_array), mcl_id))
		except:
			sys.stderr.write('Error occurred while setting fim connectivity\n')
			sys.exit(1)
		self.records += 1
	

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hd:k:t:g:a:rc", ["help", "dbname=", "schema=", "table=", "graphdir=", "avg=", "report", "commit"])
	except:
		print __doc__
		sys.exit(2)
	
	dbname = 'graphdb'
	schema = ''
	table = ''
	graphdir = ''
	avg = None
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
		elif opt in ("-t", "--table"):
			table = arg
		elif opt in ("-g", "--graphdir"):
			graphdir = arg
		elif opt in ("-a", "--avg"):
			avg = arg
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-r", "--report"):
			report = 1


	if schema and table and graphdir:
		instance = connectivity_no_splat(dbname, schema, table, graphdir, avg, report, commit)
		instance.run()
	elif schema and table and avg:
		instance = connectivity_no_splat(dbname, schema, table, graphdir, avg, report, commit)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
