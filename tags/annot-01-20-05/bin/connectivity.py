#!/usr/bin/env python
"""
Usage: connectivity.py -k SCHEMA -t TABLE [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-s ..., --source_table=...	which splat table
	-t ..., --target_table=...	which mcl table
	-c, --commit	commit the database transaction
	-r, --report	report the progress(a number)
	-h, --help              show this help
	
Examples:
	connectivity.py -k shu -c -s splat_result
	connectivity.py -k shu -c -s splat_result -t mcl_result
	connectivity.py -k sc_yh60_splat -r -s splat_result_sup_5 -t mcl_result
	
Description:
	This program computes the connectivity of splat and mcl results.
	5000 records are done in one time. Reduce the memory usage.
	This is done through postgresql's CURSOR mechanism.
	If only source_table is given, this program computes the connectivity
	of the splat patterns. If both source_table and target_table are given, the
	program computes the connectivity of mcl clusters as the latter depends on
	splat patterns.
	
"""

import sys,os,cStringIO,psycopg,getopt
from sets import Set

class compute_connectivity:
	'''
	a class for computing the connectivity of splat and mcl results.
	1000 records are done in one time. Reduce the memory usage.
	This is done through postgresql's CURSOR mechanism.
	Computing connectivity of mcl results is based on splat_id.
	Because edge_dict is required in advance.
	'''
	def __init__(self, hostname, dbname, schema, source_table, target_table, report, needcommit=0):
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.source_table = source_table
		self.target_table = target_table
		try:
			self.curs.execute("drop index %s_connectivity_idx"%self.target_table)
		except psycopg.ProgrammingError, error:
			sys.stderr.write("Warning: drop the index of connectivity error.\n")
			self.conn.rollback()
			self.curs.execute("set search_path to %s"%schema)
		self.report = int(report)
		self.needcommit = int(needcommit)
		self.vertex_dict = {}
		#a structure to store vertices of a splat pattern
		self.edge_dict = {}
		#structure to store edges of a splat pattern
		self.no_of_records = 0
		self.vertex_pool = Set()
		#store the vertices of all the splat patterns
		self.log_file = open('/tmp/connectivity.log', 'w')
					
	def run(self):
		self.curs.execute("begin")
		self.curs.execute("DECLARE crs CURSOR FOR select splat_id,edge_set from %s"%self.source_table)
		self.curs.execute("fetch 5000 from crs")
		rows = self.curs.fetchall()
		while rows:
			for row in rows:
				splat_id = row[0]
				edge_set = row[1]
				if self.target_table:
					self.mcl_atom_update(splat_id, edge_set)
				else:
					#target_table is not given, so compute only the connectivity of source_table
					self.splat_atom_update(splat_id, edge_set)

			if self.report:
				sys.stderr.write("%s%s"%("\x08"*20,self.no_of_records))
			self.curs.execute("fetch 5000 from crs")
			rows = self.curs.fetchall()

		sys.stderr.write("\tNumber of distinct genes in all splat patterns: %d\n"%len(self.vertex_pool) )
		self.log_file.write(repr(self.vertex_pool))
		if self.needcommit:
			if self.target_table:
				self.curs.execute("create index %s_connectivity_idx on %s(connectivity)"%(self.target_table, self.target_table))
			self.curs.execute("end")
			sys.stderr.write('\n\tTotal %d records updated\n'%self.no_of_records)
		else:
			self.conn.rollback()
			sys.stderr.write('\n\tNo real updates\n')
		
	def mcl_atom_update(self, splat_id, edge_set):
		'''
		Compute the connectivity of all mcl clusters pertaining to the same splat_id.
		'''
		self.dstruc_from_edge_set(edge_set)
		self.curs.execute("select mcl_id, vertex_set from %s where splat_id=%d"%\
			(self.target_table, splat_id))
		rows = self.curs.fetchall()
		for row in rows:
			mcl_id = row[0]
			vertex_set = row[1][1:-1]
			vertex_list = vertex_set.split(',')
			no_of_vertices = len(vertex_list)
			no_of_edges = 0
			for i in xrange(no_of_vertices):
				for j in xrange(i+1, no_of_vertices):
					tuple = (vertex_list[i],vertex_list[j])
					if self.edge_dict.has_key(tuple) or self.edge_dict.has_key(tuple):
						no_of_edges += 1
			connectivity = 2.0*no_of_edges/((no_of_vertices-1)*no_of_vertices)
			try:
				self.curs.execute("update %s set connectivity=%f where mcl_id=%d"% \
				(self.target_table, connectivity, mcl_id))
			except:
				sys.stderr.write('Error occurred while setting mcl connectivity\n')
				sys.exit(1)
			self.no_of_records += 1

	def splat_atom_update(self, splat_id, edge_set):
		self.dstruc_from_edge_set(edge_set)
		no_of_edges = len(self.edge_dict)
		no_of_vertices = len(self.vertex_dict)
		connectivity = 2.0*no_of_edges/((no_of_vertices-1)*no_of_vertices)
		try:
			self.curs.execute("update %s set connectivity=%f where splat_id=%d"% \
			(self.source_table, connectivity, splat_id))
		except:
			sys.stderr.write('Error occurred while setting splat connectivity\n')
			sys.exit(1)
		self.no_of_records += 1
	
	def dstruc_from_edge_set(self, edge_set):
		self.edge_dict = {}
		self.vertex_dict = {}
		edge_list = edge_set[2:-2].split('},{')
		for edge in edge_list:
			vertex_list = edge.split(',')
			vertex_list = (vertex_list[0], vertex_list[1])
			self.edge_dict[vertex_list] = 1
			vertex1 = vertex_list[0]
			vertex2 = vertex_list[1]
			self.vertex_pool.add(vertex1)
			self.vertex_pool.add(vertex2)
			if vertex1 not in self.vertex_dict:
				self.vertex_dict[vertex1] = 1
			if vertex2 not in self.vertex_dict:
				self.vertex_dict[vertex2] = 1

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hrz:d:k:s:t:c", \
			["help", "report", "hostname=", "dbname=", "schema=", "source_table=", "target_table=", "commit"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	source_table = ''
	target_table = ''
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
		elif opt in ("-s", "--source_table"):
			source_table = arg
		elif opt in ("-t", "--target_table"):
			target_table = arg
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-r", "--report"):
			report = 1

	if schema and source_table:
		instance = compute_connectivity(hostname, dbname, schema, source_table, target_table, report, commit)
		instance.run()

	else:
		print __doc__
		sys.exit(2)
