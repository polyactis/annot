#!/usr/bin/env python
"""
Usage: kmax_edge2vertex.py -k SCHEMA -s SOURCE_TABLE -t TARGET_TABLE [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-s ..., --source_table=...	which source table
	-t ..., --target_table=...	which target table
	-c, --commit	commit the database transaction
	-r, --report	report the progress(a number)
	-h, --help              show this help
	
Examples:
	kmax_edge2vertex.py -k shu -c -s kmax_result -t mcl_result_kmax
	
Description:
	transfer the kmax results stored in a splat_result-like table(source_table)
	to a mcl_result-like table(target_table).
	
"""

import sys, os, psycopg, getopt
from sets import Set

class kmax_edge2vertex:
	'''
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
		self.no_of_records = 0
		#data structure to store all vertices
		self.vertex_pool = Set()
		self.log_file = open('/tmp/kmax_edge2vertex.log', 'w')
					
	def run(self):
		#create a mcl_result-like table
		try:
			self.curs.execute("create table %s(\
				mcl_id	integer,\
				splat_id	integer,\
				vertex_set	integer[],\
				parameter	varchar,\
				connectivity	float,\
				p_value_min	float,\
				go_no_vector	integer[],\
				unknown_gene_ratio	float,\
				recurrence_array	integer[])"%self.target_table)
		except:
			sys.stderr.write("Error occurred when creating table %s\n"%self.target_table)
		self.curs.execute("begin")
		self.curs.execute("DECLARE crs CURSOR FOR select splat_id, splat_id, edge_set, connectivity,recurrence_array from %s"%self.source_table)
		self.curs.execute("fetch 5000 from crs")
		rows = self.curs.fetchall()
		while rows:
			for row in rows:
				#list() convert the tuple to list, as tuple doesn't support item assignment.
				self._kmax_edge2vertex(list(row))

			if self.report:
				sys.stderr.write("%s%s"%("\x08"*20,self.no_of_records))
			self.curs.execute("fetch 5000 from crs")
			rows = self.curs.fetchall()

		sys.stderr.write("\tNumber of distinct genes in all splat patterns: %d\n"%len(self.vertex_pool) )

		if self.needcommit:
			if self.target_table:
				self.curs.execute("create index %s_connectivity_idx on %s(connectivity)"%(self.target_table, self.target_table))
			self.curs.execute("end")
			sys.stderr.write('\n\tTotal %d records updated\n'%self.no_of_records)
		else:
			self.conn.rollback()
			sys.stderr.write('\n\tNo real updates\n')

	def _kmax_edge2vertex(self, row):
		local_vertex_pool = self.dstruc_from_edge_set(row[2])
		vertex_list = list(local_vertex_pool)
		vertex_list.sort()
		string_vertex_set = '{' + repr(vertex_list)[1:-1] + '}'
		row[2] = string_vertex_set
		try:
			self.curs.execute("insert into " + self.target_table + "(mcl_id, splat_id, vertex_set, connectivity, recurrence_array)\
				values (%d, %d, %s, %f, %s)", row)
		except:
			sys.stderr.write('Error occurred while inserting.\n')
			sys.exit(1)
		self.no_of_records += 1
	
	def dstruc_from_edge_set(self, edge_set):
		local_vertex_pool = Set()
		edge_list = edge_set[2:-2].split('},{')
		for edge in edge_list:
			vertex_list = edge.split(',')
			vertex_list = map(int, vertex_list)
			self.vertex_pool.add(vertex_list[0])
			self.vertex_pool.add(vertex_list[1])
			local_vertex_pool.add(vertex_list[0])
			local_vertex_pool.add(vertex_list[1])
		return local_vertex_pool

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

	if schema and source_table and target_table:
		instance = kmax_edge2vertex(hostname, dbname, schema, source_table, target_table, report, commit)
		instance.run()

	else:
		print __doc__
		sys.exit(2)
