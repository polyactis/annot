#!/usr/bin/env python
"""
Usage: cluster_prune.py -k SCHEMA [OPTION]

Option:
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	the database table storing the clusters, mcl_result(default)
	-p, --prune_type	prune_on_recurrence_array, default is prune_on_edge_set
	-c, --commit	commit the database transaction
	-r, --report	report the progress(a number)
	-h, --help              show this help
	
Examples:
	cluster_prune.py -k shu -c
	cluster_prune.py -k ming1 -p 1 -r
	cluster_prune.py -k shu_whole -p 0 -c -r
Description:
	THis program will discard redundant clusters.
	Criteria to discard:
	Among clusters with same vertex_set and same recurrence_pattern,
	retain the one with highest connectivity.
"""

import sys,os,psycopg,getopt
from kjbuckets import *

class cluster_prune:

	def __init__(self, dbname, schema, table, prune_type, report, needcommit=0):
		self.conn = psycopg.connect('dbname=%s'%dbname)
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.table = table
		self.prune_type = int(prune_type)
		self.report = int(report)
		self.needcommit = int(needcommit)
		self.vertex_set_dict = {}
		#store a recurrence_pattern_dict with the same vertex_set
		self.mcl_id_2_connectivity_dict ={}
		#store the mcl_id:connectivity pair
		self.mcl_id_2_recurrence_array_dict = {}
		self.bad_mcl_id_list = []
		#store the mcl_ids that are going to be deleted.
		self.good_mcl_id_list = []
		self.log_file = open("/tmp/cluster_prune.log", 'w')
		self.no_of_records = 0
		
	def dstruc_loadin(self):
		self.curs.execute("DECLARE crs CURSOR FOR select m.mcl_id, m.connectivity, s.recurrence_array from %s m, splat_result s\
			where m.splat_id=s.splat_id"%self.table)
		self.curs.execute("fetch 5000 from crs")
		rows = self.curs.fetchall()
		i = 0
		while rows:
			for row in rows:
				self.mcl_id_2_connectivity_dict[row[0]] = row[1]
				recurrence_array = row[2][1:-1].split(',')
				recurrence_array = map(int, recurrence_array)
				self.mcl_id_2_recurrence_array_dict[row[0]] = kjSet(recurrence_array)
				i += 1
			if self.report:
				sys.stderr.write("%s%s"%("\x08"*20, i))
			self.curs.execute("fetch 5000 from crs")
			rows = self.curs.fetchall()
		if self.report:
			sys.stderr.write('\n')
	
	def prune_on_edge_set(self):
		self.curs.execute("DECLARE crs_e CURSOR FOR select m.mcl_id, m.vertex_set, s.edge_set \
			from %s m, splat_result s where m.splat_id = s.splat_id"%self.table)
		self.curs.execute("fetch 5000 from crs_e")
		rows = self.curs.fetchall()
		i = 0
		while rows:
			for row in rows:
				mcl_id = row[0]
				vertex_set = row[1]
				edge_set = row[2]
				if vertex_set not in self.vertex_set_dict:
					self.vertex_set_dict[vertex_set] = {}
				if edge_set not in self.vertex_set_dict[vertex_set]:
					self.vertex_set_dict[vertex_set][edge_set] = mcl_id
					self.good_mcl_id_list.append(mcl_id)
				else:
					old_mcl_id = self.vertex_set[vertex_set][edge_set]
					self.log_file.write("%d: %d\n"%(old_mcl_id, mcl_id))
					self.mcl_id_2_recurrence_array_dict[old_mcl_id] += self.mcl_id_2_recurrence_array_dict[mcl_id]
				i += 1
			if self.report:
				sys.stderr.write("%sFetching%s"%("\x08"*20,i))
			self.curs.execute("fetch 5000 from crs_e")
			rows = self.curs.fetchall()
		if self.report:
			sys.stderr.write("\n")
		for mcl_id in self.good_mcl_id_list:
			recurrence_array = self.mcl_id_2_recurrence_array_dict[mcl_id].items()
			recurrence_array.sort()
			self.log_file.write("%d: %s\n"%(mcl_id, repr(recurrence_array)))
			self.curs.execute("update %s set recurrence_array = ARRAY%s where mcl_id =%d"%\
				(self.table, repr(recurrence_array), mcl_id))
		
	def prune_on_recurrence_array(self):
		self.curs.execute("DECLARE crs_r CURSOR FOR select mcl_id, vertex_set, recurrence_array \
			from %s where recurrence_array notnull"%self.table)
		self.curs.execute("fetch 5000 from crs_r")
		rows = self.curs.fetchall()
		i = 0
		while rows:
			for row in rows:
				mcl_id = row[0]
				vertex_set = row[1]
				recurrence_array = row[2]
				if vertex_set not in self.vertex_set_dict:
					self.vertex_set_dict[vertex_set] = {}
				if recurrence_array not in self.vertex_set_dict[vertex_set]:
					self.vertex_set_dict[vertex_set][recurrence_array] = mcl_id
				else:
					old_mcl_id = self.vertex_set[vertex_set][recurrence_array]
					if self.mcl_id_dict[mcl_id] > self.mcl_id_dict[old_mcl_id]:
						self.vertex_set_dict[vertex_set][recurrence_array] = mcl_id
						self.bad_mcl_id_list.append(old_mcl_id)
					else:
						self.bad_mcl_id_list.append(mcl_id)
				i += 1
			if self.report:
				sys.stderr.write("%sFetching%s"%("\x08"*20,i))
			self.curs.execute("fetch 5000 from crs_r")
			rows = self.curs.fetchall()
		if self.report:
			sys.stderr.write("\n")
		for mcl_id in self.bad_mcl_id_list:
			self.curs.execute("delete from %s where mcl_id=%d"%(self.table, mcl_id))
		#self.curs.executemany("delete from "+self.table+" where mcl_id=%d", self.bad_mcl_id_list)
		
	def run(self):
		self.dstruc_loadin()
		if prune_type == 0:
			self.prune_on_edge_set()
		elif prune_type == 1:
			self.prune_on_recurrence_array()
		if self.needcommit:
			self.conn.commit()
			
		sys.stderr.write("Total records deleted: %d\n"%len(self.bad_mcl_id_list))
		sys.stderr.write("Total records updated: %d\n"%len(self.good_mcl_id_list))

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hrd:k:pct:", ["help", "report", "dbname=", "schema=", "prune_type", "commit", "table="])
	except:
		print __doc__
		sys.exit(2)
	
	dbname = 'graphdb'
	schema = ''
	prune_type = 0
	commit = 0
	report = 0
	table = 'mcl_result'
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-k", "--schema"):
			schema = arg
		elif opt in ("-p", "--prune_type"):
			prune_type = 1
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-t", "--table"):
			table = arg

	if schema:
		instance = cluster_prune(dbname, schema, table, prune_type, report, commit)
		instance.run()

	else:
		print __doc__
		sys.exit(2)
