#!/usr/bin/env python
"""
Usage: cluster_prune.py -k SCHEMA [OPTION]

Option:
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	the database table storing the clusters, mcl_result(default)
	-p, --prune_type	prune_on_recurrence_array, default is prune_on_edge_set
	-s ..., --threshhold=...	the memory threshhold, default is 1e6.
	-c, --commit	commit the database transaction
	-r, --report	report the progress(a number)
	-h, --help              show this help
	
Examples:
	cluster_prune.py -k shu -c -s 3e6 -r
	cluster_prune.py -k ming1 -p -r
	cluster_prune.py -k shu_whole -p -c -r

Description:
	THis program will discard redundant clusters.
	Two Criteria to discard:
	1. clusters with same vertex_set and edge_set, choose an arbitrary one to
	retain and update its recurrence_array.
	2. Among clusters with same vertex_set and same recurrence_pattern,
	retain the one with highest connectivity.
"""

import sys,os,psycopg,getopt
from kjbuckets import *

class mcl_id_struc:
	#deprecated.
	def __init__(self):
		self.subgraph = None
		self.recurrence_array = None
		self.vertex_set = None
		self.isgood = 0

class cluster_prune:

	def __init__(self, dbname, schema, table, prune_type, threshhold, report, needcommit=0):
		self.conn = psycopg.connect('dbname=%s'%dbname)
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.table = table
		self.prune_type = int(prune_type)
		self.threshhold = float(threshhold)
		self.report = int(report)
		self.needcommit = int(needcommit)
		self.vertex_set_dict = {}
		#store a edge_set_dict or recurrence_pattern_dict with the same vertex_set
		self.mcl_id2connectivity_dict ={}
		#store the mcl_id:connectivity pair
		self.good_mcl_id_dict = {}
		#the mcl_id:recurrence_array pair dictionary
		self.bad_mcl_id_list = []
		#store the mcl_ids that are going to be deleted.
		self.log_file = open("/tmp/cluster_prune.log", 'w')
		self.no_of_goods = 0
		self.no_of_bads = 0
		self.run_no = 0
		self.remaining = 1	#1 is useful to let the first run start

	def subgraph_from_edge_set(self, vertex_set, edge_set):
		#in vertex_set, vertices are in ascending order
		#also for each edge of edge_set, two vertices are in ascending order
		graph = {}
		edge_list = edge_set[2:-2].split('},{')
		for edge in edge_list:
			graph[edge] = 1
			#edge in a graph is like '4,5'
		vertex_list = vertex_set[1:-1].split(',')
		no_of_vertices = len(vertex_list)
		subgraph = ''
		for i in range(no_of_vertices):
			for j in range(i+1, no_of_vertices):
				edge = '%s,%s'%(vertex_list[i], vertex_list[j])
				if edge in graph:
					subgraph += '{'+edge+'}'
					#a string like '{1,3}{3,4}' is used to denote a subgraph
		return subgraph
	
	def _prune_on_edge_set(self):
		del self.vertex_set_dict
		del self.good_mcl_id_dict
		del self.bad_mcl_id_list
		self.vertex_set_dict = {}
		self.good_mcl_id_dict = {}
		self.bad_mcl_id_list = []
		self.remaining = 0
		self.run_no += 1
		#data structure initialization
		sys.stderr.write("Run No.%d\n"%self.run_no)
		self.crs = 'crs%d'%self.run_no
		#in one transaction, multiple cursors to avoid cursor name collision
		self.curs.execute("DECLARE %s CURSOR FOR select m.mcl_id, m.vertex_set, s.recurrence_array, s.edge_set\
			from %s m, splat_result s where m.splat_id=s.splat_id and m.recurrence_array isnull"%(self.crs, self.table))
		#m.recurrence_array notnull means the clusters which have same vertex_set
		#and edge_set as m.mcl_id have all been merged into this mcl_id.
		self.curs.execute("fetch 5000 from %s"%self.crs)
		rows = self.curs.fetchall()
		i = 0
		while rows:
			for row in rows:
				mcl_id = row[0]
				vertex_set = row[1]
				recurrence_array = row[2][1:-1].split(',')
				edge_set = row[3]
				recurrence_array = kjSet(map(int, recurrence_array))
				subgraph = self.subgraph_from_edge_set(vertex_set, edge_set)
				if vertex_set in self.vertex_set_dict:
				#identify a cluster first by its vertex_set
					if subgraph in self.vertex_set_dict[vertex_set]:
						old_mcl_id = self.vertex_set_dict[vertex_set][subgraph]
						self.log_file.write("%d swallows %d\n"%(old_mcl_id, mcl_id))
						self.good_mcl_id_dict[old_mcl_id] += recurrence_array
						#merge recurrence_array
						self.bad_mcl_id_list.append(mcl_id)
					else:
						if len(self.good_mcl_id_dict)>self.threshhold:
							self.remaining += 1
							#leave it to the next run
						else:
							self.good_mcl_id_dict[mcl_id] = recurrence_array
							#the first one is always the good one
							self.vertex_set_dict[vertex_set][subgraph] = mcl_id
				else:
					if len(self.good_mcl_id_dict)>self.threshhold:
						self.remaining += 1
						#leave it to the next run
					else:
						self.good_mcl_id_dict[mcl_id] = recurrence_array
						#the first one is always the good one
						self.vertex_set_dict[vertex_set] = {}
						#different from above codes
						self.vertex_set_dict[vertex_set][subgraph] = mcl_id
				self.log_file.write('%s: %s\n'%(mcl_id, subgraph))
				i += 1
			if self.report:
				sys.stderr.write("%s\t%s"%("\x08"*20, i))
			self.curs.execute("fetch 5000 from %s"%self.crs)
			rows = self.curs.fetchall()
		if self.report:
			sys.stderr.write('\n')
		sys.stderr.write("\trecords to be updated: %d\n"%len(self.good_mcl_id_dict))
		sys.stderr.write("\trecords to be deleted: %d\n"%len(self.bad_mcl_id_list))
		sys.stderr.write("\tDatabase transacting...")
		for mcl_id in self.good_mcl_id_dict:
			#update the recurrence_array of the good mcl_ids
			recurrence_array = self.good_mcl_id_dict[mcl_id].items()
			recurrence_array.sort()
			self.log_file.write("%d: %s\n"%(mcl_id, repr(recurrence_array)))
			self.curs.execute("update %s set recurrence_array = ARRAY%s where mcl_id =%d"%\
				(self.table, repr(recurrence_array), mcl_id))
		#delete the bad mcl_ids
		for mcl_id in self.bad_mcl_id_list:
			self.curs.execute("delete from %s where mcl_id=%d"%(self.table, mcl_id))
		sys.stderr.write("\tDone\n")
		self.no_of_goods += len(self.good_mcl_id_dict)
		self.no_of_bads += len(self.bad_mcl_id_list)
		
	def prune_on_edge_set(self):
		i = 0
		sys.stderr.write("Pruning based on edge_set...\n")
		while self.remaining != 0:
			self._prune_on_edge_set()
		sys.stderr.write("Total records updated: %d\n"%self.no_of_goods)
		sys.stderr.write("Total records deleted: %d\n"%self.no_of_bads)
				
	def prune_on_recurrence_array(self):
		self.curs.execute("select mcl_id, connectivity from %s"%self.table)
		rows = self.curs.fetchall()
		for row in rows:
			self.mcl_id2connectivity_dict[row[0]] = row[1]
		
		self.curs.execute("DECLARE crs_r CURSOR FOR select mcl_id, vertex_set, recurrence_array \
			from %s where recurrence_array notnull"%self.table)
		self.curs.execute("fetch 5000 from crs_r")
		rows = self.curs.fetchall()
		if self.report:
			sys.stderr.write("Pruning based on recurrence_array...\n")
		i = 0
		while rows:
			for row in rows:
				mcl_id = row[0]
				vertex_set = row[1]
				recurrence_array = row[2]
				if vertex_set not in self.vertex_set_dict:
				#identify a cluster first by its vertex_set
					self.vertex_set_dict[vertex_set] = {}
				if recurrence_array not in self.vertex_set_dict[vertex_set]:
				#clusters with same vertex_set, but different recurrence_array
					self.vertex_set_dict[vertex_set][recurrence_array] = mcl_id
				else:
				#clusters with same vertex_set, also same recurrence_array
				#then compare the connectivity. 
					old_mcl_id = self.vertex_set[vertex_set][recurrence_array]
					if self.mcl_id2connectivity_dict[mcl_id] > self.mcl_id2connectivity_dict[old_mcl_id]:
						self.vertex_set_dict[vertex_set][recurrence_array] = mcl_id
						self.bad_mcl_id_list.append(old_mcl_id)
					else:
						self.bad_mcl_id_list.append(mcl_id)
					self.no_of_bads += 1
				i += 1
			if self.report:
				sys.stderr.write("%s%s"%("\x08"*20,i))
			self.curs.execute("fetch 5000 from crs_r")
			rows = self.curs.fetchall()
		if self.report:
			sys.stderr.write("\n")
		sys.stderr.write("Total records to be deleted: %d\n"%self.no_of_bads)
		sys.stderr.write("Total records to be updated: %d\n"%self.no_of_goods)
		sys.stderr.write("Database transacting...")
		for mcl_id in self.bad_mcl_id_list:
			self.curs.execute("delete from %s where mcl_id=%d"%(self.table, mcl_id))
		#self.curs.executemany("delete from "+self.table+" where mcl_id=%d", self.bad_mcl_id_list)
		sys.stderr.write("Done\n")
	
	def run(self):
		if prune_type == 0:
			self.prune_on_edge_set()
		elif prune_type == 1:
			self.prune_on_recurrence_array()
		if self.needcommit:
			self.conn.commit()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hrd:k:ps:ct:", ["help", "report", "dbname=", "schema=", "prune_type", "threshhold=", "commit", "table="])
	except:
		print __doc__
		sys.exit(2)
	
	dbname = 'graphdb'
	schema = ''
	prune_type = 0
	threshhold = 1e6
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
		elif opt in ("-s", "--threshhold"):
			threshhold = float(arg)
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-t", "--table"):
			table = arg

	if schema:
		instance = cluster_prune(dbname, schema, table, prune_type, threshhold, report, commit)
		instance.run()

	else:
		print __doc__
		sys.exit(2)
