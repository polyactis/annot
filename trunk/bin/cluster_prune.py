#!/usr/bin/env python
"""
Usage: cluster_prune.py -k SCHEMA [OPTION]

Option:
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	the database table storing the clusters, mcl_result(default)
	-p ..., --prune_type=...	0(prune_on_edge_set, default), 1(prune_on_recurrence_array)
		2(prune_on_vertex_set)
	-s ..., --threshhold=...	the memory threshhold, default is 1e6.
	-c, --commit	commit the database transaction
	-r, --report	report the progress(a number)
	-h, --help              show this help
	
Examples:
	cluster_prune.py -k shu -c -s 3e6 -r
		:prune_on_edge_set, threshhold 3e6, on table mcl_result
	cluster_prune.py -k ming1 -p 1 -r
		:prune_on_recurrence_array, on table mcl_result
	cluster_prune.py -k shu_whole -p 2 -t mcl_result2 -r -c
		:prune_on_vertex_set and on table mcl_result2

Description:
	THis program will discard redundant clusters.
	Three types of pruning:
	1. clusters with same vertex_set and edge_set, choose an arbitrary one to
	retain and update its recurrence_array.
	2. Among clusters with same vertex_set and same recurrence_pattern,
	retain the one with highest connectivity.
	3. clusters with same vertex_set will be merged. THIS should be carried upon
	a new table, like mcl_result2 and after above two prunings.
	
"""

import sys,os,psycopg,getopt
from kjbuckets import *

'''
Comments in this program is sort of confusing. Some are for the codes above.
Some are for the codes below.
'''
class mcl_id_struc:
	#deprecated.
	def __init__(self):
		self.subgraph = None
		self.recurrence_array = None
		self.vertex_set = None
		self.connectivity = None

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
		self.prune_func_dict = {0:self.prune_on_edge_set,
			1:self.prune_on_recurrence_array,
			2:self.prune_on_vertex_set}
		#a function mapping structure.
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


	def _prune_on_vertex_set(self):
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
		self.curs.execute("DECLARE %s CURSOR FOR select m.mcl_id, m.vertex_set, m0.recurrence_array, m.connectivity \
			from %s m, mcl_result m0 where m.mcl_id=m0.mcl_id and m.recurrence_array isnull"%(self.crs, self.table))
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
				connectivity = row[3]
				recurrence_array = kjSet(map(int, recurrence_array))
				if vertex_set in self.vertex_set_dict:
				#identify a cluster by its vertex_set
					old_mcl_id = self.vertex_set_dict[vertex_set][0]
					old_connectivity = self.vertex_set_dict[vertex_set][1]
					if connectivity <= old_connectivity:
						#old one is still good
						self.log_file.write("%d swallows %d\n"%(old_mcl_id, mcl_id))
						self.good_mcl_id_dict[old_mcl_id] += recurrence_array
						#merge recurrence_array
						self.bad_mcl_id_list.append(mcl_id)
					else:
						#new one is good, old one is bad.
						self.log_file.write("%d swallows %d\n"%(mcl_id, old_mcl_id))
						self.vertex_set_dict[vertex_set] = [mcl_id, connectivity]
						#store the old recurrence_array first
						old_recurrence_array = self.good_mcl_id_dict[old_mcl_id]
						#push the new one into good_mcl_id_dict
						self.good_mcl_id_dict[mcl_id] = recurrence_array
						#merge with the old.
						self.good_mcl_id_dict[mcl_id] += old_recurrence_array
						#delete the old from the good_mcl_id_dict
						del self.good_mcl_id_dict[old_mcl_id]
						#now old_mcl_id goes to the bad_mcl_id_list
						self.bad_mcl_id_list.append(old_mcl_id)
				else:
					if len(self.good_mcl_id_dict)>self.threshhold:
						self.remaining += 1
						#leave it to the next run
					else:
						self.good_mcl_id_dict[mcl_id] = recurrence_array
						#this time, the first one is not always the good one,
						#it may be replaced by one with higher connectivity
						self.vertex_set_dict[vertex_set] = [mcl_id, connectivity]
				self.log_file.write('%s: %s\n'%(mcl_id, vertex_set))
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
		
	def prune_on_vertex_set(self):
		'''
		First off, clone a table named as self.table from mcl_result. Leave the recurrence_array empty.
		'''
		#make sure no stupid choice
		if self.table == 'mcl_result':
			sys.stderr.write("Not for %s. Please choose mcl_result2.\n"%self.table)
			sys.exit(2)
		#create the table structure from mcl_result
		try:
			self.curs.execute("create table %s(like mcl_result)"%self.table)
		except psycopg.ProgrammingError, error:
			sys.stderr.write('%s\n'%error)
			sys.exit(2)
		#real cloning starts here.
		sys.stderr.write("Cloning mcl_result into %s..."%self.table)
		self.curs.execute("insert into %s(mcl_id, splat_id, vertex_set, parameter, connectivity, p_value_min,\
			go_no_vector, unknown_gene_ratio) select mcl_id, splat_id, vertex_set, parameter, connectivity,\
			p_value_min, go_no_vector, unknown_gene_ratio from mcl_result"%self.table)
		sys.stderr.write("Done\n")
		
		i = 0
		sys.stderr.write("Pruning based on vertex_set...\n")
		while self.remaining != 0:
			self._prune_on_vertex_set()
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
		self.prune_func_dict[self.prune_type]()
		if self.needcommit:
			self.conn.commit()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hrd:k:p:s:ct:", ["help", "report", "dbname=", "schema=", "prune_type=", "threshhold=", "commit", "table="])
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
			prune_type = int(arg)
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
