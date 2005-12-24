#!/usr/bin/env python
"""
Usage: filter_clusters.py -k SCHEMA -n GENE_P_TABLE -p P_GENE_TABLE
	-m MCL_TABLE -l GOOD_CLUSTER_TABLE [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-g ...,	GENE_P_TABLE
	-p ...,	P_GENE_TABLE
	-m ...,	MCL_TABLE
	-l ...,	GOOD_CLUSTER_TABLE
	-o ...,	occurrence_cutoff, 0.8(default)
	-s ...,	qsize, 400000(default)
	-n,	CLUSTER_BS_TABLE is new
	-c,	commit the database transaction
	-b,	debug version.
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	filter_clusters.py -k mm_fim_97 -g gene_p_mm_fim_97m4x200_e5_a60
	-p p_gene_mm_fim_97m4x200_e5 -m mcl_mm_fim_97m4x200 -l good_clusters
	
Description:
	Program to filter clusters using gene_p_table. In fact, only retain clusters
	which are used to do function prediction.

"""

import sys, os, getopt
sys.path += [os.path.expanduser('~/script/annot/bin')]
from codense.common import db_connect, p_gene_id_set_from_gene_p_table
from sets import Set
from Queue import Queue
from threading import *
if sys.version_info[:2] < (2, 3):       #python2.2 or lower needs some extra
	from python2_3 import *
	
class submit_clusters(Thread):
	def __init__(self, cluster_queue, hostname='zhoudb', dbname='graphdb', schema=None, gene_p_table=None,\
		p_gene_table=None, mcl_table=None, good_cluster_table=None,  occurrence_cutoff=0.8, \
		new_table=0, commit=0, debug=0, report=0):
		Thread.__init__(self)
		self.cluster_queue = cluster_queue
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.gene_p_table = gene_p_table
		self.p_gene_table = p_gene_table
		self.mcl_table = mcl_table
		self.good_cluster_table = good_cluster_table
		self.occurrence_cutoff = float(occurrence_cutoff)
		self.new_table = int(new_table)
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
		
	def get_mcl_id2unknown_ratio(self, curs, p_gene_table, p_gene_id_set):
		"""
			10-29-05 use p_gene_id_set to screen p_gene_table
		"""
		sys.stderr.write("Getting mcl_id2unknown_ratio...\n")
		mcl_id2unknown_ratio = {}
		curs.execute("DECLARE crs CURSOR FOR select p.p_gene_id, p.mcl_id, p.unknown_cut_off, \
			p.go_no, p.p_value_cut_off from %s p"%(p_gene_table))
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		counter = 0
		real_counter = 0
		while rows:
			for row in rows:
				p_gene_id, mcl_id, unknown_ratio, go_no, p_value = row
				if p_gene_id in p_gene_id_set:
					if mcl_id not in mcl_id2unknown_ratio:
						mcl_id2unknown_ratio[mcl_id] = [unknown_ratio, Set()]
					prediction_tuple = (go_no, p_value)
					mcl_id2unknown_ratio[mcl_id][1].add(prediction_tuple)
					real_counter += 1
				counter += 1
			if self.report:
				sys.stderr.write("%s%s\t%s"%('\x08'*15, counter, real_counter))
			if self.debug and counter ==5000:
				break
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		curs.execute("Close crs")	#close the CURSOR
		sys.stderr.write("%s clusters. Done.\n"%(len(mcl_id2unknown_ratio)))
		return mcl_id2unknown_ratio
	
	def create_good_cluster_table(self, curs, good_cluster_table):
		sys.stderr.write("Creating table %s...\n"%good_cluster_table)
		curs.execute("create table %s(\
			id	serial primary key,\
			mcl_id	integer,\
			splat_id	integer,\
			vertex_set	integer[],\
			recurrence_array	integer[],\
			recurrence	integer,\
			connectivity	float,\
			unknown_ratio	float,\
			size	integer,\
			go_no_list	integer[],\
			p_value_list	float[],\
			comment	varchar)"%good_cluster_table)
		sys.stderr.write("Done.\n")
	
	def submit_good_clusters(self, curs, cluster_queue, good_cluster_table, mcl_id2unknown_ratio, occurrence_cutoff):
		"""
		10-03-05
			fix a bug related to cluster size
		"""
		row = cluster_queue.get()
		while row != -1:
			mcl_id, splat_id, vertex_set, recurrence_array, connectivity = row
			if mcl_id in mcl_id2unknown_ratio:
				unknown_ratio, prediction_tuple_set = mcl_id2unknown_ratio[mcl_id]
				vertex_list = vertex_set[1:-1].split(',')	#10-03-05 this is a list, to calculate the correct cluster size, len(vertex_list)
				recurrence_array = recurrence_array[1:-1].split(',')
				recurrence_array = map(float, recurrence_array)
				for i in range(len(recurrence_array)):
					if recurrence_array[i] >=occurrence_cutoff:
						recurrence_array[i] = 1
					else:
						recurrence_array[i] = 0
				recurrence = sum(recurrence_array)
				go_no_list = []
				p_value_list = []
				for prediction_tuple in prediction_tuple_set:
					go_no_list.append(prediction_tuple[0])
					p_value_list.append(prediction_tuple[1])
				curs.execute("insert into %s(mcl_id, splat_id, vertex_set, recurrence_array, recurrence,\
					connectivity, unknown_ratio, size, go_no_list, p_value_list) values(\
					%s, %s, '%s', '{%s}', %s, %s, %s, %s, '{%s}', '{%s}')"%\
					(good_cluster_table, mcl_id, splat_id, vertex_set, repr(recurrence_array)[1:-1], recurrence,\
					connectivity, unknown_ratio, len(vertex_list), repr(go_no_list)[1:-1], repr(p_value_list)[1:-1]))
			row = cluster_queue.get()
		
	def run(self):
		"""
		10-29-05 call p_gene_id_set_from_gene_p_table()
		"""
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		if self.new_table:
			self.create_good_cluster_table(curs, self.good_cluster_table)
		p_gene_id_set = p_gene_id_set_from_gene_p_table(curs, self.gene_p_table)
		mcl_id2unknown_ratio = self.get_mcl_id2unknown_ratio(curs, self.p_gene_table, p_gene_id_set)
		self.submit_good_clusters(curs, self.cluster_queue, self.good_cluster_table, mcl_id2unknown_ratio, self.occurrence_cutoff)
		if self.commit:
			curs.execute("end")
	
class fetch_clusters(Thread):
	def __init__(self, cluster_queue, hostname='zhoudb', dbname='graphdb', schema=None, gene_p_table=None,\
		p_gene_table=None, mcl_table=None, good_cluster_table=None,  occurrence_cutoff=0.8, \
		new_table=0, commit=0, debug=0, report=0):
		Thread.__init__(self)
		self.cluster_queue = cluster_queue
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.gene_p_table = gene_p_table
		self.p_gene_table = p_gene_table
		self.mcl_table = mcl_table
		self.good_cluster_table = good_cluster_table
		self.occurrence_cutoff = float(occurrence_cutoff)
		self.new_table = int(new_table)
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
	
	def fetch_good_clusters(self, curs, mcl_table, cluster_queue):
		sys.stderr.write("Fetching the contents of good clusters...\n")
		curs.execute("DECLARE crs CURSOR FOR select mcl_id, splat_id, vertex_set, recurrence_array, connectivity\
			from %s"%mcl_table)
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		counter = 0
		while rows:
			for row in rows:
				cluster_queue.put(row)
					
				counter += 1
			if self.report:
				sys.stderr.write("%s%s qs:%s"%('\x08'*20, counter,cluster_queue.qsize()))
			if self.debug and counter>=200000:
				break
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		cluster_queue.put(-1)	#The exit signal
	
	def run(self):
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		self.fetch_good_clusters(curs, self.mcl_table, self.cluster_queue)
	

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:g:p:m:l:s:o:ncbr", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	gene_p_table = None
	p_gene_table = None
	mcl_table = None
	good_cluster_table = None
	occurrence_cutoff = 0.8
	qsize = 400000
	new_table = 0
	commit = 0
	debug = 0
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
		elif opt in ("-g",):
			gene_p_table = arg
		elif opt in ("-p",):
			p_gene_table = arg
		elif opt in ("-m",):
			mcl_table = arg
		elif opt in ("-l",):
			good_cluster_table = arg
		elif opt in ("-o",):
			occurrence_cutoff = float(arg)
		elif opt in ("-s",):
			qsize = int(arg)
		elif opt in ("-n",):
			new_table = 1
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-b",):
			debug = 1
		elif opt in ("-r",):
			report = 1
	if schema and gene_p_table and p_gene_table and mcl_table and good_cluster_table:
		q = Queue(qsize)
		instance1 = submit_clusters(q, hostname, dbname, schema, gene_p_table, p_gene_table,\
			mcl_table, good_cluster_table, occurrence_cutoff, new_table, commit, debug, report)
		instance2 = fetch_clusters(q, hostname, dbname, schema, gene_p_table, p_gene_table,\
			mcl_table, good_cluster_table, occurrence_cutoff, new_table, commit, debug, report)
		instance1.start()
		instance2.start()
		instance1.join()
		instance2.join()
	else:
		print __doc__
		sys.exit(2)
