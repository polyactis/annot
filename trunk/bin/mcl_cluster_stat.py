#!/usr/bin/env python
"""
Usage: mcl_cluster_stat.py -k SCHEMA [OPTIONS] OFNAME

Option:
	OFNAME is the file to store the stat table
	-z ..., --hostname=...	the hostname, zhoudb (default)
	-d ..., --dbname=...	the database name, graphdb (default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	mcl_result (default)
	-p ..., --p_value_cut_off=...	0.001 (default)
	-n ..., --connectivity_cut_off=...	0.8 (default), minimum connectivity of a mcl cluster
	-r, --report	report the progress (a number)
	-h, --help              show this help

Examples:
	mcl_cluster_stat.py -k sc_yh60_splat -t mcl_result_sup_3_2 mcl_cluster_stat

Description:
	A program to dissect the composition of mcl clusters.
"""

import sys, os, psycopg, getopt, csv
from sets import Set

class mcl_cluster_struc:
	def __init__(self):
		self.gene_set = Set()
		self.cluster_size_list = []
		self.connectivity_list = []
		self.p_value_list = []
		self.no_of_good_clusters = 0

class mcl_cluster_stat:
	def __init__(self, hostname, dbname, schema, table, p_value_cut_off, \
		connectivity_cut_off, report, ofname):
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.table = table
		self.p_value_cut_off = float(p_value_cut_off)
		self.connectivity_cut_off = float(connectivity_cut_off)
		self.report = int(report)
		self.writer = csv.writer(open(ofname, 'w'), delimiter='\t')
		
		self.mcl_struc_dict = {}
		self.no_of_records = 0
		
	def dstruc_loadin(self):
		sys.stderr.write("Loading Data STructure...")

		sys.stderr.write("Done\n")
		
	def run(self):
		self.curs.execute("DECLARE crs CURSOR FOR select vertex_set, connectivity, p_value_min, \
				array_upper(recurrence_array, 1) from %s where p_value_min notnull"%(self.table))
		
		self.curs.execute("fetch 5000 from crs")
		rows = self.curs.fetchall()
		while rows:
			for row in rows:
				self._stat(row)
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, self.no_of_records))
			
			self.curs.execute("fetch 5000 from crs")
			rows = self.curs.fetchall()
		
		self.output()

	def _stat(self, row):
		self.no_of_records += 1
		vertex_set = row[0][1:-1].split(',')
		vertex_set = map(int, vertex_set)
		connectivity = row[1]
		p_value_min = row[2]
		recurrence = row[3]
		if recurrence not in self.mcl_struc_dict:
			self.mcl_struc_dict[recurrence] = mcl_cluster_struc()
		for gene_no in vertex_set:
			self.mcl_struc_dict[recurrence].gene_set.add(gene_no)
		self.mcl_struc_dict[recurrence].cluster_size_list.append(len(vertex_set))
		self.mcl_struc_dict[recurrence].connectivity_list.append(connectivity)
		self.mcl_struc_dict[recurrence].p_value_list.append(p_value_min)
		if p_value_min <= self.p_value_cut_off and connectivity >= self.connectivity_cut_off:
			self.mcl_struc_dict[recurrence].no_of_good_clusters += 1

	def output(self):
		recurrence_list = self.mcl_struc_dict.keys()
		recurrence_list.sort()
		rows = {}
		rows[0] = ['recurrence']
		rows[1] = ['number of clusters']
		rows[2] = ['number of genes']
		rows[3] = ['avg cluster_size']
		rows[4] = ['avg p_value']
		rows[5] = ['avg connectivity']
		rows[6] = ['good clusters']
		for recurrence in self.mcl_struc_dict:
			unit = self.mcl_struc_dict[recurrence]
			rows[0].append(recurrence)
			no_of_clusters = len(unit.cluster_size_list)
			rows[1].append(no_of_clusters)
			rows[2].append(len(unit.gene_set))
			rows[3].append('%3.2f'%(sum(unit.cluster_size_list)/float(no_of_clusters)))
			rows[4].append('%1.6f'%(sum(unit.p_value_list)/float(no_of_clusters)))
			rows[5].append('%1.3f'%(sum(unit.connectivity_list)/float(no_of_clusters)))
			rows[6].append('%s(%2.3f%%)'%(unit.no_of_good_clusters, unit.no_of_good_clusters*100/float(no_of_clusters)))
		for i in range(7):
			self.writer.writerow(rows[i])
		


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema=", "table=", "p_value_cut_off=",\
		"connectivity_cut_off=", "report"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:p:n:r", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = 'mcl_result'
	p_value_cut_off = 0.001
	connectivity_cut_off = 0.8
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
		elif opt in ("-p", "--p_value_cut_off"):
			p_value_cut_off = float(arg)
		elif opt in ("-n", "--connectivity_cut_off"):
			connectivity_cut_off = float(arg)
		elif opt in ("-r", "--report"):
			report = 1
		
	if schema and len(args) == 1:
		instance = mcl_cluster_stat(hostname, dbname, schema, table, p_value_cut_off,\
			connectivity_cut_off, report, args[0])
		instance.dstruc_loadin()
		instance.run()
	else:
		print __doc__
		sys.exit(2)
