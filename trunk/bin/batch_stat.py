#!/usr/bin/env python
"""
Usage: batch_stat.py -d DATABASENAME [OPTION]

Option:
	-d ..., --dbname=...	the database name
	-c ..., --commit=...	1 or 0(default) specifies commit or not
	-t ..., --tag=...		a tag to identify all the results, required if commit
	-h, --help              show this help
	
Examples:
	batch_stat.py -d mdb
	batch_stat.py -d mdb -t sc_known -c 1
"""

import sys, os, cStringIO, psycopg, getopt
from gene_stat import gene_stat

class batch_stat:
	def __init__(self, dbname, tag, needcommit=0):
		self.dbname = dbname
		self.conn = psycopg.connect('dbname=%s'%self.dbname)
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to graph")
		self.connectivity_list = [0.8,0.9]
		self.limit_list = [5,10,15,20]
		self.p_value_cut_off_list = [0.1, 0.01, 0.001, 0.0005, 0.0001, 0.00001]
		self.stat_data = []
		self.tag = tag
		self.needcommit = int(needcommit)
		
	def run(self):
		for connectivity in self.connectivity_list:
			for limit in self.limit_list:
				for p_value_cut_off in self.p_value_cut_off_list:
					instance = gene_stat(self.dbname, p_value_cut_off, limit, connectivity)
					instance.dstruc_loadin()
					instance.run()
					self.stat_data.append([connectivity, limit, p_value_cut_off, \
					 instance.tp, instance.tn, instance.fp, instance.fn])
					del instance
		if self.needcommit:
			self.submit()
			
	def submit(self):
		for item in self.stat_data:
			curs.execute("insert into stat_plot_data(connectivity, limit_of_cluster, p_value_cut_off,\
			tp, tn, fp, fn, tag) values (%1.2f, %d, %1.5f, %d, %d, %d, %d, '%s')"%\
			(item[0], item[1], item[2], item[3], item[4], item[5], item[6], self.tag))
	
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hd:t:c:", ["help", "dbname=", "tag=", "commit="])
	except:
		print __doc__
		sys.exit(2)
	
	dbname = ''
	commit = 0
	tag = ''
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-t", "--tag"):
			tag = arg
		elif opt in ("-c", "--commit"):
			commit = int(arg)
	
	if dbname:
		if commit == 1 and tag == '':
			print __doc__
			sys.exit(2)
		instance = batch_stat(dbname, tag, commit)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
