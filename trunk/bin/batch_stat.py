#!/usr/bin/env python
"""
Usage: batch_stat.py -k SCHEMA [OPTION]

Option:
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --tag=...		a tag to identify all the results, required if commit
	-w, --wu	Wu's strategy(Default is Jasmine's strategy)
	-c, --commit	commit the database transaction, records in table stat_plot_data
	-h, --help              show this help
	
Examples:
	batch_stat.py -k shu
	batch_stat.py -k shu -t sc_known -c
"""

import sys, os, cStringIO, psycopg, getopt
from gene_stat import gene_stat

class batch_stat:
	def __init__(self, dbname, schema, tag, wu=0, needcommit=0):
		self.dbname = dbname
		self.conn = psycopg.connect('dbname=%s'%self.dbname)
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.connectivity_list = [0.7, 0.8, 0.9]
		self.limit_list = [0]
		self.p_value_cut_off_list = [0.001, 0.0005, 0.0001, 0.00005, 0.00001]
		self.stat_data = []
		self.tag = tag
		self.wu = int(wu)
		self.needcommit = int(needcommit)
		
	def run(self):
		for connectivity in self.connectivity_list:
			for limit in self.limit_list:
				for p_value_cut_off in self.p_value_cut_off_list:
					if limit == 0 and p_value_cut_off == 0.001 and connectivity == 0.8:
						instance = gene_stat(self.dbname, schema, p_value_cut_off, limit, connectivity, self.wu, 0, 1)
					else:
						instance = gene_stat(self.dbname, schema, p_value_cut_off, limit, connectivity, self.wu)
					instance.dstruc_loadin()
					instance.run()
					self.stat_data.append([connectivity, limit, p_value_cut_off, \
					 instance.tp, instance.tn, instance.fp, instance.fn])
					del instance
		self.submit()
			
	def submit(self):
		for item in self.stat_data:
			self.curs.execute("insert into graph.stat_plot_data(connectivity, limit_of_cluster, p_value_cut_off,\
			tp, tn, fp, fn, tag) values (%1.2f, %d, %1.5f, %d, %d, %d, %d, '%s')"%\
			(item[0], item[1], item[2], item[3], item[4], item[5], item[6], self.tag))
		if self.needcommit:
			self.conn.commit()
		
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hd:k:t:wc", ["help", "dbname=", "schema=", "tag=", "wu", "commit"])
	except:
		print __doc__
		sys.exit(2)
	
	dbname = 'graphdb'
	schema = ''
	tag = ''
	wu =0 
	commit = 0
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-k", "--schema"):
			schema = arg
		elif opt in ("-t", "--tag"):
			tag = arg
		elif opt in ("-w", "--wu"):
			wu = 1
		elif opt in ("-c", "--commit"):
			commit = 1
	
	if schema:
		if commit == 1 and tag == '':
			print __doc__
			sys.exit(2)
		instance = batch_stat(dbname, schema, tag, wu, commit)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
