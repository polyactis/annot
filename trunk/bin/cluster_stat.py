#!/usr/bin/env python
"""
Usage: cluster_stat.py -k SCHEMA [OPTION]

Option:
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-w, --wu	apply Wu's strategy(Default is Jasmine's strategy)
	-c, --commit	commit the database transaction
	-r, --report	report the progress(a number)
	-h, --help              show this help
	
Examples:
	cluster_stat.py -k shu -c
	cluster_stat.py -k ming1
	
Description:
	Program for computing cluster p-value vectors(leave-one-out).
	The length of p_value_vector is subject to the number of functional categories.
	The difference between two strategies:
	validatoin:	Jasmine's only works on known-gene clusters.
		Wu's works on whole-gene clusters.
	prediciton:  Jasmine's only works on clusters with only one unknown gene.
		Wu's works on all clusters with any number of unknown genes and needs
		a cut_off for unknown functional class.
"""

import sys,os,psycopg,pickle,getopt
from rpy import r

class cluster_stat:
	def __init__(self, dbname, schema, report=0, wu=0, needcommit=0):
		self.conn = psycopg.connect('dbname=%s'%dbname)
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.report = int(report)
		self.wu = int(wu)
		self.needcommit = int(needcommit)
		self.global_go_id_to_no_dict = {}
		self.global_go_no_to_size_dict = {}
		self.global_gene_to_go_dict = {}
		self.no_of_records = 0
		self.logfile = open('/tmp/cluster_stat.log','w')

	def dstruc_loadin(self):
		if self.wu:
			self.curs.execute("select go_id,go_no,array_upper(gene_array,1) from go")
		else:
			self.curs.execute("select go_id,go_no,array_upper(gene_array,1) from go where go_no!=0")
		rows = self.curs.fetchall()
		for row in rows:
			self.global_go_id_to_no_dict[row[0]] = row[1]
			self.global_go_no_to_size_dict[row[1]] = row[2]
		self.no_of_functions = len(self.global_go_no_to_size_dict)
		if self.wu:
			self.curs.execute("select gene_no,go_functions from gene")
		else:
			self.curs.execute("select gene_no,go_functions from gene where known=TRUE")
			
		rows = self.curs.fetchall()
		for row in rows:
			self.global_gene_to_go_dict[row[0]] = []
			go_functions_list = row[1][1:-1].split(',')
			for go_no in go_functions_list:
				self.global_gene_to_go_dict[row[0]].append(int(go_no))
		self.no_of_genes = len(self.global_gene_to_go_dict)
		
	def run(self):
		self.curs.execute("begin")
		self.curs.execute("DECLARE crs CURSOR FOR select mcl_id,vertex_set,connectivity from mcl_result")
		self.curs.execute("fetch 500 from crs")
		rows = self.curs.fetchall()
		while rows:
			for row in rows:
				mcl_id = row[0]
				vertex_set = row[1]
				connectivity = row[2]
				self._cluster_stat(mcl_id, vertex_set, connectivity)
				
			if self.report:
				sys.stderr.write("%s%s"%("\x08"*20,self.no_of_records))
			self.curs.execute("fetch 500 from crs")
			rows = self.curs.fetchall()
		if self.needcommit:
			self.curs.execute("end")
			sys.stderr.write('\n\tTotal %d records.\n'%self.no_of_records)
		else:
			self.conn.rollback()
			sys.stderr.write('\n\tNo real updates\n')
		
	def _cluster_stat(self, mcl_id, vertex_set, connectivity):
		vertex_list_all = vertex_set[1:-1].split(',')
		vertex_list = []
		for vertex in vertex_list_all:
			if int(vertex) in self.global_gene_to_go_dict:
			#this filter will only be useful when Jasmine's strategy is applied to whole gene-set(unknown included)
				vertex_list.append(int(vertex))
		cluster_size = len(vertex_list)
		p_value_vector = [1] * self.no_of_functions
		self.local_go_no_dict_construct(vertex_list)
		for gene_no in vertex_list:
			self.go_no_dict_adjust(gene_no)
			for go_no in self._local_go_no_dict:
				x = self._local_go_no_dict[go_no]
				m = self._global_go_no_dict[go_no]
				n = self.no_of_genes -1 - m
				k = cluster_size-1
				p_value = r.phyper(x-1,m,n,k,lower_tail = r.FALSE)
				self.logfile.write('%d %d %d %d %d %d %d %f\n'%\
					(mcl_id,gene_no,go_no,x,m,n,k,p_value))
				if self.wu:
					p_value_vector[go_no] = p_value
				else:
					p_value_vector[go_no-1] = p_value
			if self.needcommit:
				self.curs.execute("insert into cluster_stat(mcl_id, leave_one_out, p_value_vector, connectivity)\
				values(%d, %d, ARRAY%s, %8.6f)"%(mcl_id, gene_no, repr(p_value_vector), connectivity))
			self.no_of_records += 1

	def local_go_no_dict_construct(self, vertex_list):
		'''
		construct a local go_no:size dictionary for a specific cluster.
		'''
		self.local_go_no_dict = {}
		for gene_no in vertex_list:
			go_no_list = self.global_gene_to_go_dict[gene_no]
			for go_no in go_no_list:
				if go_no in self.local_go_no_dict:
					self.local_go_no_dict[go_no] += 1
				else:
					self.local_go_no_dict[go_no] = 1

		
	def go_no_dict_adjust(self, gene_no):
		'''
		After one gene is left out, both local and global go_no:size dictionaries
		need to be adjusted.
		'''
		self._local_go_no_dict = self.local_go_no_dict.copy()
		self._global_go_no_dict = self.global_go_no_to_size_dict.copy()
		go_no_list = self.global_gene_to_go_dict[gene_no]
		for go_no in go_no_list:
			self._local_go_no_dict[go_no] -= 1
			self._global_go_no_dict[go_no] -= 1
			if self._local_go_no_dict[go_no] == 0:
				del self._local_go_no_dict[go_no]
			if self._global_go_no_dict[go_no] == 0:
				del self._global_go_no_dict[go_no]


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hrd:k:cw", ["help", "report", "dbname=", "schema=", "commit", "wu"])
	except:
		print __doc__
		sys.exit(2)
	
	dbname = 'graphdb'
	schema = ''
	commit = 0
	report = 0
	wu = 0
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-k", "--schema"):
			schema = arg
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-w", "--wu"):
			wu = 1

	if schema:
		instance = cluster_stat(dbname, schema, report, wu, commit)
		instance.dstruc_loadin()
		instance.run()

	else:
		print __doc__
		sys.exit(2)
