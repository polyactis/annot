#!/usr/bin/env python
"""
Usage: mcl_result_stat.py -k SCHEMA [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	mcl_result(default) or fim_result
	-b, --bonferroni	bonferroni correction(default, ignore it)
	-w, --wu	apply Wu's strategy(default, ignore it)
	-c, --commit	commit the database transaction
	-r, --report	report the progress(a number)
	-h, --help              show this help
	
Examples:
	mcl_result_stat.py -k shu -c
	mcl_result_stat.py -k ming1
	mcl_result_stat.py -k shu_whole -c -t fim_result

Description:
	Program for computing cluster p-value vectors.
"""

import sys,os,psycopg,pickle,getopt
from rpy import r

class mcl_result_stat:
	def __init__(self, hostname, dbname, schema, table, bonferroni=0, report=0, wu=1, needcommit=0):
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.table = table
		self.bonferroni = int(bonferroni)
		self.report = int(report)
		self.wu = int(wu)
		self.needcommit = int(needcommit)
		self.global_go_id_to_no_dict = {}
		self.global_go_no_to_size_dict = {}
		self.global_gene_to_go_dict = {}
		self.no_of_records = 0
		self.logfile = open('/tmp/mcl_result_stat.log','w')
		self.cluster_memory = {}
		self.to_db = []
		
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
		self.curs.execute("DECLARE crs CURSOR FOR select mcl_id,vertex_set from %s"%self.table)
		self.curs.execute("fetch 5000 from crs")
		rows = self.curs.fetchall()
		while rows:
			for row in rows:
				mcl_id = row[0]
				vertex_set = row[1]
				self._cluster_stat(mcl_id, vertex_set)
				
			if self.report:
				sys.stderr.write("%s%s"%("\x08"*20,self.no_of_records))
			self.curs.execute("fetch 5000 from crs")
			rows = self.curs.fetchall()
		if self.needcommit:
			sys.stderr.write("\nDatabase updating...")
			for item in self.to_db:
				self.curs.execute("update %s set p_value_min=%f, go_no_vector=ARRAY%s, unknown_gene_ratio=%f\
					where mcl_id=%d"%(self.table, item[0], item[1], item[2], item[3]))
			self.curs.execute("end")
			sys.stderr.write("Done.\n")
			sys.stderr.write('\n\tTotal %d records.\n'%self.no_of_records)
		else:
			self.conn.rollback()
			sys.stderr.write('\n\tNo real updates\n')
		
	def _cluster_stat(self, mcl_id, vertex_set):
		if vertex_set in self.cluster_memory:
			entry = self.cluster_memory[vertex_set]
			p_value_min = entry[0]
			go_no_vector = entry[1]
			unknown_gene_ratio = entry[2]
			self.to_db.append([p_value_min, go_no_vector, unknown_gene_ratio, mcl_id])
			
			self.no_of_records += 1
			return
		else:
			_cluster_memroy = {}
		vertex_list = vertex_set[1:-1].split(',')
		vertex_list = map(int, vertex_list)
		cluster_size = len(vertex_list)
		self.local_go_no_dict_construct(vertex_list)
		if 0 in self.local_go_no_dict:
			unknown_gene_ratio = self.local_go_no_dict[0]/float(cluster_size)
		else:
			unknown_gene_ratio = 0
		if self.local_go_no_dict == {}:
			self.logfile.write('%d %s: local_go_no_dict empty\n'%(mcl_id, repr(vertex_set)))
			return
		for go_no in self.local_go_no_dict:
			if self.wu:
			# code after 'or' deals with the situation that Jasmine's strategy is applied to whole gene-set(unknown included)
				x = self.local_go_no_dict[go_no]
				m = self.global_go_no_to_size_dict[go_no]
				n = self.no_of_genes - m
				k = cluster_size
			else:
				pass
			if self.bonferroni:
				p_value = r.phyper(x-1,m,n,k,lower_tail = r.FALSE)*len(self.local_go_no_dict)
			else:
				p_value = r.phyper(x-1,m,n,k,lower_tail = r.FALSE)
			self.logfile.write('%d %d %d %d %d %d %f %f\n'%(mcl_id,go_no,x,m,n,k,p_value, unknown_gene_ratio))
			if p_value in _cluster_memroy:
				_cluster_memroy[p_value].append(go_no)
			else:
				_cluster_memroy[p_value] = [go_no]
		p_value_vector = _cluster_memroy.keys()
		if p_value_vector == []:
			self.logfile.write('%d %s: all vertices belong to population singleton classes\n'%(mcl_id, repr(vertex_set)))
			return
		p_value_min = min(p_value_vector)
		go_no_vector = _cluster_memroy[p_value_min]

		self.no_of_records += 1
		self.cluster_memory[vertex_set] =[p_value_min, go_no_vector, unknown_gene_ratio]
		self.to_db.append([p_value_min, go_no_vector, unknown_gene_ratio, mcl_id])
		
	def local_go_no_dict_construct(self, vertex_list):
		'''
		construct a local go_no:size dictionary for a specific cluster.
		'''
		self.local_go_no_dict = {}
		for gene_no in vertex_list:
			go_no_list = self.global_gene_to_go_dict[gene_no]
			for go_no in go_no_list:
				if self.global_go_no_to_size_dict[go_no] > 1:
				#population singleton is discarded
					if go_no in self.local_go_no_dict:
						self.local_go_no_dict[go_no] += 1
					else:
						self.local_go_no_dict[go_no] = 1



if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hrz:d:k:t:bcw", ["help", "report", "hostname=", "dbname=", "table=", "schema=", "bonferroni", "commit", "wu"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = 'mcl_result'
	bonferroni = 1
	commit = 0
	report = 0
	wu = 1
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
		elif opt in ("-b", "--bonferroni"):
			bonferroni = 1
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-w", "--wu"):
			wu = 1

	if schema:
		instance = mcl_result_stat(hostname, dbname, schema, table, bonferroni, report, wu, commit)
		instance.dstruc_loadin()
		instance.run()

	else:
		print __doc__
		sys.exit(2)
