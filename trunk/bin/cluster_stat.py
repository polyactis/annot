#!/usr/bin/env python
import sys,os,psycopg,pickle
from rpy import r

class cluster_stat:
	def __init__(self, dbname,organism, needcommit=0):
		self.conn = psycopg.connect('dbname=%s'%dbname)
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to graph")
		self.needcommit = int(needcommit)
		self.organism = organism
		self.global_go_id_to_no_dict = {}
		self.global_go_no_to_size_dict = {}
		self.global_gene_to_go_dict = {}
		self.no_of_records = 0
		self.logfile = open('/tmp/cluster_stat.log','w')

	def dstruc_loadin(self):
		self.curs.execute("select go_id,go_no,array_upper(gene_array,1) from graph.sc_go")
		rows = self.curs.fetchall()
		for row in rows:
			self.global_go_id_to_no_dict[row[0]] = row[1]
			self.global_go_no_to_size_dict[row[1]] = row[2]
		self.no_of_functions = len(self.global_go_no_to_size_dict)
		self.curs.execute("select gene_no,go_functions from graph.sc_gene where known=TRUE")
		rows = self.curs.fetchall()
		for row in rows:
			self.global_gene_to_go_dict[row[0]] = []
			go_functions_list = row[1][1:-1].split(',')
			for go_no in go_functions_list:
				self.global_gene_to_go_dict[row[0]].append(int(go_no))
		self.no_of_genes = len(self.global_gene_to_go_dict)
		
	def run(self):
		self.curs.execute("begin")
		self.curs.execute("DECLARE crs CURSOR FOR select mcl_id,vertex_set from mcl_result")
		self.curs.execute("fetch 500 from crs")
		rows = self.curs.fetchall()
		while rows:
			for row in rows:
				mcl_id = row[0]
				vertex_set = row[1]
				self._cluster_stat(mcl_id, vertex_set)
				
			sys.stderr.write("%s%s"%("\x08"*20,self.no_of_records))
			self.curs.execute("fetch 500 from crs")
			rows = self.curs.fetchall()
		if self.needcommit:
			self.curs.execute("end")
			sys.stderr.write('\n\tTotal %d records.\n'%self.no_of_records)
		else:
			self.conn.rollback()
			sys.stderr.write('\n\tNo real updates\n')
		
	def _cluster_stat(self, mcl_id, vertex_set):
		vertex_list = vertex_set[1:-1].split(',')
		cluster_size = len(vertex_list)
		for i in range(cluster_size):
			vertex_list[i] = int(vertex_list[i])
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
				p_value_vector[go_no-1] = p_value
			if self.needcommit:
				self.curs.execute("insert into sc_cluster_stat(mcl_id, leave_one_out, p_value_vector)\
				values(%d, %d, ARRAY%s)"%(mcl_id, gene_no, repr(p_value_vector)))
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
	def helper():
		sys.stderr.write('\
	argv[1] is the database name\n\
	argv[2] is two-letter abbreviation for an organism.\n\
	argv[3] is 1 or 0 indicating whether to commit or not. Default is 0.\n')
	
	if len(sys.argv) == 3:
		instance = cluster_stat(sys.argv[1], sys.argv[2])
		instance.dstruc_loadin()
		instance.run()
	elif len(sys.argv) == 4:
		instance = cluster_stat(sys.argv[1], sys.argv[2],sys.argv[3])
		instance.dstruc_loadin()
		instance.run()
	else:
		helper()
