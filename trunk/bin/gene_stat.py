#!/usr/bin/env python

import sys, os,psycopg
from rpy import *
import numarray.ma as ma


class gene_stat:
	def __init__(self, dbname, p_value_cut_off, connectivity_cut_off=0.9,needcommit=0):
		self.conn = psycopg.connect('dbname=%s'%dbname)
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to graph")
		self.needcommit = int(needcommit)
		self.p_value_cut_off = float(p_value_cut_off)
		self.connectivity_cut_off = float(connectivity_cut_off)
		self.tp = 0.0
		self.tn = 0.0
		self.fp = 0.0
		self.fn = 0.0
		self.known_genes_dict = {}
		self.no_of_records = 0
		
	def dstruc_loadin(self):
		self.curs.execute("select gene_no,go_functions from sc_gene where known=TRUE")
		rows = self.curs.fetchall()
		for row in rows:
			self.known_genes_dict[row[0]] = row[1]
		
	def run(self):
		for gene_no in self.known_genes_dict:
			mcl_id_list = []
			self.curs.execute("select c.mcl_id, c.p_value_vector \
				from sc_cluster_stat c, mcl_result m where c.leave_one_out=%d \
				and c.mcl_id=m.mcl_id and m.connectivity>=%f"%(gene_no, self.connectivity_cut_off))
			rows = self.curs.fetchall()
			init_mask = ma.zeros(83)
			for row in rows:
				mcl_id_list.append(row[0])
				p_value_vector = row[1][1:-1].split(',')
				for i in range(83):
					if float(p_value_vector[i]) <= self.p_value_cut_off:
						init_mask[i] = 1
			self._gene_stat(gene_no, mcl_id_list, init_mask)
			self.no_of_records += 1
			sys.stderr.write('%s%s'%('\x08'*20, self.no_of_records))
		if self.needcommit:
			self.conn.commit()
		sys.stderr.write('\n\tTotal genes: %d\n'%self.no_of_records)
		sys.stderr.write('\tSensitvity: %f\n'%(self.tp/(self.tp+self.fn)))
		sys.stderr.write('\tSpecificity: %f\n'%(self.tn/(self.fp+self.tn)))
			
	def _gene_stat(self, gene_no, mcl_id_list, init_mask):
		tp = 0
		tn = 0
		fp = 0
		fn = 0
		go_functions_vector = self.known_genes_dict[gene_no]
		for i in range(83):
			if go_functions_vector[i] == '0':
				if init_mask[i] == 0:
					tn += 1
					self.tn += 1
				else:
					fp += 1
					self.fp += 1
			else:
				if init_mask[i] == 0:
					fn += 1
					self.fn += 1
				else:
					tp += 1
					self.tp += 1
		if self.needcommit:
			self.curs.execute("insert into sc_gene(gene_no,cluster_array,tp,tn,fp,fn)\
			values(%d, ARRAY%s, %d,%d,%d,%d)"%(gene_no,repr(mcl_id_list),tp,tn,fp,fn))


if __name__ == '__main__':
	def helper():
		sys.stderr.write('\
	argv[1] is the database name\n\
	argv[2] is p_value_cut_off.\n\
	argv[3] is the connectivity_cut_off.Default is 0.9.\n\
	argv[4] is 1 or 0 indicating whether to commit or not. Default is 0.\n')
	
	if len(sys.argv) == 3:
		instance = gene_stat(sys.argv[1], sys.argv[2])
		instance.dstruc_loadin()
		instance.run()
	elif len(sys.argv) == 4:
		instance = gene_stat(sys.argv[1], sys.argv[2],sys.argv[3])
		instance.dstruc_loadin()
		instance.run()
	elif len(sys.argv) == 5:
		instance = gene_stat(sys.argv[1], sys.argv[2],sys.argv[3],sys.argv[4])
		instance.dstruc_loadin()
		instance.run()
	else:
		helper()
