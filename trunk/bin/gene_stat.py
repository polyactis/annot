#!/usr/bin/env python

import sys, os,psycopg
from rpy import *
#import numarray.ma as ma


class gene_stat:
	def __init__(self, dbname, p_value_cut_off, limit=10, connectivity_cut_off=0.9,needcommit=0):
		self.conn = psycopg.connect('dbname=%s'%dbname)
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to graph")
		self.limit = int(limit)
		self.needcommit = int(needcommit)
		self.p_value_cut_off = float(p_value_cut_off)
		self.connectivity_cut_off = float(connectivity_cut_off)
		self.tp = 0.0
		self.tn = 0.0
		self.fp = 0.0
		self.fn = 0.0
		self.known_genes_dict = {}
		self.no_of_records = 0
		self.sc_gene_to_be_updated = []
		self.log_file = open('/tmp/gene_stat.log','w')
		
	def dstruc_loadin(self):
		self.curs.execute("select gene_no,go_functions from sc_gene where known=TRUE")
		rows = self.curs.fetchall()
		for row in rows:
			go_functions_list = row[1][1:-1].split(',')
			self.known_genes_dict[row[0]] = []
			for go_no in go_functions_list:
				self.known_genes_dict[row[0]].append(int(go_no))
	def run(self):
		for gene_no in self.known_genes_dict:
			mcl_id_list = []
			self.log_file.write('%d '%gene_no)
			self.curs.execute("select c.mcl_id, c.p_value_vector \
				from sc_cluster_stat c, mcl_result m where c.leave_one_out=%d \
				and c.mcl_id=m.mcl_id and m.connectivity>=%f limit %d"%(gene_no, self.connectivity_cut_off,self.limit))
			rows = self.curs.fetchall()
			p_functions_dict = {}
			for row in rows:
				mcl_id_list.append(row[0])
				p_value_vector = row[1][1:-1].split(',')
				for i in range(83):
					if float(p_value_vector[i]) <= self.p_value_cut_off:
						p_functions_dict[i+1] = 1
			if rows:
				self._gene_stat(gene_no, mcl_id_list, p_functions_dict)
				self.no_of_records += 1
				sys.stderr.write('%s%s'%('\x08'*10, self.no_of_records))
		if self.needcommit:
			for entry in self.sc_gene_to_be_updated:
				if entry[6]:
					self.curs.execute("update sc_gene set cluster_array=ARRAY%s,\
					tp=%d,tn=%d,fp=%d,fn=%d,p_functions=ARRAY%s where gene_no=%d"%\
					(repr(entry[1]),entry[2],entry[3],entry[4],entry[5],repr(entry[6]),entry[0]))
				else:
					self.curs.execute("update sc_gene set cluster_array=ARRAY%s,\
					tp=%d,tn=%d,fp=%d,fn=%d where gene_no=%d"%\
					(repr(entry[1]),entry[2],entry[3],entry[4],entry[5],entry[0]))
			self.conn.commit()
		sys.stderr.write('\n\tTotal genes: %d\n'%self.no_of_records)
		sys.stderr.write('\tSensitvity: %f\n'%(self.tp/(self.tp+self.fn)))
		sys.stderr.write('\tSpecificity: %f\n'%(self.tn/(self.fp+self.tn)))
		sys.stderr.write('\tFalse Positive Ratio: %f\n'%(self.fp/(self.tp+self.fp)))
			
	def _gene_stat(self, gene_no, mcl_id_list, p_functions_dict):
		tp = 0
		tn = 0
		fp = 0
		fn = 0
		p_functions = p_functions_dict.keys()
		for go_no in self.known_genes_dict[gene_no]:
			if go_no in p_functions_dict:
				tp += 1
				del p_functions_dict[go_no]
			else:
				fn += 1
		fp = len(p_functions_dict)
		tn = 83 - (tp+fp+fn)
		self.tp += tp
		self.tn += tn
		self.fp += fp
		self.fn += fn
		self.log_file.write('%d %d %d %d %s %s\n'%(tp,tn,fp,fn,repr(p_functions),repr(mcl_id_list)))
		self.sc_gene_to_be_updated.append([gene_no,mcl_id_list,tp,tn,fp,fn,p_functions])


if __name__ == '__main__':
	def helper():
		sys.stderr.write('\
	argv[1] is the database name\n\
	argv[2] is p_value_cut_off.\n\
	argv[3] specifies the maximum number of clusters related to one gene.\n\
		Default is 10.\n\
	argv[4] is the connectivity_cut_off.Default is 0.9.\n\
	argv[5] is 1 or 0 indicating whether to commit or not. Default is 0.\n')
	
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
	elif len(sys.argv) == 6:
		instance = gene_stat(sys.argv[1], sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
		instance.dstruc_loadin()
		instance.run()
	else:
		helper()
