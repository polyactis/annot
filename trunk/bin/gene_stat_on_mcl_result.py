#!/usr/bin/env python
"""
Usage: gene_stat_on_mcl_result.py -k SCHEMA -p P_VALUE_CUT_OFF [OPTION]

Option:
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-p ..., --p_value_cut_off=...	p_value_cut_off
	-u ..., --unknown_cut_off=...	unknown_gene_ratio_cut_off, 0.25(default)
	-l ..., --limit=...,	IGNORE it, interface relics.
	-n ..., --connectivity_cut_off=...	0.8(default), minimum connectivity of a mcl cluster
	-w, --wu	IGNORE it. interface relics.
	-r, --report	report the progress(a number)
	-c, --commit	commit the database transaction, records in table gene.
	-h, --help              show this help

Examples:
	gene_stat_on_mcl_result.py -k shu -p 0.001
	gene_stat_on_mcl_result.py -k shu -p 0.001 -n 0.7
	gene_stat_on_mcl_result.py -k shu -p 0.001 -n 0.7 -u 0.20 -r -c

Description:

"""

import sys, os, psycopg, getopt
from rpy import *

class gene_prediction:
	# class holding prediction information of a gene
	def __init__(self):
		self.tp = 0
		self.tp1 = 0
		self.tn = 0
		self.fp = 0
		self.fn = 0
		self.p_functions_dict = {}
		self.mcl_id_list = []

class gene_stat_on_mcl_result:
	def __init__(self, dbname, schema, p_value_cut_off, unknown_cut_off, limit=0, connectivity_cut_off=0.8, wu=0, report=0, needcommit=0):
		self.conn = psycopg.connect('dbname=%s'%dbname)
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.p_value_cut_off = float(p_value_cut_off)
		self.unknown_cut_off = float(unknown_cut_off)
		self.limit = int(limit)
		self.connectivity_cut_off = float(connectivity_cut_off)
		self.wu = int(wu)
		self.report = int(report)
		self.needcommit = int(needcommit)
		self.tp = 0.0
		self.tp1 = 0.0
		self.tn = 0.0
		self.fp = 0.0
		self.fn = 0.0
		self.known_genes_dict = {}
		self.no_of_records = 0
		self.log_file = open('/tmp/gene_stat_on_mcl_result.log','w')
		self.gene_prediction_dict = {}
		self.gono_goindex_dict = {}
		self.no_of_p_known = 0
		
	def dstruc_loadin(self):
		self.curs.execute("select gene_no,go_functions from gene where known=TRUE")
		rows = self.curs.fetchall()
		for row in rows:
			go_functions_list = row[1][1:-1].split(',')
			self.known_genes_dict[row[0]] = []
			for go_no in go_functions_list:
				self.known_genes_dict[row[0]].append(int(go_no))
		
		self.curs.execute("select count(go_no) from go")
		rows = self.curs.fetchall()
		self.no_of_functions = rows[0][0]
		
		self.curs.execute("select go_no, go_index from go")
		rows = self.curs.fetchall()
		for row in rows:
			go_index_list = row[1][2:-2].split('","')
			self.gono_goindex_dict[row[0]] = go_index_list
		
	def run(self):
		if self.limit != 0:
			sys.stderr.write("gene_stat_on_mcl_result now doesn't support limit\n")
			sys.exit(2)

		self.curs.execute("DECLARE crs CURSOR FOR select mcl_id, vertex_set, p_value_min, go_no_vector, unknown_gene_ratio\
				from mcl_result where connectivity>=%f and p_value_min notnull"%(self.connectivity_cut_off))
		self.curs.execute("fetch 5000 from crs")
		rows = self.curs.fetchall()
		while rows:
			for row in rows:
				#if self.known_genes_dict.has_key(row[1]):
				# for clusters containing unknown genes
				self._gene_stat(row)
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, self.no_of_records))
			
			self.curs.execute("fetch 5000 from crs")
			rows = self.curs.fetchall()
		
		self.final()
		if self.needcommit:
			self.curs.execute("select gene_no from gene")
			rows = self.curs.fetchall()
			for row in rows:
				gene_no = row[0]
				if gene_no in self.gene_prediction_dict:
					entry = self.gene_prediction_dict[gene_no]
					p_functions = entry.p_functions_dict.keys()
					if gene_no in self.known_genes_dict:
						if p_functions:
							self.curs.execute("update gene set cluster_array=ARRAY%s,\
							tp=%d,tp1=%d,tn=%d,fp=%d,fn=%d,p_functions=ARRAY%s where gene_no=%d"%\
							(repr(entry.mcl_id_list),entry.tp,entry.tp1,entry.tn,entry.fp,entry.fn,repr(p_functions),gene_no))
						else:
							self.curs.execute("update gene set cluster_array=ARRAY%s,\
							tp=%d,tp1=%d,tn=%d,fp=%d,fn=%d where gene_no=%d"%\
							(repr(entry.mcl_id_list),entry.tp,entry.tp1,entry.tn,entry.fp,entry.fn,gene_no))
					else:
						if p_functions:
							self.curs.execute("update gene set cluster_array=ARRAY%s,\
							p_functions=ARRAY%s where gene_no=%d"%\
							(repr(entry.mcl_id_list),repr(p_functions),gene_no))
						else:
							self.curs.execute("update gene set cluster_array=ARRAY%s \
							where gene_no=%d"%\
							(repr(entry.mcl_id_list),gene_no))
				else:
				#cleanup for other genes predicted previously
					self.curs.execute("update gene set cluster_array=null, tp=null, tp1=null, tn=null,\
						fp=null, fn=null, p_functions=null where gene_no=%d"%gene_no)
				
			self.curs.execute("end")
		sys.stderr.write('\n\tp_value_cut_off:%f unknown_cut_off:%f connectivity_cut_off:%f\n'%(self.p_value_cut_off, self.unknown_cut_off, self.connectivity_cut_off))
		sys.stderr.write('\tTotal genes: %d\n'%len(self.gene_prediction_dict))
		sys.stderr.write('\tTotal known genes: %d\n'%self.no_of_p_known)
		sys.stderr.write('\tWhole Sensitvity: %f\n'%((self.tp+self.tp1)/(self.tp+self.tp1+self.fn)))
		sys.stderr.write('\tTP0: %d  TP1: %d\n'%(self.tp, self.tp1))
		sys.stderr.write('\tSpecificity: %f\n'%(self.tn/(self.fp+self.tn)))
		sys.stderr.write('\tFalse Positive Ratio: %f\n'%(self.fp/(self.tp+self.tp1+self.fp)))
			
	def _gene_stat(self, row):
		mcl_id = row[0]
		vertex_set = row[1][1:-1].split(',')
		vertex_set = map(int, vertex_set)
		p_value_min = row[2]
		go_no_vector = row[3][1:-1].split(',')
		go_no_vector = map(int, go_no_vector)
		unknown_gene_ratio = row[4]
		if p_value_min>self.p_value_cut_off or unknown_gene_ratio>self.unknown_cut_off:
			return
		self.log_file.write('%d, %f, %f \n'%(mcl_id, p_value_min, unknown_gene_ratio))
		self.no_of_records += 1
		for gene_no in vertex_set:
			if gene_no not in self.gene_prediction_dict:
				item = gene_prediction()
				self.gene_prediction_dict[gene_no] = item
			self.gene_prediction_dict[gene_no].mcl_id_list.append(mcl_id)
			for go_no in go_no_vector:
				self.gene_prediction_dict[gene_no].p_functions_dict[go_no] = 1

	def final(self):
		for gene_no in self.gene_prediction_dict:
			L0_dict = {}
			L1_dict = {}
			good_p_functions_dict = {}
			entry = self.gene_prediction_dict[gene_no]
			p_functions_dict = entry.p_functions_dict.copy()
			p_functions = p_functions_dict.keys()
			if gene_no not in self.known_genes_dict:
				self.log_file.write('unknown: %d %s %s\n'%(gene_no, repr(p_functions),repr(entry.mcl_id_list)))
				continue
			k_functions_dict = self.known_genes_dict[gene_no]
			self.no_of_p_known += 1
			for p_go_no in p_functions_dict:
				for k_go_no in k_functions_dict:
					if self.is_L0(p_go_no, k_go_no):
						L0_dict[k_go_no] = 1
						good_p_functions_dict[p_go_no] = 1
					elif self.is_L1(p_go_no, k_go_no):
						L1_dict[k_go_no] = 1
						good_p_functions_dict[p_go_no] = 1
			for go_no in L0_dict:
				if go_no in L1_dict:
					del L1_dict[go_no]
			entry.tp = len(L0_dict)
			entry.tp1 = len(L1_dict)
			entry.fn = len(k_functions_dict) - entry.tp - entry.tp1
			entry.fp = len(p_functions_dict) - len(good_p_functions_dict)
			entry.tn = self.no_of_functions - (entry.tp+entry.tp1+entry.fp+entry.fn)

			self.tp += entry.tp
			self.tp1 += entry.tp1
			self.tn += entry.tn
			self.fp += entry.fp
			self.fn += entry.fn
			self.log_file.write('known: %d %d %d %d %d %d %s %s\n'%(gene_no, entry.tp,entry.tp1,entry.tn,entry.fp,entry.fn,repr(p_functions),repr(entry.mcl_id_list)))
	
	def is_L0(self, p_go_no, k_go_no):
		k_go_index_list = self.gono_goindex_dict[k_go_no]
		p_go_index_list = self.gono_goindex_dict[p_go_no]
		for k_index in k_go_index_list:
			for p_index in p_go_index_list:
				k_index_m = k_index + ','
				p_index_m = p_index + ','
				if p_index_m.find(k_index_m) == 0:
					self.log_file.write('%d is L0 of %d:: %s %s\n'%(p_go_no, k_go_no, p_index_m, k_index_m))
					return 1
		return 0
		
	def is_L1(self, p_go_no, k_go_no):
		k_go_index_list = self.gono_goindex_dict[k_go_no]
		p_go_index_list = self.gono_goindex_dict[p_go_no]
		for k_index in k_go_index_list:
			for p_index in p_go_index_list:
				pos = k_index.rfind(',')
				k_index_m = k_index[:pos]
				if k_index_m == p_index:
					self.log_file.write("%d is direct parent of %d:: %s %s\n"%(p_go_no, k_go_no, p_index, k_index_m))
					return 1
				pos = p_index.rfind(',')
				p_index_m = p_index[:pos]
				if k_index_m == p_index_m:
					self.log_file.write("%d and %d are siblings:: %s %s\n"%(p_go_no, k_go_no, p_index_m, k_index_m))
					return 1
		return 0
		
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "dbname=", "schema=", "p_value_cut_off=","unknown_cut_off=", "limit=", "connectivity_cut_off=", "wu", "report", "commit"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hd:k:p:u:l:n:wrc", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	dbname = 'graphdb'
	schema = ''
	p_value_cut_off = None
	limit = 0
	connectivity_cut_off = 0.8
	wu = 0
	report = 0
	commit = 0
	unknown_cut_off = 0.25
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-k", "--schema"):
			schema = arg
		elif opt in ("-p", "--p_value_cut_off"):
			p_value_cut_off = float(arg)
		elif opt in ("-u", "--unknown_cut_off"):
			unknown_cut_off = float(arg)
		elif opt in ("-l", "--limit"):
			limit = int(arg)
		elif opt in ("-n", "--connectivity_cut_off"):
			connectivity_cut_off = float(arg)
		elif opt in ("-w", "--wu"):
			wu = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-c", "--commit"):
			commit = 1
	
	if schema and p_value_cut_off:
		instance = gene_stat_on_mcl_result(dbname, schema, p_value_cut_off, unknown_cut_off, limit, connectivity_cut_off, wu, report, commit)
		instance.dstruc_loadin()
		instance.run()
	else:
		print __doc__
		sys.exit(2)
