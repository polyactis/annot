#!/usr/bin/env python
"""
Usage: gene_stat.py -k SCHEMA -p P_VALUE_CUT_OFF [OPTION]

Option:
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-p ..., --p_value_cut_off=...	p_value_cut_off
	-l ..., --limit=...,	maximum number of clusters related to one gene, 10(default)
	-n ..., --connectivity_cut_off=...	0.8(default), minimum connectivity of a mcl cluster
	-w, --wu	Wu's strategy(Default is Jasmine's strategy)
	-r, --report	report the progress(a number)
	-c, --commit	commit the database transaction, records in table gene.
	-h, --help              show this help

Examples:
	gene_stat.py -k shu -p 0.001
	gene_stat.py -k shu -l 0 -p 0.001 -n 0.7

Description:
	This program is mainly for validation purpose. Run after cluster_stat.py.
	Implements both Jasmine and Wu's strategy.
"""

import sys, os, psycopg, getopt
from rpy import *


class gene_stat:
	def __init__(self, dbname, schema, p_value_cut_off, limit=10, connectivity_cut_off=0.8, wu=0, report=0, needcommit=0):
		self.conn = psycopg.connect('dbname=%s'%dbname)
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.p_value_cut_off = float(p_value_cut_off)
		self.limit = int(limit)
		self.connectivity_cut_off = float(connectivity_cut_off)
		self.wu = int(wu)
		self.report = int(report)
		self.needcommit = int(needcommit)
		self.tp = 0.0
		self.tn = 0.0
		self.fp = 0.0
		self.fn = 0.0
		self.known_genes_dict = {}
		self.no_of_records = 0
		self.gene_to_be_updated = []
		self.log_file = open('/tmp/gene_stat.log','w')
		
	def dstruc_loadin(self):
		self.curs.execute("select gene_no,go_functions from gene where known=TRUE")
		rows = self.curs.fetchall()
		for row in rows:
			go_functions_list = row[1][1:-1].split(',')
			self.known_genes_dict[row[0]] = []
			for go_no in go_functions_list:
				self.known_genes_dict[row[0]].append(int(go_no))
		if self.wu:
			self.curs.execute("select count(go_no) from go")
		else:
			self.curs.execute("select count(go_no) from go where go_no!=0")
		rows = self.curs.fetchall()
		self.no_of_functions = rows[0][0]

	def run(self):
		for gene_no in self.known_genes_dict:
			mcl_id_list = []
			self.log_file.write('%d '%gene_no)
			if self.limit == 0:
				self.curs.execute("select mcl_id, p_value_vector \
					from cluster_stat where leave_one_out=%d \
					and connectivity>=%f"%(gene_no, self.connectivity_cut_off))
			else:
				self.curs.execute("select mcl_id, p_value_vector \
					from cluster_stat where leave_one_out=%d \
					and connectivity>=%f limit %d"%(gene_no, self.connectivity_cut_off, self.limit))
			rows = self.curs.fetchall()
			p_functions_dict = {}
			for row in rows:
				mcl_id_list.append(row[0])
				p_value_vector = row[1][1:-1].split(',')
				for i in range(self.no_of_functions):
					if float(p_value_vector[i]) <= self.p_value_cut_off:
						if self.wu:
							p_functions_dict[i] = 1
						else:
							p_functions_dict[i+1] = 1
			if rows:
				self._gene_stat(gene_no, mcl_id_list, p_functions_dict)
				self.no_of_records += 1
				if self.report:
					sys.stderr.write('%s%s'%('\x08'*20, self.no_of_records))
		if self.needcommit:
			for entry in self.gene_to_be_updated:
				if entry[6]:
					self.curs.execute("update gene set cluster_array=ARRAY%s,\
					tp=%d,tn=%d,fp=%d,fn=%d,p_functions=ARRAY%s where gene_no=%d"%\
					(repr(entry[1]),entry[2],entry[3],entry[4],entry[5],repr(entry[6]),entry[0]))
				else:
					self.curs.execute("update gene set cluster_array=ARRAY%s,\
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
		tn = self.no_of_functions - (tp+fp+fn)
		self.tp += tp
		self.tn += tn
		self.fp += fp
		self.fn += fn
		self.log_file.write('%d %d %d %d %s %s\n'%(tp,tn,fp,fn,repr(p_functions),repr(mcl_id_list)))
		self.gene_to_be_updated.append([gene_no,mcl_id_list,tp,tn,fp,fn,p_functions])

class gene_prediction:
	def __init__(self):
		self.tp = 0
		self.tn = 0
		self.fp = 0
		self.fn = 0
		self.p_functions_dict = {}
		self.mcl_id_list = []


class gene_stat_new:
	def __init__(self, dbname, schema, p_value_cut_off, limit=0, connectivity_cut_off=0.8, wu=0, report=0, needcommit=0):
		self.conn = psycopg.connect('dbname=%s'%dbname)
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.p_value_cut_off = float(p_value_cut_off)
		self.limit = int(limit)
		self.connectivity_cut_off = float(connectivity_cut_off)
		self.wu = int(wu)
		self.report = int(report)
		self.needcommit = int(needcommit)
		self.tp = 0.0
		self.tn = 0.0
		self.fp = 0.0
		self.fn = 0.0
		self.known_genes_dict = {}
		self.no_of_records = 0
		self.log_file = open('/tmp/gene_stat.log','w')
		self.gene_prediction_dict = {}
		
	def dstruc_loadin(self):
		self.curs.execute("select gene_no,go_functions from gene where known=TRUE")
		rows = self.curs.fetchall()
		for row in rows:
			go_functions_list = row[1][1:-1].split(',')
			self.known_genes_dict[row[0]] = []
			for go_no in go_functions_list:
				self.known_genes_dict[row[0]].append(int(go_no))
		if self.wu:
			self.curs.execute("select count(go_no) from go")
		else:
			self.curs.execute("select count(go_no) from go where go_no!=0")
		rows = self.curs.fetchall()
		self.no_of_functions = rows[0][0]

	def run(self):
		if self.limit != 0:
			sys.stderr.write("gene_stat_new now doesn't support limit\n")
			sys.exit(2)

		self.curs.execute("DECLARE crs CURSOR FOR select mcl_id, leave_one_out, p_value_vector \
				from cluster_stat where connectivity>=%f"%(self.connectivity_cut_off))
		self.curs.execute("fetch 5000 from crs")
		rows = self.curs.fetchall()
		while rows:
			for row in rows:
				if self.known_genes_dict.has_key(row[1]):
				# for clusters containing unknown genes
					self._gene_stat(row)
					self.no_of_records += 1
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, self.no_of_records))
			
			self.curs.execute("fetch 5000 from crs")
			rows = self.curs.fetchall()
		
		self.final()
		if self.needcommit:
			for gene_no in self.gene_prediction_dict:
				entry = self.gene_prediction_dict[gene_no]
				p_functions = entry.p_functions_dict.keys()
				if p_functions:
					self.curs.execute("update gene set cluster_array=ARRAY%s,\
					tp=%d,tn=%d,fp=%d,fn=%d,p_functions=ARRAY%s where gene_no=%d"%\
					(repr(entry.mcl_id_list),entry.tp,entry.tn,entry.fp,entry.fn,repr(p_functions),gene_no))
				else:
					self.curs.execute("update gene set cluster_array=ARRAY%s,\
					tp=%d,tn=%d,fp=%d,fn=%d where gene_no=%d"%\
					(repr(entry.mcl_id_list),entry.tp,entry.tn,entry.fp,entry.fn,gene_no))
			self.curs.execute("end")
		sys.stderr.write('\n\tTotal genes: %d\n'%self.no_of_records)
		sys.stderr.write('\tSensitvity: %f\n'%(self.tp/(self.tp+self.fn)))
		sys.stderr.write('\tSpecificity: %f\n'%(self.tn/(self.fp+self.tn)))
		sys.stderr.write('\tFalse Positive Ratio: %f\n'%(self.fp/(self.tp+self.fp)))
			
	def _gene_stat(self, row):
		gene_no = row[1]
		if gene_no not in self.gene_prediction_dict:
			item = gene_prediction()
			self.gene_prediction_dict[gene_no] = item
		self.gene_prediction_dict[gene_no].mcl_id_list.append(row[0])
		p_value_vector = row[2][1:-1].split(',')
		for i in range(self.no_of_functions):
			if float(p_value_vector[i]) <= self.p_value_cut_off:
				if self.wu:
					self.gene_prediction_dict[gene_no].p_functions_dict[i] = 1
				else:
					self.gene_prediction_dict[gene_no].p_functions_dict[i+1] = 1
		
	def final(self):
		for gene_no in self.gene_prediction_dict:
			p_functions_dict = self.gene_prediction_dict[gene_no].p_functions_dict.copy()
			p_functions = p_functions_dict.keys()
			entry = self.gene_prediction_dict[gene_no]
			for go_no in self.known_genes_dict[gene_no]:
				if go_no in p_functions_dict:
					entry.tp += 1
					del p_functions_dict[go_no]
				else:
					entry.fn += 1
			entry.fp = len(p_functions_dict)
			entry.tn = self.no_of_functions - (entry.tp+entry.fp+entry.fn)

			self.tp += entry.tp
			self.tn += entry.tn
			self.fp += entry.fp
			self.fn += entry.fn
			self.log_file.write('%d %d %d %d %d %s %s\n'%(gene_no, entry.tp,entry.tn,entry.fp,entry.fn,repr(p_functions),repr(entry.mcl_id_list)))


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "dbname=", "schema=", "p_value_cut_off=", "limit=", "connectivity_cut_off=", "wu", "report", "commit"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hd:k:p:l:n:wrc", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	dbname = 'graphdb'
	schema = ''
	p_value_cut_off = None
	limit = 10
	connectivity_cut_off = 0.8
	wu = 0
	report = 0
	commit = 0
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
		instance = gene_stat_new(dbname, schema, p_value_cut_off, limit, connectivity_cut_off, wu, report, commit)
		instance.dstruc_loadin()
		instance.run()
	else:
		print __doc__
		sys.exit(2)
