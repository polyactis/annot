#!/usr/bin/env python
"""
Usage: gene_stat.py -k SCHEMA -p P_VALUE_CUT_OFF [OPTION]

Option:
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-p ..., --p_value_cut_off=...	p_value_cut_off
	-u ..., --unknown_cut_off=...	unknown_cut_off, 0.85(default), for wu's
	-l ..., --limit=...,	OBSOLETE, accepted for backwards compatibility
	-n ..., --connectivity_cut_off=...	0.8(default), minimum connectivity of a mcl cluster
	-w, --wu	Wu's strategy(Default is Jasmine's strategy)
	-r, --report	report the progress(a number)
	-c, --commit	commit the database transaction, records in table gene.
	-h, --help              show this help

Examples:
	gene_stat.py -k shu -p 0.001 -w
	gene_stat.py -k shu -p 0.001 -n 0.7 -w -r
	gene_stat.py -k shu -p 0.001 -n 0.7 -u 0.80 -w

Description:
	This program is mainly for validation purpose. Run after cluster_stat.py.
	Implements both Jasmine and Wu's strategy.
	The unknown_cut_off is different from the one in gene_stat_on_mcl_result.
	The former is unknown class p_value cut_off.
	The latter is unknown genes ratio.
"""

import sys, os, psycopg, getopt
from rpy import *
from graphlib import Graph
from sets import Set
from gene_stat_on_mcl_result import gene_stat_on_mcl_result, gene_prediction


class gene_stat(gene_stat_on_mcl_result):
	def __init__(self, dbname, schema, p_value_cut_off, unknown_cut_off, limit=0, \
		connectivity_cut_off=0.8, wu=0, report=0, needcommit=0):
		gene_stat_on_mcl_result.__init__(self, dbname, schema, p_value_cut_off, \
			unknown_cut_off, limit, connectivity_cut_off, wu, report, needcommit)

	def dstruc_loadin(self):
		sys.stderr.write("Loading Data STructure...")
		#setup self.known_genes_dict
		self.curs.execute("select gene_no,go_functions from gene where known=TRUE")
		rows = self.curs.fetchall()
		for row in rows:
			go_functions_list = row[1][1:-1].split(',')
			self.known_genes_dict[row[0]] = []
			for go_no in go_functions_list:
				self.known_genes_dict[row[0]].append(int(go_no))
		
		#get the non-obsolete biological_process GO DAG		
		self.curs.execute("select count(go_no) from go")
		rows = self.curs.fetchall()
		self.no_of_functions = rows[0][0]
		self.curs.execute("select t2t.term1_id, t2t.term2_id, t1.acc, t2.acc from \
			go.term2term t2t, go.term t1, go.term t2 where t2t.term1_id=t1.id and \
			t2t.term2_id=t2.id and t1.is_obsolete=0 and t2.is_obsolete=0 and \
			t1.term_type='biological_process' and t2.term_type='biological_process' ")
		rows = self.curs.fetchall()
		for row in rows:
		#setup the go_graph structure
			self.go_graph.add_edge(row[2], row[3])
		
		#setup self.go_no2go_id
		self.curs.execute("select go_no, go_id from go")
		rows = self.curs.fetchall()
		for row in rows:
			self.go_no2go_id[row[0]] = row[1]
		
		#setup self.no_of_functions
		if self.wu:
			self.curs.execute("select count(go_no) from go")
		else:
			self.curs.execute("select count(go_no) from go where go_no!=0")
		rows = self.curs.fetchall()
		self.no_of_functions = rows[0][0]
		sys.stderr.write("Done\n")

	def run(self):
		if self.limit != 0:
			sys.stderr.write("gene_stat now doesn't support limit\n")
			sys.exit(2)

		self.curs.execute("DECLARE crs CURSOR FOR select mcl_id, leave_one_out, p_value_vector \
				from cluster_stat where connectivity>=%f"%(self.connectivity_cut_off))
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
		#Database updating is too slow. Do it only if needcommit.
			self.submit()
		sys.stderr.write('\n\tp_value_cut_off:%f unknown_cut_off:%f connectivity_cut_off:%f\n'%(self.p_value_cut_off, self.unknown_cut_off, self.connectivity_cut_off))
		sys.stderr.write('\tTotal genes: %d\n'%len(self.gene_prediction_dict))
		sys.stderr.write('\tTotal known genes: %d\n'%self.no_of_p_known)
		sys.stderr.write('\tWhole Sensitvity: %f\n'%((self.tp+self.tp1)/(self.tp+self.tp1+self.fn)))
		sys.stderr.write('\tTP0: %d  TP1: %d  TN: %d  FP: %d  FN: %d\n'%(self.tp, self.tp1, self.tn, self.fp, self.fn))
		sys.stderr.write('\tTP0_M: %d  TP1_M: %d  FP_M: %d\n'%(self.tp_m, self.tp1_m, self.fp_m))
		sys.stderr.write('\tSpecificity: %f\n'%(self.tn/(self.fp+self.tn)))
		sys.stderr.write('\tFalse Positive Ratio: %f\n'%(self.fp_m/(self.tp_m+self.tp1_m+self.fp_m)))
			
	def _gene_stat(self, row):
		p_value_vector = row[2][1:-1].split(',')
		#transform into float type
		p_value_vector = map(float, p_value_vector)
		min_p_value =min(p_value_vector)
		if self.wu:
			if float(p_value_vector[0]) <= self.unknown_cut_off:
			#too many unknown genes
				return
		if min_p_value > self.p_value_cut_off:
		#none of the predicted functions in this cluster is significant
			return
		self.no_of_records += 1
		mcl_id = row[0]
		gene_no = row[1]
		if gene_no not in self.gene_prediction_dict:
			item = gene_prediction()
			self.gene_prediction_dict[gene_no] = item
		self.gene_prediction_dict[gene_no].mcl_id_list.append(mcl_id)				
		for i in range(self.no_of_functions):
			if p_value_vector[i] == min_p_value:
				if self.wu:
				#index 0 corresponds to go_no 0.
					go_no = i
				else:
				#index 0 corresponds to go_no 1
					go_no = i+1
				if go_no not in self.gene_prediction_dict[gene_no].p_functions_dict:
				#value in p_functions_dict stores the number of associated clusters.
					self.gene_prediction_dict[gene_no].p_functions_dict[go_no] = 1
				else:
					self.gene_prediction_dict[gene_no].p_functions_dict[go_no] += 1

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
	unknown_cut_off = 0.85
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
		instance = gene_stat(dbname, schema, p_value_cut_off, unknown_cut_off, limit, connectivity_cut_off, wu, report, commit)
		instance.dstruc_loadin()
		instance.run()
	else:
		print __doc__
		sys.exit(2)
