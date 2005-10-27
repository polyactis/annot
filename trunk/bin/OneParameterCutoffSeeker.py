#!/usr/bin/env python
"""
Usage: p_gene_lm.py -k SCHEMA -t P_GENE_TABLE -l LM_TABLE [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ...,	the p_gene table
	-l ...,	the lm_table to store the linear_model results, needed if needcommit
	-a ...,	0.5(default)
	-j ...,	how to judge predicted functions, 0(default), 1, 2
	-w ...,	which parameter, see below, (0, default).
	-c, --commit	commit this database transaction
	-r, --report	report flag
	-u, --debug debug flag
	-h, --help              show this help

Examples:
	
Description:
	Calculate the parameter cutoff corresponding to an accuracy_cut_off.
	Which parameter:
"""

import sys, os, getopt
sys.path += [os.path.expanduser('~/script/annot/bin')]
from codense.common import db_connect
from p_gene_lm import p_gene_lm
from heapq import heappush, heappop

def OneParameterCutoffSeeker:
	def __init__(self, hostname=None, dbname=None, schema=None, p_gene_table=None, \
		lm_table=None, accuracy_cut_off=0, judger_type=0, which=0, commit=0, report=0, debug=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.p_gene_table = p_gene_table
		self.lm_table = lm_table
		self.accuracy_cut_off = float(accuracy_cut_off)
		self.judger_type = int(judger_type)
		self.which = int(which)
		self.commit = int(commit)
		self.report = int(report)
		self.debug = int(debug)
		self.which_dict = {0: 'p_value_cut_off',
			1: 'recurrence_cut_off',
			2: 'connectivity_cut_off',
			3: 'cluster_size_cut_off',
			4: 'edge_gradient'}
	
	def get_prediction_heap(self, curs, p_gene_table, is_correct_dict, judger_type, which_dict, which):
		"""
		10-27-05
		"""
		sys.stderr.write("Getting prediction_heap...")
		prediction_heap = []
		curs.execute("DECLARE crs CURSOR FOR select p.%s, p.%s, p.gene_no, p.go_no from %s p"\
			%(which_dict[which], is_correct_dict[judger_type], p_gene_table))
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		while rows:
			for row in rows:
				heappush(prediction_heap, row)
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		sys.stderr.write("Got prediction_heap.\n")
	
	def get_sorted_score_acc_list(self, prediction_heap):
		"""
		10-27-05
		"""
		
	
	def get_cutoff(self, sorted_score_acc_list):
	
	def run(self):
