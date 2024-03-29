#!/usr/bin/env python
"""
Usage: gene_stat_on_mcl_result.py -k SCHEMA -p P_VALUE_CUT_OFF [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	 Ignore it, interface compatibility.
	-m ..., --mcl_table=...	mcl_result(default)
	-g ..., --gene_table=...	table to store the stat results, p_gene(default), needed if commit
	-p ..., --p_value_cut_off=...	p_value_cut_off
	-u ..., --unknown_cut_off=...	unknown_gene_ratio_cut_off, 0.25(default)
	-l ..., --limit=...,	IGNORE it, interface relics.
	-n ..., --connectivity_cut_off=...	0.8(default), minimum connectivity of a mcl cluster
	-y ..., --recurrence_cut_off=...	5(default), minimum recurrences
	-x ..., --cluster_size_cut_off=...	20(default), maximum cluster size
	-w, --wu	IGNORE it. interface compatibility.
	-r, --report	report the progress(a number)
	-c, --commit	commit the database transaction, records in table gene.
	-h, --help              show this help

Examples:
	gene_stat_on_mcl_result.py -k shu -p 0.001
	gene_stat_on_mcl_result.py -k shu -p 0.001 -n 0.7
	gene_stat_on_mcl_result.py -k shu -p 0.001 -n 0.7 -u 0.20 -r -c -m fim_result -g p_gene_fim_result

Description:

"""

import sys, os, psycopg, getopt
from graphlib import Graph
from sets import Set

class gene_prediction:
	# class holding prediction information of a gene
	def __init__(self):
		self.tp = {}
		self.tp1 = {}
		self.tn = 0
		self.fp = {}
		self.fn = {}
		self.p_functions_dict = {}
		self.mcl_id_list = []

class gene_stat_on_mcl_result:
	def __init__(self, hostname, dbname, schema, table, mcl_table, p_value_cut_off, unknown_cut_off, \
		limit, connectivity_cut_off, recurrence_cut_off, cluster_size_cut_off, wu, report=0, needcommit=0, gene_table='p_gene'):
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.table = table
		self.mcl_table = mcl_table
		self.p_value_cut_off = float(p_value_cut_off)
		self.unknown_cut_off = float(unknown_cut_off)
		self.limit = int(limit)
		self.connectivity_cut_off = float(connectivity_cut_off)
		self.recurrence_cut_off = int(recurrence_cut_off)
		self.cluster_size_cut_off = int(cluster_size_cut_off)
		self.wu = int(wu)
		self.report = int(report)
		self.needcommit = int(needcommit)
		self.gene_table = gene_table
		self.tp = 0.0
		self.tp_m = 0.0
		self.tp1 = 0.0
		self.tp1_m = 0.0
		self.tn = 0.0
		self.fp = 0.0
		self.fp_m =0.0
		self.fn = 0.0
		#mapping between gene_no and go_no list
		self.known_genes_dict = {}
		self.no_of_records = 0
		self.log_file = open('/tmp/gene_stat_on_mcl_result.log','w')
		self.gene_prediction_dict = {}
		#self.gono_goindex_dict = {}
		self.no_of_p_known = 0
		#GO DAG
		self.go_graph = Graph.Graph()
		#mapping between go_no and go_id
		self.go_no2go_id = {}
				
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

		sys.stderr.write("Done\n")
		
	def run(self):
		if self.limit != 0:
			sys.stderr.write("gene_stat_on_mcl_result now doesn't support limit\n")
			sys.exit(2)

		self.curs.execute("DECLARE crs CURSOR FOR select mcl_id, vertex_set, p_value_min, go_no_vector, unknown_gene_ratio\
				from %s where connectivity>=%f and p_value_min notnull and array_upper(recurrence_array, 1)>=%d\
				and array_upper(vertex_set, 1)<=%d"%(self.mcl_table, self.connectivity_cut_off, self.recurrence_cut_off, self.cluster_size_cut_off))
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
		self.stat_output()
		
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
		self.no_of_records += 1
		for gene_no in vertex_set:
			if gene_no not in self.gene_prediction_dict:
				item = gene_prediction()
				self.gene_prediction_dict[gene_no] = item
			self.gene_prediction_dict[gene_no].mcl_id_list.append(mcl_id)
			for go_no in go_no_vector:
				if go_no not in self.gene_prediction_dict[gene_no].p_functions_dict:
					self.gene_prediction_dict[gene_no].p_functions_dict[go_no] = 1
				else:
					self.gene_prediction_dict[gene_no].p_functions_dict[go_no] += 1

	def final(self):
		#get the dominant function among all the functions predicted for one gene
		for gene_no in self.gene_prediction_dict:
			entry = self.gene_prediction_dict[gene_no]
			#the biggest support is what the dominant function needs
			dominant_support = max(entry.p_functions_dict.values())
			#delete those functions whose support is less than the dominant support. There's still possibility that >1 functions are kept. 
			for go_no in entry.p_functions_dict.keys():
				if entry.p_functions_dict[go_no] < dominant_support:
					del entry.p_functions_dict[go_no]
		for gene_no in self.gene_prediction_dict:
			#initialize three data structures
			L0_dict = {}
			L1_dict = {}
			#L0 or L1 are all good known functions
			good_k_functions_dict = {}
			#each entry is a gene_prediction class
			entry = self.gene_prediction_dict[gene_no]
			p_functions_dict = entry.p_functions_dict.copy()
			p_functions = p_functions_dict.keys()
			if gene_no not in self.known_genes_dict:
				self.log_file.write('unknown: %d %s %s\n'%(gene_no, repr(p_functions_dict),repr(entry.mcl_id_list)))
				continue
			k_functions_list = self.known_genes_dict[gene_no]
			self.no_of_p_known += 1
			#compare the known go_no with the predicted go_no to see if they match
			#L0 level or L1 level
			for p_go_no in p_functions_dict:
				for k_go_no in k_functions_list:
					if self.is_L0(p_go_no, k_go_no):
						L0_dict[p_go_no] = p_functions_dict[p_go_no]
						good_k_functions_dict[k_go_no] = 1
					elif self.is_L1(p_go_no, k_go_no):
						L1_dict[p_go_no] = p_functions_dict[p_go_no]
						good_k_functions_dict[k_go_no] = 1
			#if one L1 is also counted as L0, remove it from L1_dict.
			#this case happens when one known function has both L0 and L1 matches in the predicted functions
			for go_no in L0_dict:
				if go_no in L1_dict:
					del L1_dict[go_no]
			entry.tp = L0_dict
			entry.tp1 = L1_dict
			for k_go_no in k_functions_list:
				if k_go_no not in good_k_functions_dict:
					entry.fn[k_go_no] = 1
			for p_go_no in p_functions_dict:
				if p_go_no not in L0_dict and p_go_no not in L1_dict:
					entry.fp[p_go_no] = p_functions_dict[p_go_no]
			
			entry.tn = self.no_of_functions - (len(entry.tp)+len(entry.tp1)+len(entry.fp)+len(entry.fn))
			self.tp += len(entry.tp)
			self.tp_m += sum(entry.tp.values())
			self.tp1 += len(entry.tp1)
			self.tp1_m += sum(entry.tp1.values())
			self.tn += entry.tn
			self.fp += len(entry.fp)
			self.fp_m += sum(entry.fp.values())
			self.fn += sum(entry.fn.values())
			self.log_file.write('known: %d %s %s %d %s %s %s %s\n'%(gene_no, repr(entry.tp),repr(entry.tp1),\
				entry.tn,repr(entry.fp),repr(entry.fn),repr(p_functions),repr(entry.mcl_id_list)))
	
	def is_L0(self, p_go_no, k_go_no):
		k_go_id = self.go_no2go_id[k_go_no]
		p_go_id = self.go_no2go_id[p_go_no]
		k_go_family = Set(self.go_graph.forw_bfs(k_go_id))
		if p_go_id in k_go_family:
			self.log_file.write('%d is L0 of %d:: %s %s\n'%(p_go_no, k_go_no, p_go_id, k_go_id))
			return 1
		else:
			return 0
		
	def is_L1(self, p_go_no, k_go_no):
		k_go_id = self.go_no2go_id[k_go_no]
		p_go_id = self.go_no2go_id[p_go_no]
		k_go_inc_nbrs = Set(self.go_graph.inc_nbrs(k_go_id))
		if p_go_id in k_go_inc_nbrs:
			self.log_file.write("%d is direct parent of %d:: %s %s\n"%(p_go_no, k_go_no, p_go_id, k_go_id))	
			return 1
		for k_go_inc_nbr in k_go_inc_nbrs:
			k_go_inc_nbr_out_nbrs = Set(self.go_graph.out_nbrs(k_go_inc_nbr))
			if p_go_id in k_go_inc_nbr_out_nbrs:
				self.log_file.write("%d and %d are siblings:: %s %s\n"%(p_go_no, k_go_no, p_go_id, k_go_id))			
				return 1
		return 0
	
	def list_stringlist(self, list):
		return '{' + repr(list)[1:-1] + '}'
	
	def submit(self):
		sys.stderr.write("Database transacting...")
		if self.gene_table!='p_gene':
			#create the table if it's not 'p_gene'
			self.curs.execute("create table %s(\
				gene_no integer,\
				cluster_array integer[],\
				tp integer[],\
				tp1 integer[],\
				tn integer,\
				fp integer[],\
				fn integer[],\
				p_functions integer[]\
				)"%self.gene_table)
		self.curs.execute("select gene_no from gene")
		rows = self.curs.fetchall()
		for row in rows:
			gene_no = row[0]
			if gene_no in self.gene_prediction_dict:
				entry = self.gene_prediction_dict[gene_no]
				p_functions = entry.p_functions_dict.keys()
				string_tp = self.list_stringlist(entry.tp.keys())
				string_tp1 = self.list_stringlist(entry.tp1.keys())
				string_fp = self.list_stringlist(entry.fp.keys())
				string_fn = self.list_stringlist(entry.fn.keys())
				if gene_no in self.known_genes_dict:
					self.curs.execute("insert into %s(gene_no, cluster_array, tp, tp1, tn, fp, fn, p_functions)\
						values(%d, ARRAY%s, '%s', '%s', %d, '%s', '%s', ARRAY%s)"%\
						(self.gene_table, gene_no, repr(entry.mcl_id_list),string_tp,string_tp1,\
						entry.tn, string_fp, string_fn,repr(p_functions)))
				else:
					self.curs.execute("insert into %s(gene_no, cluster_array, p_functions)\
						values(%d, ARRAY%s, ARRAY%s)"%\
						(self.gene_table, gene_no, repr(entry.mcl_id_list),repr(p_functions)))
		if self.needcommit:				
			self.curs.execute("end")	
		sys.stderr.write("done.\n")

	def stat_output(self):
		sys.stderr.write('\n\tp_value_cut_off:%f unknown_cut_off:%f connectivity_cut_off:%f\n'%(self.p_value_cut_off, self.unknown_cut_off, self.connectivity_cut_off))
		sys.stderr.write('\trecurrence_cut_off:%d cluster_size_cut_off:%d\n'%(self.recurrence_cut_off, self.cluster_size_cut_off))
		sys.stderr.write('\tTotal genes: %d\n'%len(self.gene_prediction_dict))
		sys.stderr.write('\tTotal known genes: %d\n'%self.no_of_p_known)
		sys.stderr.write("\tBased on functions:\n")
		sys.stderr.write('\t\tTP0: %d  TP1: %d  TN: %d  FP: %d  FN: %d\n'%(self.tp, self.tp1, self.tn, self.fp, self.fn))
		if (self.tp+self.tp1+self.fn) == 0:
			sys.stderr.write('\t\tSensitvity: Null\n')
		else:
			sys.stderr.write('\t\tSensitvity: %f\n'%((self.tp+self.tp1)/(self.tp+self.tp1+self.fn)))
		if (self.fp+self.tn) == 0:
			sys.stderr.write('\t\tSpecificity: Null\n')
		else:
			sys.stderr.write('\t\tSpecificity: %f\n'%(self.tn/(self.fp+self.tn)))
		if (self.tp+self.tp1+self.fp) == 0:
			sys.stderr.write('\t\tFalse Positive Ratio: Null\n')
		else:
			sys.stderr.write('\t\tFalse Positive Ratio: %f\n'%(self.fp/(self.tp+self.tp1+self.fp)))
		sys.stderr.write("\tBased on clusters:\n")
		sys.stderr.write('\t\tTP0_M: %d  TP1_M: %d  FP_M: %d\n'%(self.tp_m, self.tp1_m, self.fp_m))
		if (self.tp_m+self.tp1_m+self.fp_m) == 0:
			sys.stderr.write('\t\tFalse Positive Ratio: Null\n')
		else:
			sys.stderr.write('\t\tFalse Positive Ratio: %f\n'%(self.fp_m/(self.tp_m+self.tp1_m+self.fp_m)))

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema=", "table=", "mcl_table=", "p_value_cut_off=",\
		"unknown_cut_off=", "limit=", "connectivity_cut_off=", "recurrence_cut_off=", "cluster_size_cut_off=", "wu", "report", "commit", "gene_table="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:m:p:u:l:n:y:x:wrcg:", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = 'mcl_result'
	mcl_table = 'mcl_result'
	p_value_cut_off = None
	unknown_cut_off = 0.25
	limit = 0
	connectivity_cut_off = 0.8
	recurrence_cut_off = 5
	cluster_size_cut_off = 20
	wu = 0
	report = 0
	commit = 0
	gene_table = 'p_gene'
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
		elif opt in ("-m", "--mcl_table"):
			mcl_table = arg
		elif opt in ("-p", "--p_value_cut_off"):
			p_value_cut_off = float(arg)
		elif opt in ("-u", "--unknown_cut_off"):
			unknown_cut_off = float(arg)
		elif opt in ("-l", "--limit"):
			limit = int(arg)
		elif opt in ("-n", "--connectivity_cut_off"):
			connectivity_cut_off = float(arg)
		elif opt in ("-y", "--recurrence_cut_off"):
			recurrence_cut_off = int(arg)
		elif opt in ("-x", "--cluster_size_cut_off"):
			cluster_size_cut_off = int(arg)
		elif opt in ("-w", "--wu"):
			wu = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-g", "--gene_table"):
			gene_table = arg
	
	if schema and p_value_cut_off:
		instance = gene_stat_on_mcl_result(hostname, dbname, schema, table, mcl_table, p_value_cut_off,\
			unknown_cut_off, limit, connectivity_cut_off, recurrence_cut_off, cluster_size_cut_off, wu, report, commit, gene_table)
		instance.dstruc_loadin()
		instance.run()
	else:
		print __doc__
		sys.exit(2)
