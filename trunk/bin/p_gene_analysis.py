#!/usr/bin/env python
"""
Usage: p_gene_analysis.py -k SCHEMA -p P_VALUE_CUT_OFF [OPTION] STAT_TABLE_FILE

Option:
	STAT_TABLE_FILE is the file where the output goes.
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	cluster_stat(default)
	-m ..., --mcl_table=...	mcl_result(default), mcl_result table corresponding to above table.
	-g ..., --gene_table=...	table to store the stat results, p_gene(default), needed if commit
	-p ..., --p_value_cut_off=...	p_value_cut_off
	-j ..., --judger_type=...	how to judge predicted functions, 0(default), 1, 2
	-r, --report	report the progress(a number)
	-c, --commit	commit the database transaction, records down the go_no2accuracy.
	-a ..., --accuracy_cut_off=...	the accuracy_cut_off to be based for p_value adjusting, 0(default)
		NOTICE: 0 means using p_value_cut_off method
	-b, --debug	enable debugging, no debug by default
	-h, --help              show this help

Examples:
	p_gene_analysis.py -k sc_54 -g p_gene_e5 -p 0.001 stat_table
	p_gene_analysis.py -k sc_54 -p 0.01 -c -j 2  -g p_gene_cluster_stat2 stat_table

Description:
	02-28-05
	this is the second part of gene_stat_plot.py, which works on p_gene table
	the first part goes to gene_stat.py
"""

import sys, os, psycopg, getopt, csv, math
from sets import Set
from rpy import r
from numarray import *

class prediction_space_attr:
	"""
	data structure for prediction_space2attr
	"""
	def __init__(self):
		self.correct_predictions = 0.0
		self.known_predictions = 0.0
		self.unknown_predictions = 0.0
		self.correct_predictions_pair = 0.0
		self.known_predictions_pair = 0.0
		self.unknown_predictions_pair = 0.0
		self.prediction_pair_dict = {}
		self.tuple_list = []	#add for p_gene_lm.py
		self.p_value_cut_off = None	#add for p_gene_lm.py

class prediction_pair_attr:
	"""
	data structure for prediction_pair2attr
	"""
	def __init__(self):
		self.is_correct = -1
		self.supporting_clusters = 0
		self.go_no_mapped = -1

class accuracy_struc:
	#data structure for self.go_no2accuracy
	def __init__(self):
		self.unknown = 0.0
		self.known = 0.0
		self.correct = 0.0
		self.ratio = 0

class p_gene_analysis:
	"""
	dstruc_loadin()
	run()
		--data_fetch()
			--_p_gene_analysis()
		--overview_stats()
			--return_known_unknown_gene_sets()
		--go_no_accuracy()
		--table_output()
	"""
	def __init__(self, hostname, dbname, schema, table, mcl_table, p_value_cut_off, report=0, \
		judger_type=0, needcommit=0, gene_table='p_gene', \
		stat_table_fname=None, debug=0, accuracy_cut_off=0):
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.schema = schema
		self.curs.execute("set search_path to %s"%schema)
		self.table = table
		self.mcl_table = mcl_table
		self.p_value_cut_off = float(p_value_cut_off)
		self.report = int(report)
		self.judger_type = int(judger_type)
		self.needcommit = int(needcommit)
		self.gene_table = gene_table
		self.stat_table_fname = stat_table_fname
		#debugging flag
		self.debug = int(debug)
		self.accuracy_cut_off = float(accuracy_cut_off)
		
		#open a file
		if self.stat_table_fname:
			self.stat_table_f = csv.writer(open(self.stat_table_fname, 'w'), delimiter='\t')

		#an is_correct dictionary used in database fetch
		self.is_correct_dict = {0: 'is_correct',
			1: 'is_correct_L1',
			2: 'is_correct_lca'}
		
		#a counter
		self.no_of_records = 0
		
		#data structures to be loaded in
		self.known_genes_dict = {}
		
		#data structures to be created
		self.prediction_space2attr = {}
		self.prediction_pair2attr = {}
		self.go_no2accuracy = {}
		
	def dstruc_loadin(self):
		sys.stderr.write("Loading Data STructure...")
		
		#setup self.known_genes_dict
		self.curs.execute("select gene_no,go_functions from gene where known=TRUE")
		rows = self.curs.fetchall()
		for row in rows:
			go_functions_list = row[1][1:-1].split(',')
			self.known_genes_dict[row[0]] = Set()
			for go_no in go_functions_list:
				self.known_genes_dict[row[0]].add(int(go_no))
		
		sys.stderr.write("Done\n")

	def run(self):
		self.data_fetch()
		self.overview_stats()
		self.go_no_accuracy()
		self.table_output()

	def data_fetch(self):
		"""
		02-21-05
			--_p_gene_analysis()
		"""
		self.curs.execute("DECLARE crs CURSOR FOR select gene_no, go_no, mcl_id, %s, avg_p_value, \
			recurrence_cut_off,connectivity_cut_off, depth_cut_off from %s"\
			%(self.is_correct_dict[self.judger_type], self.gene_table))
		
		self.curs.execute("fetch 5000 from crs")
		rows = self.curs.fetchall()
		while rows:
			for row in rows:
				self._p_gene_analysis(row)
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, self.no_of_records))
			
			self.curs.execute("fetch 5000 from crs")
			rows = self.curs.fetchall()

	def _p_gene_analysis(self, row):
		gene_no = row[0]
		go_no = row[1]
		mcl_id = row[2]
		is_correct = row[3]
		p_value = row[4]
		recurrence = row[5]
		connectivity = row[6]
		depth_cut_off = row[7]
		
		if p_value > self.p_value_cut_off:
			return
		
		prediction_space = (recurrence, connectivity)
		prediction_pair = (gene_no, go_no)
		"""
		1 deal with self.prediction_space2attr
		"""
		if prediction_space not in self.prediction_space2attr:
			#initalize its attribute if it doesn't exist
			self.prediction_space2attr[prediction_space] = prediction_space_attr()
		#pass the value to a value, ease programming
		unit = self.prediction_space2attr[prediction_space]
		
		unit.correct_predictions += (is_correct+1)/2	#increment when is_correct=1
		unit.known_predictions += int(is_correct>=0)	#increment when is_correct=0 or 1
		unit.unknown_predictions += -(is_correct-1)/2	#increment when is_correct = -1
		if prediction_pair not in unit.prediction_pair_dict:
			#pair only counts the first time
			unit.prediction_pair_dict[prediction_pair] = [p_value]
			unit.correct_predictions_pair += (is_correct+1)/2
			unit.known_predictions_pair += int(is_correct>=0)
			unit.unknown_predictions_pair += -(is_correct-1)/2
		else:
			unit.prediction_pair_dict[prediction_pair].append(p_value)
		"""
		2. deal with self.prediction_pair2attr
		"""
		if prediction_pair not in self.prediction_pair2attr:
			#initalize its attribute if first time
			self.prediction_pair2attr[prediction_pair] = prediction_pair_attr()
		
		unit_pair = self.prediction_pair2attr[prediction_pair]
		unit_pair.is_correct = is_correct
		unit_pair.supporting_clusters += 1

	def overview_stats(self):
		'''
		02-21-05
			get some global statistics
			--return_known_unknown_gene_sets()
		'''
		total_correct_predictions = 0.0	#0.0 means float, 0 means int. int/int is still an int.
		total_known_predictions = 0.0
		total_unknown_predictions = 0.0
		for (prediction_space,unit) in self.prediction_space2attr.iteritems():
			total_correct_predictions += unit.correct_predictions
			total_known_predictions += unit.known_predictions
			total_unknown_predictions += unit.unknown_predictions
		
		total_correct_predictions_pair = 0.0
		total_known_predictions_pair = 0.0
		total_unknown_predictions_pair = 0.0
		for (prediction_pair,unit) in self.prediction_pair2attr.iteritems():
			is_correct = unit.is_correct
			total_correct_predictions_pair += (is_correct+1)/2	#increment when is_correct=1
			total_known_predictions_pair += int(is_correct>=0)	#increment when is_correct=0 or 1
			total_unknown_predictions_pair += -(is_correct-1)/2	#increment when is_correct = -1
		
		if total_known_predictions != 0:
			total_accuracy = total_correct_predictions/total_known_predictions
		else:
			total_accuracy = -1
		if total_known_predictions_pair !=0:
			total_accuracy_pair = total_correct_predictions_pair/total_known_predictions_pair
		else:
			total_accuracy_pair = -1
		
		(set1, set2) = self.return_known_unknown_gene_sets(self.prediction_pair2attr.keys())
		no_of_known_genes = len(set1)
		no_of_unknown_genes = len(set2)
		
		#output
		self.stat_table_f.writerow(['accuracy', total_accuracy])
		self.stat_table_f.writerow(['known_predictions', total_known_predictions])
		self.stat_table_f.writerow(['unknown_predictions', total_unknown_predictions])
		self.stat_table_f.writerow(['accuracy_pair', total_accuracy_pair])
		self.stat_table_f.writerow(['known_predictions-pair', total_known_predictions_pair])
		self.stat_table_f.writerow(['unknown_predictions-pair', total_unknown_predictions_pair])
		
		self.stat_table_f.writerow(['known genes', no_of_known_genes])
		self.stat_table_f.writerow(['unknown genes', no_of_unknown_genes])
	
	def return_known_unknown_gene_sets(self, prediction_pair_list):
		known_gene_set = Set()
		unknown_gene_set = Set()
		for prediction_pair in prediction_pair_list:
			gene_no = prediction_pair[0]
			if gene_no in self.known_genes_dict:
				known_gene_set.add(gene_no)
			else:
				unknown_gene_set.add(gene_no)
		return (known_gene_set, unknown_gene_set)
	

	def go_no_accuracy(self):
		"""
		this function computes the prediction accuracy for each go function
		"""
		for prediction_pair in self.prediction_pair2attr:
			go_no = prediction_pair[1]
			is_correct = self.prediction_pair2attr[prediction_pair].is_correct
			if go_no not in self.go_no2accuracy:
				self.go_no2accuracy[go_no] = accuracy_struc()
			self.go_no2accuracy[go_no].correct += (is_correct+1)/2	#increment when is_correct=1
			self.go_no2accuracy[go_no].known += int(is_correct>=0)	#increment when is_correct=0 or 1
			self.go_no2accuracy[go_no].unknown += -(is_correct-1)/2	#increment when is_correct = -1
		
		if self.needcommit:
			#this table is easy for linking.
			go_acc_table = '%s_acc_%s'%(self.gene_table, self.is_correct_dict[self.judger_type])
			self.curs.execute("create table %s(\
				go_no	integer,\
				accuracy	float,\
				known	integer,\
				unknown	integer\
				)"%go_acc_table)
		
		#finally, compute all the ratios and output
		go_no_list = self.go_no2accuracy.keys()
		go_no_list.sort()
		self.stat_table_f.writerow([])
		self.stat_table_f.writerow(['go_no', 'accuracy', 'known', 'unknown'])
		for go_no in go_no_list:
			unit = self.go_no2accuracy[go_no]
			if unit.known !=0:
				unit.ratio = unit.correct/unit.known
			else:
				unit.ratio = -1
			self.stat_table_f.writerow([go_no, unit.ratio, unit.known, unit.unknown])
			if self.needcommit:
				self.curs.execute("insert into %s(go_no, accuracy, known, unknown)\
					values(%s, %f, %d, %d)"%(go_acc_table,\
					go_no, unit.ratio, unit.known, unit.unknown))
		
		if self.needcommit:
			self.curs.execute("end")

	def table_output(self):
		'''
		
		'''

		header = ['recurrence', 'connectivity', 'accuracy', 'known_predictions', 'unknown_predictions',\
		'accuracy_pair', 'known_predictions_pair', 'unknown_predictions_pair', 'known genes', 'unknown genes']
		self.stat_table_f.writerow(header)
		prediction_space_list = self.prediction_space2attr.keys()
		prediction_space_list.sort()
		for prediction_space in prediction_space_list:
			unit = self.prediction_space2attr[prediction_space]
			recurrence = prediction_space[0]
			connectivity = prediction_space[1]
			if unit.known_predictions != 0:
				accuracy = unit.correct_predictions/unit.known_predictions
			else:
				accuracy = -1
			if unit.known_predictions_pair != 0:
				accuracy_pair = unit.correct_predictions_pair/unit.known_predictions_pair
			else:
				accuracy_pair = -1
			(set1, set2) = self.return_known_unknown_gene_sets(unit.prediction_pair_dict.keys())
			no_of_known_genes = len(set1)
			no_of_unknown_genes = len(set2)
			row = [recurrence, connectivity, accuracy, unit.known_predictions, unit.unknown_predictions,\
				accuracy_pair, unit.known_predictions_pair, unit.unknown_predictions_pair, no_of_known_genes,\
				no_of_unknown_genes]
			self.stat_table_f.writerow(row)
		


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema=", "table=", "mcl_table=", "p_value_cut_off=",\
		"judger_type=", "report", "commit", "gene_table=", "debug", "accuracy_cut_off="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:m:p:j:rcg:ba:", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = 'cluster_stat'
	mcl_table = 'mcl_result'
	p_value_cut_off = None
	judger_type = 0
	report = 0
	commit = 0
	gene_table = 'p_gene'
	debug = 0
	accuracy_cut_off = 0
	stat_table_fname = None
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
		elif opt in ("-j", "--judger_type"):
			judger_type = int(arg)
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-g", "--gene_table"):
			gene_table = arg
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-a", "--accuracy_cut_off="):
			accuracy_cut_off = float(arg)
	if len(args) == 1:
		stat_table_fname = args[0]
			
	if schema and p_value_cut_off and stat_table_fname:
		instance = p_gene_analysis(hostname, dbname, schema, table, mcl_table, p_value_cut_off,\
			report, judger_type, commit, gene_table, stat_table_fname, debug, accuracy_cut_off)
		instance.dstruc_loadin()
		instance.run()
	else:
		print __doc__
		sys.exit(2)
