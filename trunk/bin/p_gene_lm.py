#!/usr/bin/env python
"""
Usage: p_gene_lm.py -k SCHEMA -t TABLE [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	the p_gene table
	-l ,,,, --lm_table=...	the lm_table to store the linear_model results, needed if needcommit
	-a ..., --accuracy_cut_off=...	0.5(default)
	-j ..., --judger_type=...	how to judge predicted functions, 0(default), 1, 2
	-m ..., --min_data_points=...	the minimum data points to do linear-model fitting, 5(default)
	-v ..., --valid_space=...	the min number of known_predictions in one prediction space, 20(default)
	-c, --commit	commit this database transaction
	-r, --report	report flag
	-u, --debug debug flag
	-h, --help              show this help

Examples:
	p_gene_lm.py -k sc_54 -t p_gene_repos_2_e5 -j 1 -r
	p_gene_lm.py -k sc_54 -t p_gene_repos_2_e5 -l p_gene_repos_2_e5_lm -j 1 -r -c
	
Description:
	02-28-05
	linear model fitting of the p_value_cut_off ~ recurrence+connectivity
	
	Output contains the p_value_cut_off and other information for each go-no
	as well as the linear_model fitting results.
	
	Output is dumped on sys.stdout.
"""

import sys, os, psycopg, getopt, csv
from p_gene_analysis import prediction_space_attr
from module_cc.linear_model import linear_model
from codense.common import *
from numarray import *

class p_gene_lm:
	"""
	02-28-05
		get the data from p_gene table
		for each go, for each (recurrence, connectivity), compute its p_value_cut_off based on accuracy_cut_off
		for each go, linear-model fitting of p_value_cut_off ~ recurrence + connectivity
		submit the coefficients and chisq of linear-model to a lm_table.
		
		parameter min_data_points (no of known predictions) controls the validity of each
		prediction space(recurrence, connectivity).
		
		Some go-nos don't have enough data to do linear-model fitting, simply ignore. But lm_table
		will contain an average model for these go-nos to look up in the later prediction stage.
	
	"""
	def __init__(self, hostname=None, dbname=None, schema=None, table=None, \
		lm_table=None, accuracy_cut_off=0, judger_type=0, min_data_points=5, \
		needcommit=0, report=0, debug=0, valid_space=20):
		
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.table = table
		self.lm_table = lm_table
		self.accuracy_cut_off = float(accuracy_cut_off)
		self.judger_type = int(judger_type)
		self.min_data_points = int(min_data_points)
		self.needcommit = int(needcommit)
		self.report = int(report)		
		self.debug = int(debug)
		self.valid_space = int(valid_space)

		self.go_no2prediction_space = {}
		#an is_correct dictionary used in database fetch
		self.is_correct_dict = {0: 'is_correct',
			1: 'is_correct_L1',
			2: 'is_correct_lca'}
		
		self.go_no2lm_results = {}

		
	def init(self):
		self.lm_instance = linear_model()
		
	def data_fetch(self, curs, table):
		"""
		02-28-05
			borrowed from p_gene_analysis.py
			
			--prediction_space_setup()
		"""
		curs.execute("DECLARE crs CURSOR FOR select gene_no, go_no, mcl_id, %s, avg_p_value, \
			recurrence_cut_off,connectivity_cut_off, depth_cut_off from %s"\
			%(self.is_correct_dict[self.judger_type], table))
		no_of_records = 0
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		while rows:
			for row in rows:
				self.prediction_space_setup(row)
				no_of_records += 1
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, no_of_records))
			
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()

	def prediction_space_setup(self, row):
		"""
		02-28-05
			borrowed from _p_gene_analysis() of p_gene_analysis.py
		"""
		gene_no = row[0]
		go_no = row[1]
		mcl_id = row[2]
		is_correct = row[3]
		p_value = row[4]
		recurrence = row[5]
		connectivity = row[6]
		depth_cut_off = row[7]
		
		if go_no not in self.go_no2prediction_space:
			self.go_no2prediction_space[go_no] = {}
		#pass the value to ease programming
		prediction_space2attr = self.go_no2prediction_space[go_no]
		prediction_space = (recurrence, connectivity)
		"""
		1 deal with self.prediction_space2attr
		"""
		
		if prediction_space not in prediction_space2attr:
			#initalize its attribute if it doesn't exist
			prediction_space2attr[prediction_space] = prediction_space_attr()
		#pass the value to a value, ease programming
		unit = prediction_space2attr[prediction_space]
		
		unit.correct_predictions += (is_correct+1)/2	#increment when is_correct=1
		unit.known_predictions += int(is_correct>=0)	#increment when is_correct=0 or 1
		unit.unknown_predictions += -(is_correct-1)/2	#increment when is_correct = -1
		unit.tuple_list.append([p_value, is_correct])	#the tuple_list where the p_value_cut_off comes
	
	def prediction_space2p_value_cut_off(self, accuracy_cut_off):
		"""
		02-28-05
			for each go_no, for each prediction_space, compute its p_value_cut_off
			
			--p_value_outof_accuracy_cut_off()
		"""
		for go_no in self.go_no2prediction_space:
			for (prediction_space, attr) in self.go_no2prediction_space[go_no].iteritems():
				attr.p_value_cut_off = self.p_value_outof_accuracy_cut_off(attr.tuple_list, accuracy_cut_off)
			
	def p_value_dictionary_setup(self, prediction_list):
		"""
		02-28-05
			borrowed from gene_stat_plot.py
		"""
		p_value_dictionary = {}
		for ls in prediction_list:
			p_value = ls[0]
			if p_value in p_value_dictionary:
				p_value_dictionary[p_value].append( ls[1:])
			else:
				p_value_dictionary[p_value] = [ls[1:]]
		return p_value_dictionary

	def p_value_outof_accuracy_cut_off(self, tuple_list, accuracy_cut_off):
		"""
		02-28-05
			also borrowed from gene_stat_plot.py,
			change the column no. of is_correct from 3 to 0
		
			--p_value_dictionary_setup()
		"""
		if self.debug:
			print "\t\t###Enter p_value_outof_accuracy_cut_off()###"
		total_clusters = 0.0
		good_clusters = 0.0
		#if no good p_value_cut_off to get accuracy, it's to be 0
		p_value_cut_off = 0.0
		p_value2cumu_accuracy = {}
		p_value_dictionary = self.p_value_dictionary_setup(tuple_list)
		p_value_list = p_value_dictionary.keys()
		p_value_list.sort()
		if self.debug:
			print 'p_value_list sorted:%s'%repr(p_value_list)
			
		for p_value in p_value_list:
			prediction_array = array(p_value_dictionary[p_value])
			#is_correct = 0 or 1 is a prediction of known genes
			no_of_clusters = sum(greater_equal(prediction_array[:,0], 0))
			#is_correct = 1 is a good prediction
			no_of_good_clusters = sum(greater(prediction_array[:,0],0))
			total_clusters += no_of_clusters
			good_clusters += no_of_good_clusters
			
			if total_clusters == 0:
				p_value2cumu_accuracy[p_value] = -1
				if good_clusters !=0:
					sys.stderr.write("ERROR: good_clusters are more than total_clusters \n")
			else:
				p_value2cumu_accuracy[p_value] = good_clusters/total_clusters
			if self.debug:
				print "p_value: %s"%p_value
				print "prediction_array: %s"%prediction_array
				print "no_of_clusters: %s"%no_of_clusters
				print "no_of_good_clusters: %s"%no_of_good_clusters
				print "total_clusters:%s"%total_clusters
				print "good_clusters: %s"%good_clusters
				print "accuracy: %s"%(good_clusters/total_clusters)
				raw_input("pause:")
		p_value_list.reverse()
		for p_value in p_value_list:
			if p_value2cumu_accuracy[p_value] >= accuracy_cut_off:
				p_value_cut_off = p_value
				if self.debug:
					print "accuracy: %s"%p_value2cumu_accuracy[p_value]
				break
		if self.debug:
			print "p_value_cut_off:  %s"%p_value_cut_off
			print "\t\t##leave p_value_outof_accuracy_cut_off()"
		return p_value_cut_off
	
	def p_value_cut_off_output(self, outf):
		"""
		02-28-05
			output the self.go_no2prediction_space
		"""
		writer = csv.writer(outf, delimiter='\t')
		for (go_no, prediction_space2attr) in self.go_no2prediction_space.iteritems():
			writer.writerow(['go_no', go_no])
			writer.writerow(["recurrence", "connectivity", "p_value", "correct_predictions", "known_predictions", "unknown_predictions"])
			for (prediction_space, attr) in prediction_space2attr.iteritems():
				writer.writerow([prediction_space[0], prediction_space[1], attr.p_value_cut_off, \
				attr.correct_predictions, attr.known_predictions, attr.unknown_predictions])
	
	def lm_fit(self, lm_instance, curs=None, lm_table=None):
		"""
		02-28-05
			linear model fitting here
			
			--data_prepare
			--submit
		"""
		for go_no in self.go_no2prediction_space:
			y_list, x_2d_list = self.data_prepare(self.go_no2prediction_space[go_no])
			if len(y_list) >=self.min_data_points and y_list!=[0]*len(y_list):
				#all 0 in y_list means accuracy can't be achieved for this go_no
				lm_instance.prepare_data(y_list, x_2d_list)
				lm_instance.run()
				coeff_list = lm_instance.coefficients()
				chisq_tuple = lm_instance.chisq_return()
				lm_instance.cleanup()
				#append the chisq to coeff_list
				coeff_list.append(chisq_tuple[0])
				self.go_no2lm_results[go_no] = coeff_list
				
	def data_prepare(self, prediction_space2attr):
		"""
		02-28-05
			prepare y_list(p_value_cut_off list) and (1, recurrence, connectivity) matrix, x_2d_list
			
		"""
		y_list = []
		x_2d_list = []
		for (prediction_space, attr) in prediction_space2attr.iteritems():
			if attr.known_predictions >= self.valid_space:
				y_list.append(attr.p_value_cut_off)
				x_2d_list.append([1, prediction_space[0], prediction_space[1]])		#(1, recurrence, connectivity)
		return (y_list, x_2d_list)
	
	def lm_results_output(self, outf, go_no2lm_results):
		"""
		02-28-05
			output the self.go_no2lm_results
		"""
		writer = csv.writer(outf, delimiter='\t')
		writer.writerow(['go_no', 'intercept', 'coeff1', 'coeff2', 'chisq'])
		for (go_no, coeff_list) in go_no2lm_results.iteritems():
			writer.writerow([go_no, coeff_list[0], coeff_list[1], coeff_list[2], coeff_list[3]])
	
	def lm_table_create(self, curs, lm_table):
		try:
			curs.execute("create table %s(\
				go_no	integer,\
				intercept	float,\
				coeff1	float,\
				coeff2	float,\
				chisq	float)"%lm_table)
		except:
			sys.stderr.write("Error occurred when creating table %s\n"%lm_table)	

	def submit(self, curs, lm_table, go_no2lm_results):
		"""
		02-28-05
			submit the linear model parameters to the database
		"""
		sys.stderr.write("Database transacting...")
		try:
			for (go_no, coeff_list) in go_no2lm_results.iteritems():
				curs.execute("insert into %s(go_no, intercept, coeff1, coeff2, chisq)\
					values (%d, %f, %f, %f, %f)"%\
					(lm_table, go_no, coeff_list[0], coeff_list[1], coeff_list[2], coeff_list[3]) )
		except:
			sys.stderr.write('Error occurred when inserting pattern. Aborted.\n')
			sys.exit(1)
		sys.stderr.write("done.\n")

	def run(self):
		"""
		02-24-05
			--init
			--db_connect
			--data_fetch
			--lm_fit
		
		"""
		#some additional initialization
		self.init()
		
		(conn, curs) = db_connect(self.hostname, self.dbname, self.schema)
		curs.execute("begin")
		if self.needcommit:
			if self.lm_table:
				self.lm_table_create(curs, self.lm_table)
			else:
				sys.stderr.write("Please the lm_table to commit.\n")
				sys.exit(127)
			
		self.data_fetch(curs, self.table)
		self.prediction_space2p_value_cut_off(self.accuracy_cut_off)
		self.p_value_cut_off_output(sys.stdout)
		self.lm_fit(self.lm_instance)
		self.lm_results_output(sys.stdout, self.go_no2lm_results)
		if self.needcommit:
			self.submit(curs, self.lm_table, self.go_no2lm_results)
			curs.execute("end")

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)	
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:l:a:j:m:v:cru", ["help", "hostname=", \
			"dbname=", "schema=", "table=", "lm_table=", "accuracy_cut_off=", "judger_type=",\
			"min_data_points=", "valid_space=", "commit", "report", "debug"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = None
	lm_table = None
	accuracy_cut_off = 0.5
	judger_type = 0
	min_data_points = 5
	valid_space = 20
	commit = 0
	report = 0
	debug = 0
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
		elif opt in ("-l", "--lm_table"):
			lm_table = arg
		elif opt in ("-a", "--accuracy_cut_off"):
			accuracy_cut_off = float(arg)
		elif opt in ("-j", "--judger_type"):
			judger_type = int(arg)
		elif opt in ("-m", "--min_data_points"):
			min_data_points = int(arg)
		elif opt in ("-v", "--valid_space"):
			valid_space = int(arg)
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-u", "--debug"):
			debug = 1
	if schema and table:
		instance = p_gene_lm(hostname, dbname, schema, table, lm_table, accuracy_cut_off,\
			judger_type, min_data_points, commit, report, debug, valid_space)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
