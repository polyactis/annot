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
	-l ,,,, --lm_table=...	the lm_table to store the linear_model results, needed if needcommit
	-n ..., --gene_p_table=...	the table to store p_gene_id:p_value_cut_off, which are real predictions
	-p ..., --p_value_cut_off=...	p_value_cut_off. 0 means using the linear model from lm_table
	-j ..., --judger_type=...	how to judge predicted functions, 0(default), 1, 2
	-x ..., --recurrence_gap_size=...	2(default)
	-y ..., --connectivity_gap_size=...	2(default)
	-r, --report	report the progress(a number)
	-c, --commit	commit the database transaction, records down the go_no2accuracy.
	-a ..., --accuracy_cut_off=...	the accuracy_cut_off to be based for p_value adjusting, 0(default)
		NOTICE: 0 means using p_value_cut_off method
	-b, --debug	enable debugging, no debug by default
	-h, --help              show this help

Examples:
	#Use the p_value_cut_off
		p_gene_analysis.py -k sc_54 -g p_gene_e5 -p 0.001 stat_table
		p_gene_analysis.py -k sc_54 -p 0.01 -c -j 2  -g p_gene_cluster_stat2 stat_table -n gene_p
	
	
	#Don't use p_value_cut_off, get the linear model from table lm_p_gene_repos_2_e5_v40
		p_gene_analysis.py -k sc_54 -p 0 -j 1 -g p_gene_repos_2_e5
			-l lm_p_gene_repos_2_e5_v40 stat_table
	
	#Don't use p_value_cut_off, get the linear model from table lm_p_gene_repos_2_e5_v40
	#commit the changes to table gene_p_repos_2_e5.
		p_gene_analysis.py -k sc_54 -p 0 -j 1 -g p_gene_repos_2_e5 -c -n gene_p_repos_2_e5
			-l lm_p_gene_repos_2_e5_v40 stat_table

Description:
	02-28-05
	This is the second part of gene_stat_plot.py, which works on p_gene table(overhauled)
	The first part goes to gene_stat.py.
	Input is p_gene table. Output includes a file include some stats and a gene_p table(if commit)
	storing the good predictions under the judger_type and other factors(p-value or linear_model)
	03-01-05
		add a functionality to get p_value_cut_off from linear_model
	03-01-05
		add the gene-p table to record real predictions (p_gene_id: p_value_cut_off)
"""

import sys, os, psycopg, getopt, csv, math
from sets import Set
from numarray import *
from codense.common import db_connect

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
	02-28-05
		alter the program towards reusage
	02-28-05
		this is the second part of gene_stat_plot.py, which works on p_gene table
		the first part goes to gene_stat.py
	03-01-05
		_p_gene_analysis(): add a functionality to get p_value_cut_off from linear_model
	03-02-05
		submit(): gene_p table enlarged, p_gene_id_src is for gene_p_map_redundancy.py
	03-08-05
		add two more parameters, recurrence_gap_size and connectivity_gap_size (make them explicit)

	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, \
		table=None, mcl_table=None, p_value_cut_off=0.01, report=0, \
		judger_type=0, needcommit=0, gene_table='p_gene', lm_table=None, \
		stat_table_fname=None, debug=0, accuracy_cut_off=0, gene_p_table=None,\
		recurrence_gap_size=2, connectivity_gap_size=2):
		"""
		03-08-05
			p_value_cut_off=0.01, otherwise float(None) doesn't work.
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema		
		self.table = table
		self.mcl_table = mcl_table
		self.p_value_cut_off = float(p_value_cut_off)
		self.report = int(report)
		self.judger_type = int(judger_type)
		self.needcommit = int(needcommit)
		self.gene_table = gene_table
		self.lm_table = lm_table
		self.stat_table_fname = stat_table_fname
		#debugging flag
		self.debug = int(debug)
		self.accuracy_cut_off = float(accuracy_cut_off)
		self.gene_p_table = gene_p_table
		#the gap between two recurrences
		self.recurrence_gap_size = int(recurrence_gap_size)
		self.connectivity_gap_size = int(connectivity_gap_size)
		
	def init(self):	
		#open a file
		if self.stat_table_fname:
			self.stat_table_f = csv.writer(open(self.stat_table_fname, 'w'), delimiter='\t')

		#an is_correct dictionary used in database fetch
		self.is_correct_dict = {0: 'is_correct',
			1: 'is_correct_L1',
			2: 'is_correct_lca'}
		
		#a counter
		self.no_of_records = 0
		
		#data structures to be created
		self.go_no2prediction_space = {}
		self.prediction_space2attr = {}
		self.prediction_pair2attr = {}
		self.go_no2accuracy = {}
		self.go_no2lm_results = {}
		self.general_lm_results = None
		self.gene_p_list = []
		
	def get_known_genes_dict(self, curs):
		sys.stderr.write("Getting known_genes_dict...")
		known_genes_dict = {}
		curs.execute("select gene_no,go_functions from gene where known=TRUE")
		rows = curs.fetchall()
		for row in rows:
			go_functions_list = row[1][1:-1].split(',')
			known_genes_dict[row[0]] = Set()
			for go_no in go_functions_list:
				known_genes_dict[row[0]].add(int(go_no))
		sys.stderr.write("Done\n")
		return known_genes_dict

	def get_go_no2lm_results(self, curs, lm_table):
		"""
		03-01-05
			for each go_no, get its linear model parameters
			also return 2d matrix for a general linear model(see get_general_lm_results)
		03-27-05
			p_gene_lm.py changed the lm_table format.(logistic regression)
		07-01-05
			add coeff4
		"""
		sys.stderr.write("Getting linear model parameters...")
		go_no2lm_results = {}
		lm_results_2d_list = []
		curs.execute("select go_no, intercept, coeff1, coeff2, coeff3, coeff4, score_cut_off from %s"%lm_table)
		rows = curs.fetchall()
		for row in rows:
			go_no = row[0]
			coeff_list = row[1:]		#intercept coeff1 coeff2 coeff3 coeff4 score_cut_off
			go_no2lm_results[go_no] = coeff_list
			lm_results_2d_list.append(coeff_list)
		sys.stderr.write("Done\n")		
		return (go_no2lm_results, lm_results_2d_list)
		
	def get_general_lm_results(self, lm_results_2d_list):
		"""
		03-01-05
			return (intercept, coeff1, coeff2) by averaging the known parameters
		03-27-05
			p_gene_lm.py changed the lm_table format.(logistic regression)
		07-01-05
			add coeff4
		"""
		no_of_models = len(lm_results_2d_list)
		lm_results_array = array(lm_results_2d_list)
		intercept = sum(lm_results_array[:,0])/no_of_models
		coeff1 = sum(lm_results_array[:,1])/no_of_models
		coeff2 = sum(lm_results_array[:,2])/no_of_models
		coeff3 = sum(lm_results_array[:,3])/no_of_models
		coeff4 = sum(lm_results_array[:,4])/no_of_models
		score_cut_off = sum(lm_results_array[:,-1])/no_of_models
		return (intercept, coeff1, coeff2, coeff3, coeff4, score_cut_off)
	
	def prediction_accepted(self, go_no, property_list):
		"""
		03-01-05
			if go_no in self.go_no2lm_results, go and get it
			otherwise, use the general linear model
		03-27-05
			change to return (is_accepted, score)
		07-01-05
			add coeff4
		"""
		if go_no in self.go_no2lm_results:
			lm_results = self.go_no2lm_results[go_no]
			#intercept + coeff1*(-lg(p_value)) + coeff2*recurrence + coeff3*connectivity
			score = lm_results[0] + lm_results[1]*property_list[0] + lm_results[2]*property_list[1] +\
				lm_results[3]*property_list[2] + lm_results[4]*property_list[3]
			is_accepted = (score>=lm_results[-1])
		else:
			score = self.general_lm_results[0] + self.general_lm_results[1]*property_list[0]+\
				self.general_lm_results[2]*property_list[1]+self.general_lm_results[3]*property_list[2] +\
				self.general_lm_results[4]*property_list[3]
			is_accepted = (score>=self.general_lm_results[-1])
		return (is_accepted, score)
		


	def data_fetch(self, curs, gene_table):
		"""
		02-21-05
			--_p_gene_analysis()
		07-01-05
			add cluster_size_cut_off in fetch
		"""
		sys.stderr.write("Setting up prediction_space and prediction_pair...\n")
		curs.execute("DECLARE crs CURSOR FOR select gene_no, go_no, mcl_id, %s, avg_p_value, \
			recurrence_cut_off,connectivity_cut_off, depth_cut_off,p_gene_id, cluster_size_cut_off from %s"\
			%(self.is_correct_dict[self.judger_type], gene_table))
		
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		while rows:
			for row in rows:
				self._p_gene_analysis(row)	
				self.no_of_records+=1

			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, self.no_of_records))
			
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		sys.stderr.write("Done\n")		
	
	def _p_gene_analysis(self, row):
		"""
		02-28-05
			setup two dictionaries, self.prediction_space2attr and self.prediction_pair2attr
		03-01-05
			use the hardcode p_value_cut_off or get it from linear_model
			depend on whether self.p_value_cut_off==0
		03-06-05
			seperate the data of different function categories, into go_no2prediction_space
		03-08-05
			Take the floor of recurrence and connectivity
		03-27-05
			is_accepted is used to judge whether the prediction is good or not.
		07-01-05
			add cluster_size
		"""
		gene_no = row[0]
		go_no = row[1]
		mcl_id = row[2]
		is_correct = row[3]
		p_value = row[4]
		recurrence = row[5]
		connectivity = row[6]
		depth_cut_off = row[7]
		p_gene_id = row[8]
		cluster_size = row[9]
		
		if self.p_value_cut_off == 0:
			#the model is based on -log(p_value).
			if p_value ==0:
				p_value = 1e-8
			(is_accepted, score) = self.prediction_accepted(go_no, [-math.log(p_value), recurrence, connectivity, cluster_size])
		else:
			is_accepted = (p_value <= self.p_value_cut_off)
			score = p_value
		
		if is_accepted:
			self.gene_p_list.append([p_gene_id, score])
		else:
			return

		#take the floor of the recurrence
		recurrence = int(math.floor(recurrence/self.recurrence_gap_size)*self.recurrence_gap_size)
		#take the floor of the connectivity *10
		connectivity = int(math.floor(connectivity*10/self.connectivity_gap_size)*self.connectivity_gap_size)
		
		if go_no not in self.go_no2prediction_space:
			self.go_no2prediction_space[go_no] = {}
		#pass the value to ease programming
		local_prediction_space2attr = self.go_no2prediction_space[go_no]
		
		prediction_space = (recurrence, connectivity)
		prediction_pair = (gene_no, go_no)
		
		"""
		0. deal with local_prediction_space2attr
		"""
		if prediction_space not in local_prediction_space2attr:
			#initalize its attribute if it doesn't exist
			local_prediction_space2attr[prediction_space] = prediction_space_attr()
		#pass the value to a value, ease programming
		unit = local_prediction_space2attr[prediction_space]
		unit.correct_predictions += (is_correct+1)/2	#increment when is_correct=1
		unit.known_predictions += int(is_correct>=0)	#increment when is_correct=0 or 1
		unit.unknown_predictions += -(is_correct-1)/2	#increment when is_correct = -1
		unit.tuple_list.append([p_value, is_correct])	#the tuple_list where the p_value_cut_off comes

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

	def overview_stats(self, stat_table_f):
		'''
		02-21-05
			get some global statistics
			--return_known_unknown_gene_sets()
		'''
		sys.stderr.write("Outputing overview stats...")
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
		stat_table_f.writerow(['accuracy', total_accuracy])
		stat_table_f.writerow(['known-predictions', total_known_predictions])
		stat_table_f.writerow(['unknown-predictions', total_unknown_predictions])
		stat_table_f.writerow(['accuracy-pair', total_accuracy_pair])
		stat_table_f.writerow(['known-predictions-pair', total_known_predictions_pair])
		stat_table_f.writerow(['unknown-predictions-pair', total_unknown_predictions_pair])
		
		self.stat_table_f.writerow(['known genes', no_of_known_genes])
		self.stat_table_f.writerow(['unknown genes', no_of_unknown_genes])
		sys.stderr.write("Done\n")		
	
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
	

	def go_no_accuracy(self, prediction_pair2attr, stat_table_f, curs):
		"""
		this function computes the prediction accuracy for each go function
		"""
		sys.stderr.write("Computing and outputing go_no accuracy...")
		for prediction_pair in prediction_pair2attr:
			go_no = prediction_pair[1]
			is_correct = prediction_pair2attr[prediction_pair].is_correct
			if go_no not in self.go_no2accuracy:
				self.go_no2accuracy[go_no] = accuracy_struc()
			self.go_no2accuracy[go_no].correct += (is_correct+1)/2	#increment when is_correct=1
			self.go_no2accuracy[go_no].known += int(is_correct>=0)	#increment when is_correct=0 or 1
			self.go_no2accuracy[go_no].unknown += -(is_correct-1)/2	#increment when is_correct = -1
		
		"""
		if self.needcommit:
			#this table is easy for linking.
			go_acc_table = '%s_acc_%s'%(self.gene_table, self.is_correct_dict[self.judger_type])
			curs.execute("create table %s(\
				go_no	integer,\
				accuracy	float,\
				known	integer,\
				unknown	integer\
				)"%go_acc_table)
		"""
		
		#finally, compute all the ratios and output
		go_no_list = self.go_no2accuracy.keys()
		go_no_list.sort()
		stat_table_f.writerow([])
		stat_table_f.writerow(['go_no', 'accuracy', 'known', 'unknown'])
		for go_no in go_no_list:
			unit = self.go_no2accuracy[go_no]
			if unit.known !=0:
				unit.ratio = unit.correct/unit.known
			else:
				unit.ratio = -1
			stat_table_f.writerow([go_no, unit.ratio, unit.known, unit.unknown])
			"""
			if self.needcommit:
				curs.execute("insert into %s(go_no, accuracy, known, unknown)\
					values(%s, %f, %d, %d)"%(go_acc_table,\
					go_no, unit.ratio, unit.known, unit.unknown))
			"""
		sys.stderr.write("Done\n")		

	def table_output(self, stat_table_f, prediction_space2attr):
		'''
		
		'''

		header = ['recurrence', 'connectivity', 'accuracy', 'known_predictions', 'unknown_predictions',\
		'accuracy_pair', 'known_predictions_pair', 'unknown_predictions_pair', 'known genes', 'unknown genes']
		stat_table_f.writerow(header)
		prediction_space_list = prediction_space2attr.keys()
		prediction_space_list.sort()
		for prediction_space in prediction_space_list:
			unit = prediction_space2attr[prediction_space]
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
			stat_table_f.writerow(row)
		
	def gene_p_table_submit(self, curs, gene_p_table, gene_p_list):
		"""
		03-02-05
			gene_p table enlarged, p_gene_id_src is for gene_p_map_redundancy.py
		"""
		sys.stderr.write("Submiting gene_p_table...")
		curs.execute("create table %s(\
			gene_p_id	serial primary key,\
			p_gene_id	integer,\
			p_value_cut_off	float,\
			p_gene_id_src integer)"%gene_p_table)
			
		for (p_gene_id, p_value_cut_off) in gene_p_list:
			curs.execute("insert into %s(p_gene_id, p_value_cut_off)\
					values(%d, %f)"%(gene_p_table, p_gene_id, p_value_cut_off))
		sys.stderr.write("done.\n")


	#######03-07-05 a bunch of functions related to two posterior maneuvering of 
	#go_no2prediction_space, grouping and accumulatiing.
	#See log of 2005, section 'linear model overfitting' for detail.
	###begin 

	def return_go_no_group2prediction_space(self, go_no2prediction_space, curs, distance_table):
		"""
		03-06-05
			input: go_no2prediction_space, curs
			output: go_no_group2prediction_space
			group the go_no2prediction_space, data of several go_nos who are parent-child
			will be merged.
		"""
		go_no_list = go_no2prediction_space.keys()
		go_no_map = self.return_go_no_map(go_no_list, curs, distance_table)
		go_no_groups = self.dict_map2group(go_no_map)
		go_no_group2prediction_space = {}
		sys.stderr.write("Grouping go_no2prediction_space data...")
		for go_no_group in go_no_groups:
			unit = self.return_data_of_the_group(go_no_group, go_no2prediction_space)
			key = tuple(go_no_group)
			go_no_group2prediction_space[key] = unit
		sys.stderr.write("done.\n")
		return go_no_group2prediction_space
	
	def return_go_no_map(self, go_no_list, curs, distance_table):
		"""
		03-06-05
			input: a list of go_nos, curs
			output: a map showing which go_no corresponds to which
			
			curs is used to get the go_no2term_id and nodes pairwise distance
		"""
		sys.stderr.write("Mapping go_nos...")
		from gene_p_map_redundancy import gene_p_map_redundancy
		from codense.common import get_go_no2term_id
		borrowed_instance = gene_p_map_redundancy()
		go_no_map = {}
		go_no2term_id = get_go_no2term_id(curs)
		go_no2distance = {}
		
		for i in range(len(go_no_list)):
			go_no = go_no_list[i]
			if go_no not in go_no_map:
				#not flagged, map
				go_no_map[go_no] = go_no
				for j in range(i+1, len(go_no_list)):
					go_no2 = go_no_list[j]
					if go_no < go_no2:
						key= (go_no, go_no2)
					else:
						key = (go_no2, go_no)
					if key in go_no2distance:
						jasmine_distance = go_no2distance[key][2]
					else:
						jasmine_distance = borrowed_instance.get_distance(curs, go_no, go_no2, distance_table, go_no2distance, go_no2term_id)
					if jasmine_distance == 0:
						#jasmine_distance=0 means they are parent-child
						go_no_map[go_no2] = go_no
		sys.stderr.write("done.\n")
		return go_no_map	
		
	def dict_map2group(self, go_no_map):
		"""
		03-06-05
			input: a map of go_nos
			output: a list of lists, each list consists of go_nos who are parent-child
		"""
		go_no_reverse_map = {}
		for (go_no, go_no_src) in go_no_map.iteritems():
			if go_no_src not in go_no_reverse_map:
				go_no_reverse_map[go_no_src] = [go_no]
			else:
				go_no_reverse_map[go_no_src].append(go_no)
		return go_no_reverse_map.values()
		
	
	def return_data_of_the_group(self, go_no_list, go_no2prediction_space):
		"""
		03-06-05
			input: a list of go_nos, go_no2prediction_space
			output: the combined data(prediction_space_attr) of the go_nos from the go_no_list
		"""
		local_prediction_space2attr = {}
		for go_no in go_no_list:
			for (prediction_space,unit) in go_no2prediction_space[go_no].iteritems():
				if prediction_space not in local_prediction_space2attr:
					local_prediction_space2attr[prediction_space] = prediction_space_attr()
				#pass the value to a value, ease programming
				local_unit = local_prediction_space2attr[prediction_space]
				local_unit.correct_predictions += unit.correct_predictions
				local_unit.known_predictions += unit.known_predictions
				local_unit.unknown_predictions += unit.unknown_predictions
				local_unit.tuple_list += unit.tuple_list
		return local_prediction_space2attr

	
	def return_cumulative_prediction_space2attr(self, prediction_space2attr, recurrence_gap_size=2, connectivity_gap_size=2):
		"""
		03-06-05
			input: prediction_space2attr
			output: a cumulative version of prediction_space2attr
		"""
		cumulative_prediction_space2attr = {}
		prediction_space_keys = prediction_space2attr.keys()
		prediction_space_array = array(prediction_space_keys)
		min_recurrence = int(min(prediction_space_array[:,0]))
		max_recurrence = int(max(prediction_space_array[:,0]))
		min_connectivity = int(min(prediction_space_array[:,1]))
		max_connectivity = int(max(prediction_space_array[:,1]))
		for recurrence in range(min_recurrence, max_recurrence+recurrence_gap_size, recurrence_gap_size):
			for connectivity in range(min_connectivity, max_connectivity+connectivity_gap_size, connectivity_gap_size):
				prediction_space = (recurrence, connectivity)
				cumulative_prediction_space2attr[prediction_space] = prediction_space_attr()
				#pass the value to ease programming
				cumulative_unit = cumulative_prediction_space2attr[prediction_space] 
				for real_prediction_space in prediction_space_keys:
					if real_prediction_space[0]>=prediction_space[0] and real_prediction_space[1]>=prediction_space[1]:
						unit = prediction_space2attr[real_prediction_space]
						cumulative_unit.correct_predictions += unit.correct_predictions
						cumulative_unit.known_predictions += unit.known_predictions
						cumulative_unit.unknown_predictions += unit.unknown_predictions
						cumulative_unit.tuple_list += unit.tuple_list
		return cumulative_prediction_space2attr
	
	######end
	#######################
	
	def prediction_space_split_output(self, stat_table_f, go_no2prediction_space, recurrence_gap_size=2, connectivity_gap_size=2):
		"""
		03-08-05
			Input: go_no2prediction_space
			Output the prediction_space accuracy information go_no by go_no
		"""
		sys.stderr.write("Outputing prediction_space accuracy information go_no by go_no...")
		#make a blank line
		stat_table_f.writerow([])
		for (go_no, unit) in go_no2prediction_space.iteritems():
			stat_table_f.writerow(['go_no:', go_no])
			self.table_output(stat_table_f, unit)
			#cumulative form
			stat_table_f.writerow(["cumulative form"])
			cumulative_unit = self.return_cumulative_prediction_space2attr(unit, recurrence_gap_size, connectivity_gap_size)
			self.table_output(stat_table_f, cumulative_unit)
		sys.stderr.write("Done.\n")

	def run(self):
		"""
		02-28-05
		
		03-07-05
			implementing two posterior maneuvering of go_no2prediction_space, grouping and accumulatiing.
			See log of 2005, section 'linear model overfitting' for detail.
		
		--init()
		--db_connect()
		--IF self.p_value_cut_off==0
			--get_go_no2lm_results
			--get_general_lm_results
		
		--data_fetch()
			--_p_gene_analysis()
		--IF self.stat_table_fname
			--overview_stats()
				--return_known_unknown_gene_sets()
			--go_no_accuracy()
			--table_output()
		--IF self.gene_p_table
			--gene_p_table_submit()
		"""
		self.init()
		(conn, curs) = db_connect(self.hostname, self.dbname, self.schema)
		self.known_genes_dict = self.get_known_genes_dict(curs)
		if self.p_value_cut_off == 0:
			if self.lm_table:
				self.go_no2lm_results, lm_results_2d_list = self.get_go_no2lm_results(curs, self.lm_table)
				self.general_lm_results = self.get_general_lm_results(lm_results_2d_list)
			else:
				sys.stderr.write("p_value_cut_off==0, need the lm_table to get the linear model\n")
				sys.exit(127)
		
		self.data_fetch(curs, self.gene_table)
		if self.stat_table_fname:
			self.overview_stats(self.stat_table_f)
			self.go_no_accuracy(self.prediction_pair2attr, self.stat_table_f, curs)
			self.table_output(self.stat_table_f, self.prediction_space2attr)
			"""
			#first grouping the data of parent-child go functions
			distance_table = 'go.node_dist'
			go_no_group2prediction_space = self.return_go_no_group2prediction_space(self.go_no2prediction_space, curs, distance_table)
			#output the prediction_space go_no by go_no
			self.prediction_space_split_output(self.stat_table_f, go_no_group2prediction_space, self.recurrence_gap_size, self.connectivity_gap_size)
			"""

		if self.gene_p_table:
			self.gene_p_table_submit(curs, self.gene_p_table, self.gene_p_list)
		if self.needcommit:
			curs.execute("end")

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema=", "table=", "mcl_table=", "p_value_cut_off=",\
		"judger_type=", "report", "commit", "gene_table=", "lm_table=", "debug", "accuracy_cut_off=",\
		"gene_p_table=",  "recurrence_gap_size=", "connectivity_gap_size="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:m:p:j:rcg:l:ba:n:x:y:", long_options_list)
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
	lm_table = None
	debug = 0
	accuracy_cut_off = 0
	stat_table_fname = None
	gene_p_table = None
	recurrence_gap_size = 2
	connectivity_gap_size = 2
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
		elif opt in ("-l", "--lm_table"):
			lm_table = arg
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-a", "--accuracy_cut_off="):
			accuracy_cut_off = float(arg)
		elif opt in ("-n", "--gene_p_table="):
			gene_p_table = arg
		elif opt in ("-x", "--recurrence_gap_size"):
			recurrence_gap_size = int(arg)
		elif opt in ("-y", "--connectivity_gap_size"):
			connectivity_gap_size = int(arg)
	
	if len(args) == 1:
		stat_table_fname = args[0]
			
	if schema and p_value_cut_off!=None and stat_table_fname:
		instance = p_gene_analysis(hostname, dbname, schema, table, mcl_table, p_value_cut_off,\
			report, judger_type, commit, gene_table, lm_table, stat_table_fname, debug, \
			accuracy_cut_off, gene_p_table, recurrence_gap_size, connectivity_gap_size)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
