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


Description:

"""

import sys, os, psycopg, getopt, csv, math
from sets import Set
from rpy import r
from numarray import array
from numarray import searchsorted
from codense.common import db_connect

class prediction_space_attr:
	"""
	data structure for prediction_space2attr
	"""
	def __init__(self):
		self.is_correct_list = []
		self.is_correct_ratio_list = []
		self.p_gene_id_list = []
		self.mcl_id_set = Set()
		self.prediction_pair_dict = {}
		self.prediction_pair_ratio_dict = {}
	
	def get_average(self, list):
		return sum(list)/float(len(list))
	
	def return_is_correct_average(self):
		return self.get_average(self.prediction_pair_dict.values())
	
	def return_is_correct_ratio_average(self):
		return self.get_average(self.prediction_pair_ratio_dict.values())

def get_average(array_2d, which_column):
	"""
	03-31-05
	"""
	array_1d = array_2d[:,which_column]
	return sum(array_1d)/float(len(array_1d))

class p_gene_factor:
	"""
	03-30-05
		study the effect of different factors on prediction accuracy.
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


		#an is_correct dictionary used in database fetch
		self.is_correct_dict = {0: 'is_correct',
			1: 'is_correct_L1',
			2: 'is_correct_lca'}
		
		#more gap sizes
		self.p_value_gap_size = 10
		self.cluster_size_gap_size = 5
		self.unknown_gap_size = 2
		
		#a counter
		self.no_of_records = 0
		
		#data structures to be created
		self.go_no2prediction_space = {}
		self.prediction_space2attr = {}
		
		self.prediction_data = []
	
	def data_fetch(self, curs, gene_table):
		"""
		02-21-05
			--_p_gene_analysis()
		"""
		sys.stderr.write("Setting up prediction_space and prediction_pair...\n")
		curs.execute("DECLARE crs CURSOR FOR select gene_no, go_no, %s, avg_p_value, \
			recurrence_cut_off,connectivity_cut_off, cluster_size_cut_off, unknown_cut_off, \
			depth_cut_off,lca_list,p_gene_id,mcl_id from %s"\
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
		"""
		gene_no = row[0]
		go_no = row[1]
		is_correct = row[2]
		p_value = row[3]
		recurrence = row[4]
		connectivity = row[5]
		cluster_size = row[6]
		unknown = row[7]
		depth_cut_off = row[8]
		lca_list = row[9]
		p_gene_id = row[10]
		mcl_id = row[11]
		
		if lca_list and is_correct ==1:
			lca_list = lca_list[1:-1].split(',')
			lca_list = map(int,lca_list)
			#Xifeng's idea
			is_correct_ratio = float(self.go_no2depth[lca_list[0]]-depth_cut_off+1)/(self.go_no2depth[go_no]-depth_cut_off+1)
		elif is_correct == 0:
			is_correct_ratio = is_correct
		else:
			#no unknown genes
			return
		
		if p_value ==0:
			p_value = 1e-8
		#-log transforms p_value
		p_value = -math.log(p_value)
		
		
		self.prediction_data.append([gene_no, go_no, is_correct, is_correct_ratio, p_value, recurrence, \
			connectivity, cluster_size, unknown, p_gene_id, mcl_id])

		
	def group_data(self, data_list_2d, key_column=0, no_of_groups=6, cluster_column=-1):
		"""
		03-30-05
			output: a dictionary
			group the data based on the key_column, but each key has similar amount of clusters
			from cluster_column. idea is similar to equal.count().
		"""
		data_array = array(data_list_2d)
		cluster_list = list(data_array[:,cluster_column])
		cluster_set = Set(cluster_list)
		unit_length = len(cluster_set)/no_of_groups
		
		key2cluster_set = {}
		for i in range(len(data_array)):
			key = data_array[i,key_column]
			cluster_id = data_array[i, cluster_column]
			if key not in key2cluster_set:
				key2cluster_set[key] = Set()
			key2cluster_set[key].add(cluster_id)
		
		key_cluster_2d_list = []
		for key,cluster_set in key2cluster_set.iteritems():
			key_cluster_2d_list.append([key,cluster_set])
		key_cluster_2d_list.sort()
		
		bin_boundaries = [key_cluster_2d_list[0][0]]	#first key is already pushed in
		bin_set = Set()
		for key,cluster_set in key_cluster_2d_list:
			if len(bin_set) < unit_length:
				#the limit hasn't been reached.
				bin_set |= cluster_set
			else:
				#restart
				bin_boundaries.append(key)
				bin_set = cluster_set
		
		key2data_array = {}
		for entry in data_list_2d:
			key = entry[key_column]
			bin_key_index = searchsorted(bin_boundaries, key)	#the trick is that the range is (...]. my bin_boundaries is [...)
													#transform below
			if bin_key_index==len(bin_boundaries):
				bin_key = bin_boundaries[bin_key_index-1]
			elif bin_boundaries[bin_key_index] == key:
				bin_key = key
			else:
				bin_key = bin_boundaries[bin_key_index-1]
			if bin_key not in key2data_array:
				key2data_array[bin_key] = []
			key2data_array[bin_key].append(entry)
		
		return key2data_array
	
	def prediction_space_output(self, outf, prediction_space2attr):
		writer = csv.writer(outf, delimiter='\t')
		writer.writerow(['p_value', 'recurrence', 'connectivity', 'cluster_size', 'unknown', 'acc1', 'acc2', 'no', 'mcl_no'])
		for (prediction_space,unit) in prediction_space2attr.iteritems():
			unit_array = array(unit)

			row = list(prediction_space)
			row.append(get_average(unit_array,2))	#average is_correct
			row.append(get_average(unit_array, 3))	#average is_correct_ratio
			row.append(len(unit_array))	#no of predictions
			row.append(len(Set(unit_array[:,-1])))	#no of clusters
			writer.writerow(row)
		del writer
		
	def run(self):
		"""
		03-30-05
		"""
		self.init()
		(conn, curs) = db_connect(self.hostname, self.dbname, self.schema)
		from codense.common import  get_go_no2depth
		self.go_no2depth = get_go_no2depth(curs)
		
		self.data_fetch(curs, self.gene_table)
		local_prediction_space2attr = self.group_data(self.prediction_data,key_column=4, no_of_groups=4)	#4 is p_value
		for key, unit in local_prediction_space2attr.iteritems():
			local_prediction_space2attr_2 = self.group_data(unit, key_column=5, no_of_groups=5)	#5 is recurrence
			for key2, unit2 in local_prediction_space2attr_2.iteritems():
				self.prediction_space2attr[(key,key2)] = unit2
		#open a file
		if self.stat_table_fname:
			stat_table_f = open(self.stat_table_fname, 'w')
			self.prediction_space_output(stat_table_f, self.prediction_space2attr)

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
	p_value_cut_off = 0.001
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
			
	if schema and stat_table_fname:
		instance = p_gene_factor(hostname, dbname, schema, table, mcl_table, p_value_cut_off,\
			report, judger_type, commit, gene_table, lm_table, stat_table_fname, debug, \
			accuracy_cut_off, gene_p_table, recurrence_gap_size, connectivity_gap_size)
		instance.run()
	else:
		print __doc__
		sys.exit(2)