#!/usr/bin/env python
"""
Usage: p_gene_factor.py -k SCHEMA -p P_VALUE_CUT_OFF [OPTION] STAT_TABLE_FILE

Option:
	STAT_TABLE_FILE is the file where the output goes.
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-f ..., --netmine_fname=...	used to construct table names
		if netmine_fname is given, no need to given following three tables
	-t ..., --table=...	splat_table
	-m ..., --mcl_table=...	mcl_result(default), mcl_result table corresponding to above table.
	-g ..., --gene_table=...	table to store the stat results, p_gene(default), needed if commit
	
	-j ..., --judger_type=...	how to judge predicted functions, 0(default), 1, 2
	-w ..., --which_column=...	5(default, recurrence), which column used to group data(see below)
	-s ..., --group_size=...	10(default), the number of clusters for each group
	-r, --report	report the progress(a number) and output mcl_id and p_gene_id
	-b, --debug	enable debugging, no debug by default
	-h, --help              show this help

Examples:
	p_gene_factor.py -k mm_oxi_stress_7t1 -f fmos_7t1g1e2d40q20s200c50z0001c6
		-w 4 /tmp/stat_table_4_5
	
	p_gene_factor.py -k mm_oxi_stress_7t1 -f fmos_7t1g1e2d40q20s200c50z0001c6
		-w 4,6 -s 20,5 /tmp/stat_table_w46_s205
	
Description:
	Columns 5-tuple:
	gene_no, go_no, is_correct, is_correct_ratio, p_value
	recurrence, connectivity, cluster_size, unknown, p_gene_id
	mcl_id
	
	NOTE: p-value is -log(p-value).

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
	06-28-05
		remove lots of deprecated parameters
	06-30-05
		which_column and group_size are strings to be split into lists.
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, table=None, mcl_table=None, \
		report=0, judger_type=0, gene_table='p_gene', stat_table_fname=None, which_column='5', \
		group_size='10', debug=0):
		"""
		03-08-05
			p_value_cut_off=0.01, otherwise float(None) doesn't work.
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.table = table
		self.mcl_table = mcl_table
		self.report = int(report)
		self.judger_type = int(judger_type)
		self.gene_table = gene_table
		self.stat_table_fname = stat_table_fname
		self.which_column_list = which_column.split(',')
		self.which_column_list = map(int, self.which_column_list) 
		self.group_size_list = group_size.split(',')
		self.group_size_list = map(int, self.group_size_list)
		self.debug = int(debug)
	
	def init(self):
		#an is_correct dictionary used in database fetch
		self.is_correct_dict = {0: 'is_correct',
			1: 'is_correct_L1',
			2: 'is_correct_lca'}

		#a counter
		self.no_of_records = 0
		
		#data structures to be created
		self.go_no2prediction_space = {}
		self.prediction_space2attr = {}
		
		self.prediction_data = []
	
	def data_fetch(self, curs, splat_table, mcl_table, gene_table):
		"""
		02-21-05
			--_p_gene_analysis()
		"""
		sys.stderr.write("Setting up prediction_space and prediction_pair...\n")
		"""
		curs.execute("DECLARE crs CURSOR FOR select gene_no, go_no, %s, avg_p_value, \
			recurrence_cut_off,connectivity_cut_off, cluster_size_cut_off, unknown_cut_off, \
			depth_cut_off,lca_list,p_gene_id,mcl_id from %s"\
			%(self.is_correct_dict[self.judger_type], gene_table))
		"""
		curs.execute("DECLARE crs CURSOR FOR select p.gene_no, p.go_no, p.%s, p.avg_p_value, \
			p.recurrence_cut_off,s.connectivity, p.cluster_size_cut_off, p.unknown_cut_off, \
			p.depth_cut_off,p.lca_list,p.p_gene_id,p.mcl_id,s.no_of_edges from %s p, \
			%s s ,%s m\
			where s.splat_id=p.mcl_id and p.mcl_id=m.mcl_id"\
			%(self.is_correct_dict[self.judger_type], gene_table, splat_table, mcl_table))
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
		no_of_edges = row[12]

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

		#connectivity = no_of_edges/float(cluster_size)
		data_row = [gene_no, go_no, is_correct, is_correct_ratio, p_value, recurrence, \
			connectivity, cluster_size, unknown, p_gene_id, mcl_id]
		self.prediction_data.append(data_row)
		#if self.debug:
		#	print data_row
		#	raw_input("pause:")


	def group_data(self, data_list_2d, key_column=0, no_of_groups=6, group_size=None, cluster_column=-1):
		"""
		03-30-05
			output: a dictionary
			group the data based on the key_column, but each key has similar amount of clusters
			from cluster_column. idea is similar to equal.count().
		"""
		sys.stderr.write("Grouping data...")
		data_array = array(data_list_2d)
		cluster_list = list(data_array[:,cluster_column])
		cluster_set = Set(cluster_list)
		if group_size:
			unit_length = group_size
		else:
			unit_length = len(cluster_set)/no_of_groups
		
		#06-28-05 construct a key 2 set of cluster_id(mcl_id)'s
		key2cluster_set = {}
		for i in range(len(data_array)):
			key = data_array[i,key_column]
			cluster_id = data_array[i, cluster_column]
			if key not in key2cluster_set:
				key2cluster_set[key] = Set()
			key2cluster_set[key].add(cluster_id)
		if self.debug:
			print key2cluster_set
			raw_input("pause:")
		
		#06-28-05 convert key2cluster_set to a 2d list. and sort it based on key
		key_cluster_2d_list = []
		for key,cluster_set in key2cluster_set.iteritems():
			key_cluster_2d_list.append([key,cluster_set])
		key_cluster_2d_list.sort()
		if self.debug:
			print key_cluster_2d_list
			raw_input("pause:")
		
		#06-28-05	construct the boundaries for bin's
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
		if self.debug:
			print "The bin_boundaries is ", bin_boundaries
			raw_input("pause:")
		
		#06-28-05	construct the final data structure to return
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
		if self.debug:
			print "key2data_array is ",key2data_array
			raw_input("pause:")
		sys.stderr.write("Done.\n")
		return key2data_array
	
	def prediction_space_output(self, outf, prediction_space2attr):
		"""
		06-28-05
			gene_no, go_no, is_correct, is_correct_ratio, p_value
			recurrence, connectivity, cluster_size, unknown, p_gene_id
			mcl_id
		06-30-05
			if self.report, output the p_gene_id and mcl_id set.
			
			global structures used:
				which_column_list and group_size_list
		"""
		sys.stderr.write("Outputting...")
		writer = csv.writer(outf, delimiter='\t')
		header_row = ['gene_no','go_no','is_correct', 'is_correct_ratio', 'p_value', 'recurrence', \
			'connectivity', 'cluster_size', 'unknown', 'acc1', 'acc2', 'no', 'mcl_no']
		if len(self.which_column_list)>1 and len(self.group_size_list)>1:
			writer.writerow(['parameter1', 'parameter2', 'acc1', 'acc2', 'no', 'mcl_no'])
		else:
			writer.writerow(['parameter', 'acc1', 'acc2', 'no', 'mcl_no'])
		for (prediction_space,unit) in prediction_space2attr.iteritems():
			unit_array = array(unit)

			row = list(prediction_space)
			row.append(get_average(unit_array,2))	#average is_correct
			row.append(get_average(unit_array, 3))	#average is_correct_ratio
			#row.append(len(unit_array))	#no of predictions
			if self.report:
				row.append(repr(Set(unit_array[:,-2])))
				row.append(repr(Set(unit_array[:,-1])))
			else:
				row.append(len(Set(unit_array[:,-2])))
				row.append(len(Set(unit_array[:,-1])))	#no of clusters
			writer.writerow(row)
		sys.stderr.write("Done.\n")
		del writer

	def run(self):
		"""
		03-30-05
		
		06-30-05
			more complex data grouping via which_column_list and group_size_list
			if both lists are of length 2, 2-level grouping.
			
		--db_connect()
		--get_go_no2depth()
		--data_fetch()
		--group_data()
		if self.stat_table_fname:
			--prediction_space_output()
		"""
		self.init()
		(conn, curs) = db_connect(self.hostname, self.dbname, self.schema)
		from codense.common import  get_go_no2depth
		self.go_no2depth = get_go_no2depth(curs)
		
		self.data_fetch(curs, self.table, self.mcl_table, self.gene_table)
		local_prediction_space2attr = self.group_data(self.prediction_data,key_column=self.which_column_list[0], group_size=self.group_size_list[0])
		for key, unit in local_prediction_space2attr.iteritems():
			if len(self.which_column_list)>1 and len(self.group_size_list)>1:
				local_prediction_space2attr_2 = self.group_data(unit, key_column=self.which_column_list[1], group_size=self.group_size_list[1])
				for key2, unit2 in local_prediction_space2attr_2.iteritems():
					self.prediction_space2attr[(key,key2)] = unit2
			else:
				self.prediction_space2attr[(key,)] = unit
		stat_table_f = open(self.stat_table_fname, 'w')
		self.prediction_space_output(stat_table_f, self.prediction_space2attr)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema=", "netmine_fname=", "table=", "mcl_table=", \
		"judger_type=", "report", "gene_table=", "which_column=", "group_size=", "debug"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:f:t:m:j:rg:w:s:b", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	netmine_fname = ''
	table = 'cluster_stat'
	mcl_table = 'mcl_result'
	judger_type = 0
	report = 0
	gene_table = 'p_gene'
	debug = 0
	which_column = '5'
	group_size = '10'
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
		elif opt in ("-f", "--netmine_fname"):
			netmine_fname = arg
		elif opt in ("-t", "--table"):
			table = arg
		elif opt in ("-m", "--mcl_table"):
			mcl_table = arg
		elif opt in ("-j", "--judger_type"):
			judger_type = int(arg)
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-g", "--gene_table"):
			gene_table = arg
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-w", "--which_column"):
			which_column = arg
		elif opt in ("-s", "--group_size"):
			group_size = arg
		
	if len(args) == 1:
		stat_table_fname = args[0]
	
	if netmine_fname:
		table = 'splat_%s'%netmine_fname
		mcl_table = 'mcl_%s'%netmine_fname
		gene_table = 'p_gene_%s_e5'%netmine_fname
	
	if schema and table and mcl_table and gene_table and stat_table_fname:
		instance = p_gene_factor(hostname, dbname, schema, table, mcl_table,\
			report, judger_type, gene_table, stat_table_fname, which_column, \
			group_size, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
