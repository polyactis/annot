#!/usr/bin/env python
"""
Usage: PredictionFilterByClusterSize.py -k SCHEMA -s SPLAT_TABLE -m MCL_TABLE
	-p P_GENE_TABLE -t SPLAT_VIEW -n MCL_VIEW -q P_GENE_VIEW [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-i ...,	the old input_fname
	-j ...,	the new input_fname
	-m ...,	max cluster size, 100000 (default, no cut)
	-u ...,	unknown_gene_ratio, 1(default, no cut)
	-p ...,	p_value_cut_off, 1(default, no cut)
	-y ...,	is_correct type (2 lca, default)
	-a ...,	cluster accuracy cutoff, (0 default, no cut)
	-c,	commit the database transaction
	-b,	debug version.
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	PredictionFilterByClusterSize.py -k hs_fim_40 -i
		hs_fim_40m4x40 -j hs_fim_40m4x40ms12000m -c -m 12000
	
Description:
	Program to filter clusters by max_size. Database transaction
	is done by creating views.
	10-06-05, not just size cut. 

"""

import sys, os, getopt
sys.path += [os.path.expanduser('~/script/annot/bin')]
from codense.common import db_connect, form_schema_tables
from gene_stat import gene_stat
from sets import Set

class prediction_attributes:
	def __init__(self, row):
		self.p_gene_id = row[0]
		self.gene_no = row[1]
		self.go_no = row[2]
		self.is_correct = row[3]
		self.is_correct_l1 = row[4]
		self.is_correct_lca = row[5]
		self.avg_p_value = row[6]
		self.no_of_clusters = row[7]
		self.cluster_array = row[8]
		self.p_value_cut_off = row[9]
		self.recurrence_cut_off = row[10]
		self.connectivity_cut_off = row[11]
		self.cluster_size_cut_off = row[12]
		self.unknown_cut_off = row[13]
		self.depth_cut_off = row[14]
		self.mcl_id = row[15]
		self.lca_list = row[16]
		
		self.is_correct_dict = {0:self.is_correct, 1:self.is_correct_l1, 2:self.is_correct_lca}

class PredictionFilterByClusterSize:
	"""
	09-26-05
		upgrade the function-form to class, 
		p_gene_view is not a view anymore, a real table derived from p_gene_table
			because runtime select cluster_size_cut_off<=max_size blows the memory
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, input_fname=None,\
		jnput_fname=None, max_size=100000, unknown_gene_ratio=1, p_value_cut_off=1,\
		is_correct_type=2, acc_cut_off=0, commit=0, debug=0, report=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_fname = input_fname
		self.jnput_fname = jnput_fname
		self.max_size = int(max_size)
		self.unknown_gene_ratio = float(unknown_gene_ratio)
		self.p_value_cut_off = float(p_value_cut_off)
		self.is_correct_type = int(is_correct_type)
		self.acc_cut_off = float(acc_cut_off)
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
		
	
	def view_from_table(self, curs, table, view):
		curs.execute("CREATE OR REPLACE VIEW %s AS SELECT * FROM %s"%(view, table))
	
	def is_good_prediction(self, prediction_attr_instance, mcl_id2accuracy=None):
		"""
		10-05-05
			judge whether this prediction should be taken into consideration
		10-06-05
			add judgement by acc_cut_off
		"""
		if prediction_attr_instance.p_value_cut_off>self.p_value_cut_off or \
			prediction_attr_instance.cluster_size_cut_off>self.max_size or \
			prediction_attr_instance.unknown_cut_off>self.unknown_gene_ratio:
			return 0
		else:
			if mcl_id2accuracy:
				if mcl_id2accuracy[prediction_attr_instance.mcl_id]>=self.acc_cut_off:
					return 1
				else:
					return 0
			else:
				return 1

	
	def get_mcl_id2accuracy(self, curs, p_gene_table, crs_sentence, is_correct_type):
		"""
		10-05-05
			
		"""
		sys.stderr.write("Getting mcl_id2accuracy...\n")
		mcl_id2accuracy = {}
		curs.execute("%s"%crs_sentence)
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		counter = 0
		while rows:
			for row in rows:
				prediction_attr_instance = prediction_attributes(row)
				if self.is_good_prediction(prediction_attr_instance):
					if prediction_attr_instance.mcl_id not in mcl_id2accuracy:
						mcl_id2accuracy[prediction_attr_instance.mcl_id]  = []
					mcl_id2accuracy[prediction_attr_instance.mcl_id].append(prediction_attr_instance.is_correct_dict[is_correct_type])
				
				counter += 1
			if self.report:
				sys.stderr.write("%s%s"%('\x08'*20, counter))
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		curs.execute("close crs")
		for mcl_id, is_correct_ls in mcl_id2accuracy.iteritems():
			accuracy = sum(is_correct_ls)/float(len(is_correct_ls))
			mcl_id2accuracy[mcl_id] = accuracy
		sys.stderr.write(" %s clusters. Done.\n"%len(mcl_id2accuracy))
		return mcl_id2accuracy
	
	def submit_to_p_gene_table(self, curs, p_gene_table, p_attr_instance):
		"""
		10-05-05
		"""
		if p_attr_instance.lca_list:
			curs.execute("insert into %s(p_gene_id, gene_no, go_no, is_correct, is_correct_l1, \
			is_correct_lca, avg_p_value, no_of_clusters, cluster_array, p_value_cut_off, recurrence_cut_off, \
			connectivity_cut_off, cluster_size_cut_off, unknown_cut_off, depth_cut_off, mcl_id, lca_list)\
			values(%s, %s, %s, %s, %s,\
			%s, %s, %s, '%s', %s, %s,\
			%s, %s, %s, %s, %s, '%s')"%\
			(p_gene_table, \
			p_attr_instance.p_gene_id, p_attr_instance.gene_no, p_attr_instance.go_no, p_attr_instance.is_correct, p_attr_instance.is_correct_l1,\
			p_attr_instance.is_correct_lca, p_attr_instance.avg_p_value, p_attr_instance.no_of_clusters, p_attr_instance.cluster_array, p_attr_instance.p_value_cut_off, p_attr_instance.recurrence_cut_off,\
			p_attr_instance.connectivity_cut_off, p_attr_instance.cluster_size_cut_off, p_attr_instance.unknown_cut_off, p_attr_instance.depth_cut_off, p_attr_instance.mcl_id, p_attr_instance.lca_list))
		else:
			curs.execute("insert into %s(p_gene_id, gene_no, go_no, is_correct, is_correct_l1, \
			is_correct_lca, avg_p_value, no_of_clusters, cluster_array, p_value_cut_off, recurrence_cut_off, \
			connectivity_cut_off, cluster_size_cut_off, unknown_cut_off, depth_cut_off, mcl_id)\
			values(%s, %s, %s, %s, %s,\
			%s, %s, %s, '%s', %s, %s,\
			%s, %s, %s, %s, %s)"%\
			(p_gene_table, \
			p_attr_instance.p_gene_id, p_attr_instance.gene_no, p_attr_instance.go_no, p_attr_instance.is_correct, p_attr_instance.is_correct_l1,\
			p_attr_instance.is_correct_lca, p_attr_instance.avg_p_value, p_attr_instance.no_of_clusters, p_attr_instance.cluster_array, p_attr_instance.p_value_cut_off, p_attr_instance.recurrence_cut_off,\
			p_attr_instance.connectivity_cut_off, p_attr_instance.cluster_size_cut_off, p_attr_instance.unknown_cut_off, p_attr_instance.depth_cut_off, p_attr_instance.mcl_id))
		
	def setup_new_p_gene_table(self, curs, p_gene_table, p_gene_view, crs_sentence, mcl_id2accuracy=None, acc_cut_off=0):
		"""
		10-05-05
		10-06-05	(use acc_cut_off, not keep the most accurate one
			simple, just use is_good_prediction() to judge if it's good or not
		"""
		sys.stderr.write("Starting to setup new p_gene_table...\n")
		curs.execute("%s"%(crs_sentence))
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		counter = 0
		no_of_good_predictions = 0
		gene_no2go2prediction_dict = {}
		while rows:
			for row in rows:
				p_attr_instance = prediction_attributes(row)
				if self.is_good_prediction(p_attr_instance, mcl_id2accuracy):
					no_of_good_predictions += 1
					"""
					if acc_cut_off:
						if p_attr_instance.gene_no not in gene_no2go2prediction_dict:
							gene_no2go2prediction_dict[p_attr_instance.gene_no] = {}
						unit = gene_no2go2prediction_dict[p_attr_instance.gene_no]
						go_no = p_attr_instance.go_no
						if go_no not in unit:
							unit[go_no] = []
						current_mcl_id_accuracy = mcl_id2accuracy[p_attr_instance.mcl_id]
						if current_mcl_id_accuracy>=acc_cut_off:	#current cluster's accuracy >= cutoff
							unit[go_no].append(p_attr_instance.p_gene_id)
					else:
					"""
					self.submit_to_p_gene_table(curs, p_gene_view, p_attr_instance)
				counter += 1
			if self.report:
				sys.stderr.write("%s%s:%s"%('\x08'*20, counter, no_of_good_predictions))
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		curs.execute("close crs")
		sys.stderr.write(" %s good predictions. Done.\n"%(no_of_good_predictions))
		return gene_no2go2prediction_dict
	
	def distinct_predictions_set(self, gene_no2go2prediction_dict):
		"""
		10-05-05
		10-06-05
			the value of go2prediction_dict is just a list of p_gene_ids
		"""
		sys.stderr.write("Getting distinct predictions...\n")
		p_gene_id_set = Set()
		for gene_no, go2prediction_dict in gene_no2go2prediction_dict.iteritems():
			for go_no, ls in go2prediction_dict.iteritems():
				p_gene_id_set |= Set(ls)
		sys.stderr.write("Done.\n")
		return p_gene_id_set
	
	def submit_distinct_predictions(self, curs, new_p_gene_table, crs_sentence, p_gene_id_set):
		"""
		10-05-05
		"""
		sys.stderr.write("Submitting distinct predictions...\n")
		curs.execute("%s"%(crs_sentence))
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		no_of_distinct_p_gene_ids = 0
		counter = 0
		while rows:
			for row in rows:
				p_attr_instance = prediction_attributes(row)
				if p_attr_instance.p_gene_id in p_gene_id_set:
					no_of_distinct_p_gene_ids += 1
					self.submit_to_p_gene_table(curs, new_p_gene_table, p_attr_instance)
				counter += 1
			if self.report:
				sys.stderr.write("%s%s:%s"%('\x08'*20, counter, no_of_distinct_p_gene_ids))
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		curs.execute("close crs")
		sys.stderr.write(" %s distinct predictions. Done.\n"%no_of_distinct_p_gene_ids)
		
	def run(self):
		"""
		10-05-05
			
			--db_connect()
			--form_schema_tables()
			--form_schema_tables()
			--view_from_table()
			--view_from_table()
			
			--createGeneTable()
			if acc_cut_off:
				--get_mcl_id2accuracy()
					--is_good_prediction()
			
			--setup_new_p_gene_table()
				--is_good_prediction()
				if acc_cut_off
					--submit_to_p_gene_table()
		10-06-05(below is not called anymore)
			if acc_cut_off
				--distinct_predictions_set()
				--submit_distinct_predictions()
		"""
		(conn, curs) =  db_connect(hostname, dbname, schema)
		old_schema_instance = form_schema_tables(self.input_fname)
		new_schema_instance = form_schema_tables(self.jnput_fname)
		self.view_from_table(curs, old_schema_instance.splat_table, new_schema_instance.splat_table)
		self.view_from_table(curs, old_schema_instance.mcl_table, new_schema_instance.mcl_table)
		
		crs_sentence = 'DECLARE crs CURSOR FOR SELECT p_gene_id, gene_no, go_no, is_correct, is_correct_l1, \
			is_correct_lca, avg_p_value, no_of_clusters, cluster_array, p_value_cut_off, recurrence_cut_off, \
			connectivity_cut_off, cluster_size_cut_off, unknown_cut_off, depth_cut_off, mcl_id, lca_list  \
			from %s'%old_schema_instance.p_gene_table
		
		gene_stat_instance = gene_stat()
		gene_stat_instance.createGeneTable(curs, new_schema_instance.p_gene_table)
		if self.acc_cut_off:
			mcl_id2accuracy = self.get_mcl_id2accuracy(curs, old_schema_instance.p_gene_table, crs_sentence, self.is_correct_type)
		else:
			mcl_id2accuracy = None
		gene_no2go2prediction_dict = self.setup_new_p_gene_table(curs, old_schema_instance.p_gene_table, \
			new_schema_instance.p_gene_table, crs_sentence, mcl_id2accuracy, self.acc_cut_off)
		if commit:
			curs.execute("End")

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:j:m:u:p:y:a:cbr", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	input_fname = None
	jnput_fname = None
	max_size = 100000
	unknown_gene_ratio = 1
	p_value_cut_off = 1
	is_correct_type = 2
	acc_cut_off = 0
	commit = 0
	debug = 0
	report = 0
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
		elif opt in ("-i"):
			input_fname = arg
		elif opt in ("-j"):
			jnput_fname = arg
		elif opt in ("-m"):
			max_size = int(arg)
		elif opt in ("-u"):
			unknown_gene_ratio = float(arg)
		elif opt in ("-p"):
			p_value_cut_off = float(arg)
		elif opt in ("-y"):
			is_correct_type = int(arg)
		elif opt in ("-a"):
			acc_cut_off = float(arg)
		elif opt in ("-c"):
			commit = 1
		elif opt in ("-b"):
			debug = 1
		elif opt in ("-r"):
			report = 1
	if schema and input_fname and jnput_fname:
		instance = PredictionFilterByClusterSize(hostname, dbname, schema, input_fname, jnput_fname, max_size, \
			unknown_gene_ratio, p_value_cut_off, is_correct_type, acc_cut_off, commit, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
