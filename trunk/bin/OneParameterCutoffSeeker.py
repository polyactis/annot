#!/usr/bin/env python
"""
Usage: OneParameterCutoffSeeker.py -k SCHEMA -t P_GENE_TABLE -l LM_TABLE [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ...,	the p_gene table
	-l ...,	the lm_table to store the linear_model results, needed if needcommit
	-a ...,	0.5(default)
	-j ...,	how to judge predicted functions, 0(default), 1, 2
	-w ...,	a bit_string controlling which parameter to use, 00001(default, edge_gradient)
	-c,	commit this database transaction
	-r,	report flag
	-u,	debug flag
	-h, --help              show this help

Examples:
	
Description:
	Calculate the parameter cutoff corresponding to an accuracy_cut_off.
	Which parameter based on bit_string(e.g. 001, 1, 01, 0001, 00010):
		1: 'p_value_cut_off'; 2: 'recurrence_cut_off'; 3: 'connectivity_cut_off';
		4: 'cluster_size_cut_off'; 5: 'edge_gradient'; 6:'rpart_prob'
	In MpiStatCluster.py's resulting p_gene table, p_value_cut_off=edge_gradient.
	10-28-05 memory threshold is 1.5e7
"""

import sys, os, getopt
sys.path += [os.path.expanduser('~/script/annot/bin')]
from codense.common import db_connect
from p_gene_lm import p_gene_lm
from heapq import heappush, heappop

class OneParameterCutoffSeeker:
	def __init__(self, hostname=None, dbname=None, schema=None, p_gene_table=None, \
		lm_table=None, accuracy_cut_off=0, judger_type=0, which=0, commit=0, report=0, debug=0):
		"""
		11-09-05 add the sixth parameter: rpart_prob(which is vertex_gradient in p_gene_table)
		"""
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
			4: 'edge_gradient',
			5: 'vertex_gradient'}
	
	def get_prediction_step(self, curs, table, is_correct_dict, judger_type):
		"""
		10-27-05 get a step to control the memory, from p_gene_lm.py
			diff: 1e7 becomes 2.5e7 and the where condition
		10-28-05 downsize 2.5e7 to 1.5e7
		"""
		sys.stderr.write("Getting step...\n")
		curs.execute("select count(*) from %s where %s != -1"%(table, is_correct_dict[judger_type]))
		rows = curs.fetchall()
		no_of_entries = rows[0][0]
		step = int(no_of_entries/1.5e7) + 1
		sys.stderr.write("Got step=%s.\n"%step)
		return step
	
	def get_prediction_heap(self, curs, p_gene_table, is_correct_dict, judger_type, which_dict, which, step):
		"""
		10-27-05 get the prediction into a heap sorted by the param
		10-27-05 add step to control memory
		"""
		sys.stderr.write("Getting prediction_heap...\n")
		prediction_heap = []
		curs.execute("DECLARE crs CURSOR FOR select p.%s, p.%s, p.gene_no, p.go_no from %s p"\
			%(which_dict[which], is_correct_dict[judger_type], p_gene_table))
		curs.execute("fetch 10000 from crs")
		rows = curs.fetchall()
		counter = 0
		while rows:
			for row in rows:
				param, is_correct, gene_no, go_no = row
				if is_correct!=-1:	#unknown prediction is not needed
					if (counter%step) == 0:	#10-27-05 step control
						param = -param	#10-27-05 reverse the sign to use the min-heap as max-heap
						heappush(prediction_heap, [param, is_correct, gene_no, go_no])
					counter += 1
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, counter))
			curs.execute("fetch 10000 from crs")
			rows = curs.fetchall()
		sys.stderr.write("Got prediction_heap.\n")
		return prediction_heap
	
	def get_sorted_param_acc_list(self, prediction_heap):
		"""
		10-27-05 Calculate the accuracy for each param from highest to lowest
		"""
		sys.stderr.write("Getting sorted_param_acc_list...\n")
		sorted_param_acc_list = []
		prediction_pair2highest_param = {}
		if len(prediction_heap)>0:
			prev_param, is_correct, gene_no, go_no = heappop(prediction_heap)
			prev_param = -prev_param	#reverse the sign
			prediction_pair2highest_param[(gene_no,go_no)] = prev_param
			if is_correct==1:
				no_of_correct = 1
				no_of_wrong = 0
			elif is_correct==0:
				no_of_correct = 0
				no_of_wrong = 1
		else:
			sys.stderr.write("No prediction.\n")
			sys.exit(2)
		while len(prediction_heap)>0:
			param, is_correct, gene_no, go_no = heappop(prediction_heap)
			param = -param	#reverse the sign
			if param != prev_param:	#new param, calculate the acc for the prev_param
				acc = float(no_of_correct)/(no_of_correct+no_of_wrong)
				sorted_param_acc_list.append([prev_param, no_of_correct, no_of_wrong, acc])
				if self.debug:
					print "prev_param",prev_param,"no_of_correct",no_of_correct,"no_of_wrong",no_of_wrong,"acc",acc
					raw_input("Continue?(Y/n)")
				prev_param = param
			#following prediction_pair judgement must be after param!=prev_param judgement
			prediction_pair = (gene_no, go_no)
			if prediction_pair not in prediction_pair2highest_param:	#judge whether the prediction_pair redundant or not
				prediction_pair2highest_param[prediction_pair] = param
				if is_correct==1:
					no_of_correct += 1
				elif is_correct==0:
					no_of_wrong += 1
		#handle the last param
		acc = float(no_of_correct)/(no_of_correct+no_of_wrong)
		sorted_param_acc_list.append([prev_param, no_of_correct, no_of_wrong, acc])
		
		sys.stderr.write("Got sorted_param_acc_list.\n")
		return sorted_param_acc_list
	
	def get_cutoff(self, sorted_param_acc_list, accuracy_cut_off):
		"""
		10-27-05 From lowest to highest, seek the param giving the accuracy>=accuracy_cut_off
		"""
		sys.stderr.write("Getting cutoff...\n")
		no_of_diff_params = len(sorted_param_acc_list)
		for i in range(len(sorted_param_acc_list)):
			if sorted_param_acc_list[no_of_diff_params-1-i][-1] >= accuracy_cut_off:
				break
		if sorted_param_acc_list[no_of_diff_params-1-i][-1] >= accuracy_cut_off:
			sys.stderr.write("Got cutoff.\n")
			return sorted_param_acc_list[no_of_diff_params-1-i]
		else:
			sys.stderr.write("No param's value meeting the accuracy.\n")
			return None
	
	def run(self):
		"""
		10-27-05
			
			--db_connect()
			--get_prediction_step()
			--get_prediction_heap()
			--get_sorted_param_acc_list()
			--get_cutoff()
			--lm_table_create()
			--submit()
		"""
		p_gene_lm_instance = p_gene_lm()
		(conn, curs) = db_connect(self.hostname, self.dbname, self.schema)
		step = self.get_prediction_step(curs, self.p_gene_table, p_gene_lm_instance.is_correct_dict, \
			self.judger_type)
		prediction_heap = self.get_prediction_heap(curs, self.p_gene_table, p_gene_lm_instance.is_correct_dict, \
			self.judger_type, self.which_dict, self.which, step)
		sorted_param_acc_list = self.get_sorted_param_acc_list(prediction_heap)
		del prediction_heap	#10-27-05 release memory
		cutoff_row = self.get_cutoff(sorted_param_acc_list, self.accuracy_cut_off)
		del sorted_param_acc_list	#10-27-05 release memory
		print "cutoff_row",cutoff_row
		if self.commit and cutoff_row and self.lm_table:	#cutoff_row is not None
			p_gene_lm_instance.lm_table_create(curs, self.lm_table)
			go_no2lm_results = {}
			go_no2lm_results[-1] = [[0]*7, [1]*7, cutoff_row[0]]	#11-09-05 extend the list
			go_no2lm_results[-1][0][which+1] = 1	#the coeffcient for "which" param is 1, others are 0
			p_gene_lm_instance.submit(curs, self.lm_table, go_no2lm_results)
			curs.execute("end")

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:l:a:j:w:cru", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	p_gene_table = None
	lm_table = None
	accuracy_cut_off = 0.5
	judger_type = 0
	which = 0
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
		elif opt in ("-t",):
			p_gene_table = arg
		elif opt in ("-l",):
			lm_table = arg
		elif opt in ("-a",):
			accuracy_cut_off = float(arg)
		elif opt in ("-j",):
			judger_type = int(arg)
		elif opt in ("-w",):
			for i in range(len(arg)):	#10-27-05 find the first index which is not 0
				if arg[i] == '1':
					break
			which = i
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-r",):
			report = 1
		elif opt in ("-u",):
			debug = 1
	
	if schema and p_gene_table:
		instance = OneParameterCutoffSeeker(hostname, dbname, schema, p_gene_table, \
			lm_table, accuracy_cut_off, judger_type, which, commit, report, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
