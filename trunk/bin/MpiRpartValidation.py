#!/usr/bin/env python
"""
Usage: MpiRpartValidation.py -k -i -j

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-i ...,	fname of schema setting
	-j ...,	output_file
	-f ...,	filter type
	-y ...,	is_correct type (2 lca, default)
	-p ...,	rpart cp value list (0.01, default)
	-l ...,	loss matrix list, 0,1,1,0 (default) i.e. 0,1,1,0=0,2,1,0
	-o ...,	prior prob list, 0.5 (default)
	-s ...,	percentage of data selected to do training, 0.8(default)
	-x ...,	no_of_validations, 10(default)
	-g, 	calculate the hypergeometric p-value to replace p_value_cut_off(gradient)
	-b,	enable debug flag
	-c,	commit the database transaction(IGNORE)
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	mpirun -np 20 -machinefile ~/hostfile /usr/bin/mpipython  ~/script/annot/bin/MpiRpartValidation.py
		-k hs_fim_40 -i hs_fim_40m4x40rec0_8 -j output_file -c -r

Description:
	Program to do rpart validation. Inherit rpart_prediction.py
"""

import sys, os, getopt, csv, cPickle
sys.path += [os.path.expanduser('~/script/annot/bin')]
from Scientific import MPI
from codense.common import db_connect, form_schema_tables, \
	get_go_no2gene_no_set, get_no_of_total_genes, output_node, \
	computing_node, mpi_synchronize
from rpart_prediction import rpart_prediction
from rpy import r

class MpiRpartValidation(rpart_prediction):
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, fname=None, output_file=None, \
		filter_type=1, is_correct_type=2, rpart_cp_ls=[0.01], loss_matrix_ls=[[0,1,1,0]], prior_prob_ls=[0.5], \
		training_perc=0.8, no_of_validations=10, need_cal_hg_p_value=0, debug=0, commit=0, report=0):
		"""
		11-19-05 add no_of_validations
		"""
		rpart_prediction.__init__(self, hostname, dbname, schema, fname, fname, \
			filter_type, is_correct_type, 0.01, [0,1,1,0], None, training_perc, \
			need_cal_hg_p_value, debug, commit, report)
		self.fname = fname
		self.output_file = output_file
		self.rpart_cp_ls = rpart_cp_ls
		self.loss_matrix_ls = loss_matrix_ls
		self.prior_prob_ls  = prior_prob_ls
		self.no_of_validations = int(no_of_validations)
	
	def form_setting_ls(self, rpart_cp_ls, loss_matrix_ls, prior_prob_ls):
		setting_ls = []
		for rpart_cp in rpart_cp_ls:
			for loss_matrix in loss_matrix_ls:
				for prior_prob in prior_prob_ls:
					setting_ls.append([rpart_cp, loss_matrix, prior_prob])
		return setting_ls
	
	def get_data(self, curs, fname, filter_type, is_correct_type, need_cal_hg_p_value):
		"""
		11-19-05
			data_fetch() of rpart_prediction.py changed
			return unknown_data
		"""
		schema_instance = form_schema_tables(fname)
		
		no_of_total_genes = get_no_of_total_genes(curs)
		go_no2gene_no_set = get_go_no2gene_no_set(curs)
		unknown_prediction_ls, known_prediction_ls, unknown_data, known_data = self.data_fetch(curs, schema_instance, \
			filter_type, is_correct_type, no_of_total_genes, go_no2gene_no_set, need_cal_hg_p_value)
		del unknown_prediction_ls, known_prediction_ls
		return unknown_data, known_data
	
	
	def input_node(self, communicator, setting_ls, free_computing_nodes, report=0):
		node_rank = communicator.rank
		sys.stderr.write("Input node(%s) working...\n"%node_rank)
		counter = 0
		while setting_ls:
			setting = setting_ls.pop(0)
			communicator.send("1", communicator.size-1, 1)	#WATCH: tag is 1, to the output_node.
			free_computing_node, source, tag = communicator.receiveString(communicator.size-1, 2)
				#WATCH: tag is 2, from the output_node
			setting_pickle = cPickle.dumps(setting, -1)
			communicator.send(setting_pickle, int(free_computing_node), 0)	#WATCH: int()
			if report:
				sys.stderr.write("block %s sent to %s.\n"%(counter, free_computing_node))
			counter += 1
		#tell computing_node to exit the loop
		for node in free_computing_nodes:	#send it to the computing_node
			communicator.send("-1", node, 0)
		sys.stderr.write("Input node(%s) done\n"%(node_rank))
	
	def computing_handler(self, communicator, data, parameter_list):
		"""
		11-19-05
			add no_of_validations
			add unknown_data and do rpart_fit_and_predict on known_data and unknown_data
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		unknown_data, known_data, training_perc, no_of_validations = parameter_list
		rpart_cp, loss_matrix, prior_prob = data
		testing_acc_ls, training_acc_ls = self.rpart_validation(known_data, training_perc, rpart_cp, \
			loss_matrix, prior_prob, no_of_validations)
		
		unknown_pred, known_pred = self.rpart_fit_and_predict(unknown_data, known_data, rpart_cp, loss_matrix, prior_prob)
		result = [[rpart_cp, loss_matrix, prior_prob],\
			self.summary_stat_of_accuracy_ls(testing_acc_ls),\
			self.summary_stat_of_accuracy_ls(training_acc_ls),\
			self.cal_accuracy(unknown_data, unknown_pred),\
			self.cal_accuracy(known_data, known_pred)]
		sys.stderr.write("Node no.%s done.\n"%node_rank)
		return result
	
	def summary_stat_of_accuracy_ls(self, acc_ls):
		"""
		11-17-05
			get average, standard deviation from acc_ls
		"""
		accuracy_ls, no_of_predictions_ls, no_of_genes_ls = [], [], []
		for row in acc_ls:
			accuracy_ls.append(row[0])
			no_of_predictions_ls.append(row[2])
			no_of_genes_ls.append(row[3])
		no_of_samples = float(len(acc_ls))
		accuracy_avg = sum(accuracy_ls)/no_of_samples
		accuracy_std = r.sqrt(r.var(accuracy_ls)/no_of_samples)
		no_of_predictions_avg = sum(no_of_predictions_ls)/no_of_samples
		no_of_predictions_std = r.sqrt(r.var(no_of_predictions_ls)/no_of_samples)
		no_of_genes_avg = sum(no_of_genes_ls)/no_of_samples
		no_of_genes_std = r.sqrt(r.var(no_of_genes_ls)/no_of_samples)
		return [accuracy_avg,accuracy_std,no_of_predictions_avg,no_of_predictions_std,no_of_genes_avg,no_of_genes_std]
	
	def output_handler(self, communicator, parameter_list, data):
		"""
		11-19-05 output the result from unknown and known
		"""
		writer = parameter_list[0]
		data = cPickle.loads(data)
		writer.writerow(data[0] + ['testing']+data[1])	#testing result
		writer.writerow(data[0] + ['training']+data[2])	#training result
		writer.writerow(data[0] + ['unknown']+data[3])
		writer.writerow(data[0] + ['known']+data[4])
	
	def run(self):
		"""
		11-16-05
			
			--computing_handler()
				--is_site_confirmed()
					--get_no_of_mismatches_allowed()
					--get_no_of_mismatches_for_consensus()
						--is_good_consensus()
					--get_no_of_mismatches_for_site()
		"""
		communicator = MPI.world.duplicate()
		node_rank = communicator.rank	
		free_computing_nodes = range(1,communicator.size-1)	#exclude the last node
		if node_rank == 0:
			(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
			unknown_data, known_data = self.get_data(curs, self.fname, self.filter_type, self.is_correct_type, self.need_cal_hg_p_value)
			known_data_pickle = cPickle.dumps(known_data, -1)
			for node in free_computing_nodes:	#send it to the computing_node
				communicator.send(known_data_pickle, node, 0)
			unknown_data_pickle = cPickle.dumps(unknown_data, -1)
			for node in free_computing_nodes:	#send it to the computing_node
				communicator.send(unknown_data_pickle, node, 0)
		elif node_rank in free_computing_nodes:
			data, source, tag = communicator.receiveString(0, 0)
			known_data = cPickle.loads(data)	#take the data
			data, source, tag = communicator.receiveString(0, 0)
			unknown_data = cPickle.loads(data)	#take the data
		elif node_rank==communicator.size-1:
			writer = csv.writer(open(self.output_file, 'w'), delimiter='\t')
			#write down the header
			writer.writerow(['rpart_cp', 'loss_matrix', 'prior_prob', 'type', 'accuracy_avg','accuracy_std', 'no_of_predictions_avg',\
				'no_of_predictions_std', 'no_of_genes_avg', 'no_of_genes_std'])
			
		mpi_synchronize(communicator)
		if node_rank == 0:
			setting_ls = self.form_setting_ls(self.rpart_cp_ls, self.loss_matrix_ls, self.prior_prob_ls)
			self.input_node(communicator, setting_ls, free_computing_nodes, self.report)
		elif node_rank in free_computing_nodes:
			parameter_list = [unknown_data, known_data, self.training_perc, self.no_of_validations]
			computing_node(communicator, parameter_list, self.computing_handler, report=self.report)
		elif node_rank==communicator.size-1:
			parameter_list = [writer]
			output_node(communicator, free_computing_nodes, parameter_list, self.output_handler, self.report)
			del writer

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:j:f:y:p:l:o:s:x:gbcr", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	fname = None
	output_file = None
	filter_type = 1
	is_correct_type = 2
	rpart_cp_ls = [0.01]
	loss_matrix_ls = [[0,1,1,0]]
	prior_prob_ls = [0.5]
	training_perc = 0.8
	no_of_validations = 10
	need_cal_hg_p_value = 0
	debug = 0
	commit = 0
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
			fname = arg
		elif opt in ("-j"):
			output_file = arg
		elif opt in ("-f"):
			filter_type = int(arg)
		elif opt in ("-y"):
			is_correct_type = int(arg)
		elif opt in ("-p"):
			rpart_cp_ls = arg.split(',')
			rpart_cp_ls = map(float, rpart_cp_ls)
		elif opt in ("-l"):
			loss_matrix_ls = arg.split('=')
			for i in range(len(loss_matrix_ls)):
				loss_matrix_ls[i] = map(float, loss_matrix_ls[i].split(','))
		elif opt in ("-o"):
			prior_prob_ls = map(float, arg.split(','))
		elif opt in ("-s"):
			training_perc = float(arg)
		elif opt in ("-x"):
			no_of_validations = int(arg)
		elif opt in ("-g"):
			need_cal_hg_p_value = 1
		elif opt in ("-b"):
			debug = 1
		elif opt in ("-c"):
			commit = 1
		elif opt in ("-r"):
			report = 1
	if schema and fname and output_file:
		instance = MpiRpartValidation(hostname, dbname, schema, fname, output_file, \
			filter_type, is_correct_type, rpart_cp_ls, loss_matrix_ls, prior_prob_ls, training_perc, \
			no_of_validations, need_cal_hg_p_value, debug, commit, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
