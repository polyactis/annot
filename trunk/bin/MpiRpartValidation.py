#!/usr/bin/env mpipython
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
	-o ...,	prior prob list, 0.5 (default), 0 means proportional to the observed data.
	-s ...,	percentage of data selected to do training, 0.8(default)
	-x ...,	no_of_validations, 10(default)
	-t ...,	type, 1(default, rpart), 2(randomForest)
	-m ...,	mty list(see rpart_prediction.py's doc)
			0(default, automatically choosen by randomForest)
	-b ...,	bit_string, '1111111'(default)
	-g, 	calculate the hypergeometric p-value to replace p_value_cut_off(gradient)
	-u,	enable debug flag
	-c,	commit the database transaction(IGNORE)
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	mpirun -np 20 -machinefile ~/hostfile /usr/bin/mpipython  ~/script/annot/bin/MpiRpartValidation.py
		-k hs_fim_40 -i hs_fim_40m4x40rec0_8 -j output_file -c -r

Description:
	Program to do rpart validation. Inherit rpart_prediction.py
	For bit_string's meaning see rpart_prediction.py's doc.
	
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import getopt, csv, cPickle, random
from Scientific import MPI
from codense.common import db_connect, form_schema_tables, \
	get_go_no2gene_no_set, get_no_of_total_genes, output_node, \
	computing_node, mpi_synchronize
from rpart_prediction import rpart_prediction
import rpy
from rpy import r
from sets import Set
if sys.version_info[:2] < (2, 3):       #python2.2 or lower needs some extra
        from python2_3 import *

class MpiRpartValidation(rpart_prediction):
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, fname=None, output_file=None, \
		filter_type=1, is_correct_type=2, rpart_cp_ls=[0.01], loss_matrix_ls=[[0,1,1,0]], prior_prob_ls=[0.5], \
		training_perc=0.8, no_of_validations=10, type=1, mty_ls=[0], bit_string='11111', need_cal_hg_p_value=0, debug=0, commit=0, report=0):
		"""
		11-19-05 add no_of_validations
		03-17-06
			add type, mty_ls, bit_string
		"""
		rpart_prediction.__init__(self, hostname, dbname, schema, fname, fname, \
			filter_type, is_correct_type, 0.01, [0,1,1,0], None, training_perc, \
			type, 2, bit_string, need_cal_hg_p_value, debug, commit, report)
		self.fname = fname
		self.output_file = output_file
		self.rpart_cp_ls = rpart_cp_ls
		self.loss_matrix_ls = loss_matrix_ls
		self.prior_prob_ls  = prior_prob_ls
		self.type = int(type)
		self.mty_ls = mty_ls
		self.bit_string = bit_string
		self.no_of_validations = int(no_of_validations)
	
	def form_setting_ls(self, rpart_cp_ls, loss_matrix_ls, prior_prob_ls, no_of_validations):
		"""
		11-19-05
			use no_of_validations to multiply the setting(separate the one setting's validations
				to different nodes)
			the extra setting copy  is for a non-validation real model fitting
		11-26-05
			make prior_prob None (proportional to observed data)
		11-29-05 add back prior_prob
		"""
		setting_ls = []
		for rpart_cp in rpart_cp_ls:
			for loss_matrix in loss_matrix_ls:
				for prior_prob in prior_prob_ls:
					setting_ls += [[rpart_cp, loss_matrix, prior_prob, 0]]*no_of_validations
					#11-19-05 0 is for validation
					setting_ls.append([rpart_cp, loss_matrix, prior_prob, 1])		#11-19-05 the extra copy
					#1 is for non-validation
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
		11-19-05
			setting is changed in form_setting_ls()
			one computing is just for one validation or non-validation
		11-23-05
			separate fit and prediction
			in validation, do prediction on unknown_data, return function-prediction set
		11-26-05 unknown_data doesn't need the prediction_set
		03-17-06
			rpart_prediction.py changed interface and added randomForest
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		unknown_data, known_data, training_perc, no_of_validations, model_type, bit_string = parameter_list
		rpart_cp, loss_matrix, prior_prob, type_code = data	#11-19-05 type_code
		parameter_list_for_fit = [rpart_cp, loss_matrix, prior_prob]	#03-17-06
		if type_code==0:
			training_data, testing_data = self.sample_data(known_data, training_perc)
			fit_model = self.fit_function_dict[model_type](training_data, parameter_list_for_fit, bit_string)
			pred_testing = self.predict_function_dict[model_type](fit_model, testing_data)
			pred_training = self.predict_function_dict[model_type](fit_model, training_data)
			pred_unknown = self.predict_function_dict[model_type](fit_model, unknown_data)	#11-23-05
			#03-17-06
			if model_type==2:	#randomForest's model has its own oob prediction
				fit_model_py = fit_model.as_py(rpy.BASIC_CONVERSION)
				training_acc_ls = self.cal_accuracy(training_data,  fit_model_py['predicted'], pred_type=1)
			else:
				training_acc_ls = self.cal_accuracy(training_data, pred_training, pred_type=model_type)
			
			result = [data,\
				self.cal_accuracy(testing_data, pred_testing, pred_type=model_type),\
				training_acc_ls,\
				self.cal_accuracy(unknown_data, pred_unknown, pred_type=model_type)]	#11-23-05	#11-26-05 don't need the prediction_set
		else:
			fit_model = self.fit_function_dict[model_type](known_data, parameter_list_for_fit, bit_string)
			known_pred = self.predict_function_dict[model_type](fit_model, known_data)
			unknown_pred = self.predict_function_dict[model_type](fit_model, unknown_data)
			#03-17-06
			if model_type==2:	#randomForest's model has its own oob prediction
				fit_model_py = fit_model.as_py(rpy.BASIC_CONVERSION)
				training_acc_ls = self.cal_accuracy(known_data,  fit_model_py['predicted'], pred_type=1)
			else:
				training_acc_ls = self.cal_accuracy(known_data, known_pred, pred_type=model_type)
			
			result = [data,\
				self.cal_accuracy(unknown_data, unknown_pred, pred_type=model_type),\
				training_acc_ls]
		"""
		testing_acc_ls, training_acc_ls = self.rpart_validation(known_data, training_perc, rpart_cp, \
			loss_matrix, prior_prob, no_of_validations)
		
		unknown_pred, known_pred = self.rpart_fit_and_predict(unknown_data, known_data, rpart_cp, loss_matrix, prior_prob)
		result = [data,\
			self.summary_stat_of_accuracy_ls(testing_acc_ls),\
			self.summary_stat_of_accuracy_ls(training_acc_ls),\
			self.cal_accuracy(unknown_data, unknown_pred),\
			self.cal_accuracy(known_data, known_pred)]
		"""
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
		11-23-05
			take the unknown result from validation
		11-26-05
			remove the data structure which stores the frequencies of predictions
		03-03-06
			fix a bug here
		"""
		writer, setting2validation_stat, setting2unknown_known_acc_ls, no_of_validations = parameter_list
		data = cPickle.loads(data)
		setting = data[0]
		key = (setting[0], tuple(setting[1]), setting[2])	#transform it to tuple
		if setting[3] == 0:	#type_code == 0(validation)
			if key not in setting2validation_stat:
				setting2validation_stat[key] = [[], [], [], {}]	#1st [] is for testing, 2nd [] is for training, 3rd [] is for unknown, 4th {} is for unknown prediction
			setting2validation_stat[key][0].append(data[1])	#testing
			setting2validation_stat[key][1].append(data[2])	#training
			setting2validation_stat[key][2].append(data[3][:4])	#unknown, only the 1st 4 entries, the 5th entry is below	11-23-05
			"""
			gene_no_ls = [row[0] for row in data[3][4]]	#11-27-05 check the frequency of unknown genes in all models
			gene_no_set = Set(gene_no_ls)
			for gene_no in gene_no_set:
				if gene_no not in setting2validation_stat[key][3]:
					setting2validation_stat[key][3][gene_no] = 0
				setting2validation_stat[key][3][gene_no] += 1
			"""
			"""
			for gene_no_go_no_pair in data[3][4]:	#11-23-05	#11-26-05
				if gene_no_go_no_pair not in setting2validation_stat[key][3]:	#11-23-05 check the frequency of unknown prediction-pairs
					setting2validation_stat[key][3][gene_no_go_no_pair] = 0
				setting2validation_stat[key][3][gene_no_go_no_pair] += 1
			"""
		else:
			setting2unknown_known_acc_ls[key] = [data[1], data[2]]
		if key in setting2validation_stat and len(setting2validation_stat[key][0])==no_of_validations and key in setting2unknown_known_acc_ls:
			#11-19-05 all validations are done and non-validation is also done
			testing_sum_stat = self.summary_stat_of_accuracy_ls(setting2validation_stat[key][0])
			training_sum_stat = self.summary_stat_of_accuracy_ls(setting2validation_stat[key][1])
			unknown_sum_stat = self.summary_stat_of_accuracy_ls(setting2validation_stat[key][2])	#11-23-05
			writer.writerow(list(key) +['testing'] + testing_sum_stat)
			writer.writerow(list(key) +['training'] + training_sum_stat)
			writer.writerow(list(key) +['unknown_valid'] + unknown_sum_stat)	#11-23-05
			writer.writerow(list(key) +['unknown']+ setting2unknown_known_acc_ls[key][0])
			writer.writerow(list(key) +['known']+ setting2unknown_known_acc_ls[key][1])
	
	def run(self):
		"""
		11-16-05
		11-19-05
			use no_of_validations to multiply the setting(separate the one setting's validations
				to different nodes)
			the extra setting copy  is for a non-validation real model fitting
			
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
			"""
			#11-19-05 shuffle data to check
			index_ls = range(len(known_data))
			random.shuffle(index_ls)
			for i in range(len(index_ls)):
				index_ls[i] = known_data[i]
			known_data = index_ls
			"""
			data, source, tag = communicator.receiveString(0, 0)
			unknown_data = cPickle.loads(data)	#take the data
			"""
			#11-19-05 shuffle data to check
			index_ls = range(len(unknown_data))
			random.shuffle(index_ls)
			for i in range(len(index_ls)):
				index_ls[i] = unknown_data[i]
			unknown_data = index_ls
			"""
		elif node_rank==communicator.size-1:
			writer = csv.writer(open(self.output_file, 'w'), delimiter='\t')
			#write down the header
			writer.writerow(['rpart_cp', 'loss_matrix', 'prior_prob', 'type', 'accuracy_avg','accuracy_std', 'no_of_predictions_avg',\
				'no_of_predictions_std', 'no_of_genes_avg', 'no_of_genes_std'])
			
		mpi_synchronize(communicator)
		if node_rank == 0:
			if self.type==1:
				setting_ls = self.form_setting_ls(self.rpart_cp_ls, self.loss_matrix_ls, self.prior_prob_ls, self.no_of_validations)
			elif self.type==2:
				#randomForest replaces rpart_cp_ls with mty_ls, others are ignored later
				setting_ls = self.form_setting_ls(self.mty_ls, self.loss_matrix_ls, self.prior_prob_ls, self.no_of_validations)
			else:
				sys.stderr.write("type %s not supported.\n"%self.type)
				sys.exit(3)
			self.input_node(communicator, setting_ls, free_computing_nodes, self.report)
		elif node_rank in free_computing_nodes:
			parameter_list = [unknown_data, known_data, self.training_perc, self.no_of_validations, self.type, self.bit_string]	#03-17-06 add type, bit_string
			computing_node(communicator, parameter_list, self.computing_handler, report=self.report)
		elif node_rank==communicator.size-1:
			setting2validation_stat = {}
			setting2unknown_known_acc_ls = {}
			parameter_list = [writer, setting2validation_stat, setting2unknown_known_acc_ls, self.no_of_validations]
			output_node(communicator, free_computing_nodes, parameter_list, self.output_handler, self.report)
			#cPickle.dump([setting2validation_stat, setting2unknown_known_acc_ls], open('/home/yuhuang/MpiRpartValidation.setting2result.pickle','w'))	#11-23-05
			del writer

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:j:f:y:p:l:o:s:x:t:m:b:gucr", ["help", "hostname=", \
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
	type = 1
	mty_ls = [0]
	bit_string = '1111111'
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
		elif opt in ("-i",):
			fname = arg
		elif opt in ("-j",):
			output_file = arg
		elif opt in ("-f",):
			filter_type = int(arg)
		elif opt in ("-y",):
			is_correct_type = int(arg)
		elif opt in ("-p",):
			rpart_cp_ls = arg.split(',')
			rpart_cp_ls = map(float, rpart_cp_ls)
		elif opt in ("-l",):
			loss_matrix_ls = arg.split('=')
			for i in range(len(loss_matrix_ls)):
				loss_matrix_ls[i] = map(float, loss_matrix_ls[i].split(','))
		elif opt in ("-o",):
			prior_prob_ls = map(float, arg.split(','))
		elif opt in ("-s",):
			training_perc = float(arg)
		elif opt in ("-x",):
			no_of_validations = int(arg)
		elif opt in ("-t",):
			type = int(arg)
		elif opt in ("-m",):
			mty_ls = map(int, arg.split(','))
		elif opt in ("-b",):
			bit_string = arg
		elif opt in ("-g",):
			need_cal_hg_p_value = 1
		elif opt in ("-u",):
			debug = 1
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-r",):
			report = 1
	if schema and fname and output_file:
		instance = MpiRpartValidation(hostname, dbname, schema, fname, output_file, \
			filter_type, is_correct_type, rpart_cp_ls, loss_matrix_ls, prior_prob_ls, training_perc, \
			no_of_validations, type, mty_ls, bit_string, need_cal_hg_p_value, debug, commit, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
