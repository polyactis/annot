#!/usr/bin/env python
"""
Usage: rpart_prediction.py -k -i -j

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-i ...,	fname of setting1
	-j ...,	fname of setting2
	-f ...,	filter type
	-y ...,	is_correct type (2 lca, default)
	-p ...,	rpart cp value(0.01, default)
	-l ...,	loss matrix, 0,1,1,0 (default)
	-o ...,	prior prob for the 1st class, None (default, i.e. 0.65)
	-s ...,	no_of_buckets selected to do training, 5(default)
	-g, 	calculate the hypergeometric p-value to replace p_value_cut_off(gradient)
	-b,	enable debug flag
	-c,	commit the database transaction
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	~/script/annot/bin/rpart_prediction.py -k hs_fim_40 -i hs_fim_40m4x40rec0_8
		-j hs_fim_40m4x40e1s1__1_0l10 -c -r

Description:
	Program to do rpart fittng and prediction.
"""
import sys, os, getopt, csv, random, math
sys.path += [os.path.expanduser('~/script/annot/bin')]
from codense.common import db_connect, form_schema_tables, cal_hg_p_value, \
	get_go_no2gene_no_set, get_no_of_total_genes
from MpiPredictionFilter import prediction_attributes, MpiPredictionFilter
from numarray import array
from rpy import r, set_default_mode,NO_CONVERSION,BASIC_CONVERSION
from sets import Set

class rpart_prediction:
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, fname1=None, fname2=None, \
		filter_type=1, is_correct_type=2, rpart_cp=0.01, loss_matrix=[0,1,1,0], prior_prob=None, no_of_buckets=5,\
		need_cal_hg_p_value=0, debug=0, commit=0, report=0):
		"""
		11-09-05 add rpart_cp
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.fname1 = fname1
		self.fname2 = fname2
		self.filter_type = int(filter_type)
		self.is_correct_type = int(is_correct_type)
		self.rpart_cp = rpart_cp
		self.loss_matrix = loss_matrix
		self.prior_prob = prior_prob
		self.no_of_buckets = int(no_of_buckets)
		self.need_cal_hg_p_value = int(need_cal_hg_p_value)
		self.debug = int(debug)
		self.commit = int(commit)
		self.report = int(report)
	
	def data_fetch(self, curs, schema_instance, filter_type, is_correct_type, \
		no_of_total_genes, go_no2gene_no_set, need_cal_hg_p_value=0):
		"""
		11-09-05
			1st get the data from p_gene_table and remove redundancy given filter_type
			2nd transform the data to three lists
		11-10-05 add a chunk of code to get hg p-value(leave one out) for the prediction
			mcl_id2vertex_list might blow the memory.(?) 
		"""
		sys.stderr.write("Fetching data from old p_gene_table...\n")
		prediction_pair2instance = {}
		curs.execute("DECLARE crs CURSOR FOR SELECT p.p_gene_id, p.gene_no, p.go_no, p.is_correct, p.is_correct_l1, \
			p.is_correct_lca, p.avg_p_value, p.no_of_clusters, p.cluster_array, p.p_value_cut_off, p.recurrence_cut_off, \
			p.connectivity_cut_off, p.cluster_size_cut_off, p.unknown_cut_off, p.depth_cut_off, p.mcl_id, p.lca_list, p.vertex_gradient,\
			p.edge_gradient from %s p"%(schema_instance.p_gene_table))
		curs.execute("fetch 10000 from crs")
		rows = curs.fetchall()
		counter = 0
		real_counter = 0
		while rows:
			for row in rows:
				p_attr_instance = prediction_attributes(row, type=2)
				prediction_pair = (p_attr_instance.gene_no, p_attr_instance.go_no)
				if prediction_pair not in prediction_pair2instance:
					prediction_pair2instance[prediction_pair] = p_attr_instance
					real_counter += 1
				else:	#remove redundancy
					if filter_type==1:
						new_cmp_value = p_attr_instance.recurrence_cut_off
						old_cmp_value = prediction_pair2instance[prediction_pair].recurrence_cut_off
					elif filter_type==2:
						new_cmp_value = p_attr_instance.edge_gradient
						old_cmp_value = prediction_pair2instance[prediction_pair].edge_gradient
					elif filter_type==3:
						new_cmp_value = p_attr_instance.recurrence_cut_off+p_attr_instance.edge_gradient
						old_cmp_value = prediction_pair2instance[prediction_pair].recurrence_cut_off + prediction_pair2instance[prediction_pair].edge_gradient
					if new_cmp_value>old_cmp_value:
						prediction_pair2instance[prediction_pair] = p_attr_instance
				counter += 1
			if self.report:
				sys.stderr.write("%s%s/%s"%('\x08'*20, counter, real_counter))
			curs.execute("fetch 10000 from crs")
			rows = curs.fetchall()
		prediction_ls = []
		all_data = []
		known_data = []
		for prediction_pair,p_attr_instance in prediction_pair2instance.iteritems():
			prediction_ls.append(p_attr_instance)
			#11-10-05
			mcl_id2vertex_list = {}
			if need_cal_hg_p_value:
				mcl_id = p_attr_instance.mcl_id
				if mcl_id not in mcl_id2vertex_list:
					mcl_id2vertex_list[mcl_id] = self.get_vertex_list(curs, schema_instance, mcl_id)
				p_attr_instance.p_value_cut_off = cal_hg_p_value(p_attr_instance.gene_no, p_attr_instance.go_no,\
					mcl_id2vertex_list[mcl_id], no_of_total_genes, go_no2gene_no_set, r)
			
			is_correct = p_attr_instance.is_correct_dict[is_correct_type]
			data_row = [p_attr_instance.p_value_cut_off, p_attr_instance.recurrence_cut_off, p_attr_instance.connectivity_cut_off,\
				p_attr_instance.cluster_size_cut_off, p_attr_instance.edge_gradient, p_attr_instance.gene_no, p_attr_instance.go_no, \
				is_correct]
			all_data.append(data_row)
			if is_correct!=-1:
				known_data.append(data_row)	#to do fitting
		
		sys.stderr.write("Done fetching data.\n")
		return prediction_ls, all_data, known_data
	
	def rpart_fit_and_predict(self, all_data, known_data, rpart_cp, loss_matrix, prior_prob, bit_string='11111'):
		"""
		11-09-05
			1st use known_data to get the fit model
			2nd use the fit model to do prediction on all_data, result is prob for each class
		11-09-05 add rpart_cp
		11-17-05
			add loss_matrix, prior_prob
			return two pred
		"""
		sys.stderr.write("rpart fitting and predicting...\n")
		r.library("rpart")
		coeff_name_list = ['p_value', 'recurrence', 'connectivity', 'cluster_size', 'gradient']
		formula_list = []
		for i in range(len(bit_string)):
			if bit_string[i] == '1':
				formula_list.append(coeff_name_list[i])
		#11-17-05 transform into array
		all_data = array(all_data)
		known_data = array(known_data)
		
		set_default_mode(NO_CONVERSION)
		data_frame = r.as_data_frame({"p_value":known_data[:,0], "recurrence":known_data[:,1], "connectivity":known_data[:,2], \
			"cluster_size":known_data[:,3], "gradient":known_data[:,4], "is_correct":known_data[:,-1]})
		if prior_prob:
			prior_prob = [prior_prob, 1-prior_prob]	#get the full list
			fit = r.rpart(r("is_correct~%s"%'+'.join(formula_list)), data=data_frame, method="class", control=r.rpart_control(cp=rpart_cp),\
				parms=r.list(prior=prior_prob, loss=r.matrix(loss_matrix) ) )
		else:
			fit = r.rpart(r("is_correct~%s"%'+'.join(formula_list)), data=data_frame, method="class", control=r.rpart_control(cp=rpart_cp),\
				parms=r.list(loss=r.matrix(loss_matrix) ) )
		
		set_default_mode(BASIC_CONVERSION)
		pred_training = r.predict(fit, data_frame, type=["class"])
		del data_frame
		
		set_default_mode(NO_CONVERSION)
		all_data_frame = r.as_data_frame({"p_value":all_data[:,0], "recurrence":all_data[:,1], "connectivity":all_data[:,2], \
			"cluster_size":all_data[:,3], "gradient":all_data[:,4], "is_correct":all_data[:,-1]})
		set_default_mode(BASIC_CONVERSION)
		pred = r.predict(fit, all_data_frame, type=["class"])	#11-17-05 type=c("class")
		del all_data_frame
		sys.stderr.write("Done rpart fitting and predicting.\n")
		return pred, pred_training
	
	def get_vertex_list(self, curs, schema_instance, mcl_id):
		"""
		11-10-05 for data_fetch() to cal_hg_p_value
		"""
		curs.execute("select vertex_set from %s where id=%s"%(schema_instance.pattern_table, mcl_id))
		rows = curs.fetchall()
		vertex_list = rows[0][0]
		vertex_list = vertex_list[1:-1].split(',')
		vertex_list = map(int, vertex_list)
		return vertex_list
	
	def record_data(self, curs, MpiPredictionFilter_instance, prediction_ls, pred, schema_instance):
		"""
		11-09-05
		11-18-05 pred is type="class"
		"""
		sys.stderr.write("Recording prediction...\n")
		for i in range(len(prediction_ls)):
			p_attr_instance = prediction_ls[i]
			p_attr_instance.vertex_gradient = int(pred[repr(i+1)])	#WATCH R's index starts from 1
			MpiPredictionFilter_instance.submit_to_p_gene_table(curs, schema_instance.p_gene_table, p_attr_instance)
		sys.stderr.write("Done recoding prediction...\n")
	
	####11-17-05  for parameter tuning purpose
	def random_partition(self, total_no_of_cases, no_of_buckets):
		"""
		11-18-05
			randomly partition the data into several buckets
		"""
		case_indices = range(total_no_of_cases)
		random.shuffle(case_indices)	#shuffle it
		bucket_size = int(math.ceil(total_no_of_cases/float(no_of_buckets)))	#make it a little bit larger. WATCH int() coersion
		bucket_ls = []
		for i in range(no_of_buckets):	#create empty sets
			#WATCH: don't use [Set()]*no_of_buckets. When you add one item to one of these sets, all of them get it.
			#They are pointer-copies which all points to the same memory location.
			bucket_ls.append(Set())
		for i in range(total_no_of_cases):
			bucket_ls[i/bucket_size].add(case_indices[i])	#i/bucket_size
		return bucket_ls
	
	def sample_data(self, data,  bucket):
		"""
		11-17-05
			for parameter tuning purpose
		11-18-05
			use bucket got from random_partition(), (old way is to use training_perc to sample all data)
		"""
		total_no_of_cases = len(data)
		training_data = []
		testing_data  = []
		for i in range(total_no_of_cases):
			if i in bucket:	#bucket is separated to do validation as testing_data
				testing_data.append(data[i])
			else:
				training_data.append(data[i])
		return training_data, testing_data
	
	def cal_accuracy(self, data, pred):
		"""
		11-17-05
			for parameter tuning purpose
		"""
		no_of_total_predictions = 0
		no_of_correct_predictions = 0
		gene_set = Set()
		for i in range(len(data)):
			p_value, recurrence, connectivity, cluster_size, edge_gradient, gene_no, go_no, is_correct = data[i]
			key_in_pred = repr(i+1)	#WATCH R's index starts from 1
			if pred[key_in_pred] == '1':
				no_of_total_predictions += 1
				gene_set.add(gene_no)
				if is_correct==1:
					no_of_correct_predictions += 1
		if no_of_total_predictions>0:
			accuracy = float(no_of_correct_predictions)/no_of_total_predictions
		else:
			accuracy = 0
		return [accuracy, no_of_correct_predictions, no_of_total_predictions, len(gene_set)]
	
	def rpart_validation(self, known_data, no_of_buckets, rpart_cp, loss_matrix, prior_prob):
		"""
		11-17-05
			for parameter tuning purpose
		11-18-05 use random_partition()
			
			--random_partition()
			(loop)
				--sample_data()
				--cal_accuracy()
		"""
		training_acc_ls = []
		testing_acc_ls = []
		bucket_ls = self.random_partition(len(known_data), no_of_buckets)
		for bucket in bucket_ls:
			training_data, testing_data = self.sample_data(known_data, bucket)
			pred_testing, pred_training = self.rpart_fit_and_predict(testing_data, training_data, rpart_cp, loss_matrix, prior_prob)
			testing_acc_ls.append(self.cal_accuracy(testing_data, pred_testing))
			training_acc_ls.append(self.cal_accuracy(training_data, pred_training))
		return testing_acc_ls, training_acc_ls
		
	def run(self):
		"""
		11-09-05
		11-09-05 add rpart_cp
		11-10-05 add need_cal_hg_p_value
		
			--db_connect()
			--form_schema_tables()
			--form_schema_tables()
			--get_no_of_total_genes()
			--get_go_no2gene_no_set()
			--data_fetch()
				--get_vertex_list()
				--cal_hg_p_value()
			--rpart_fit_and_predict()
			--MpiPredictionFilter_instance....()
			--record_data()
		"""
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		old_schema_instance = form_schema_tables(self.fname1)
		new_schema_instance = form_schema_tables(self.fname2)
		
		no_of_total_genes = get_no_of_total_genes(curs)
		go_no2gene_no_set = get_go_no2gene_no_set(curs)
		
		prediction_ls, all_data, known_data = self.data_fetch(curs, old_schema_instance, self.filter_type, self.is_correct_type, \
			no_of_total_genes, go_no2gene_no_set, need_cal_hg_p_value)
		
		"""
		testing_acc_ls, training_acc_ls = self.rpart_validation(known_data, self.no_of_buckets, self.rpart_cp, \
			self.loss_matrix, self.prior_prob)
		print testing_acc_ls
		print training_acc_ls
		"""
		pred, pred_training = self.rpart_fit_and_predict(all_data, known_data, self.rpart_cp, self.loss_matrix, self.prior_prob)
		
		MpiPredictionFilter_instance = MpiPredictionFilter()
		MpiPredictionFilter_instance.view_from_table(curs, old_schema_instance.splat_table, new_schema_instance.splat_table)
		MpiPredictionFilter_instance.view_from_table(curs, old_schema_instance.mcl_table, new_schema_instance.mcl_table)
		MpiPredictionFilter_instance.view_from_table(curs, old_schema_instance.pattern_table, new_schema_instance.pattern_table)
		MpiPredictionFilter_instance.createGeneTable(curs, new_schema_instance.p_gene_table)
		self.record_data(curs, MpiPredictionFilter_instance, prediction_ls, pred, new_schema_instance)
		if self.commit:
			curs.execute("end")
		

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:j:f:y:p:l:o:s:gbcr", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	fname1 = None
	fname2 = None
	filter_type = 1
	is_correct_type = 2
	rpart_cp = 0.01
	loss_matrix = [0,1,1,0]
	prior_prob = None
	no_of_buckets = 5
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
			fname1 = arg
		elif opt in ("-j"):
			fname2 = arg
		elif opt in ("-f"):
			filter_type = int(arg)
		elif opt in ("-y"):
			is_correct_type = int(arg)
		elif opt in ("-p"):
			rpart_cp = float(arg)
		elif opt in ("-l"):
			loss_matrix = map(float, arg.split(','))
		elif opt in ("-o"):
			prior_prob = float(arg)
		elif opt in ("-s"):
			no_of_buckets = int(arg)
		elif opt in ("-g"):
			need_cal_hg_p_value = 1
		elif opt in ("-b"):
			debug = 1
		elif opt in ("-c"):
			commit = 1
		elif opt in ("-r"):
			report = 1
	if schema and fname1 and fname2:
		instance = rpart_prediction(hostname, dbname, schema, fname1, fname2, \
			filter_type, is_correct_type, rpart_cp, loss_matrix, prior_prob, no_of_buckets, \
			need_cal_hg_p_value, debug, commit, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
