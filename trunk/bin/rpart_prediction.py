#!/usr/bin/env python
"""
Usage: rpart_prediction.py -k -i -j

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-i ...,	fname of setting1
	-l ...,	lm_bit of setting1, ('00001', default)
	-j ...,	fname of setting2
	-m ...,	lm_bit of setting2, ('00001' default)
	-f ...,	filter type
	-y ...,	is_correct type (2 lca, default)
	-b,	enable debug flag
	-c,	commit the database transaction
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	~/script/annot/bin/rpart_prediction.py -k hs_fim_40 -i hs_fim_40m4x40rec0_8 -l 11100
		-j hs_fim_40m4x40e1s1__1_0l10 -m 11001 -c -r

Description:
	Program to do rpart fittng and prediction.
"""
import sys, os, getopt, csv
sys.path += [os.path.expanduser('~/script/annot/bin')]
from codense.common import db_connect, form_schema_tables
from MpiPredictionFilter import prediction_attributes, MpiPredictionFilter
from numarray import array
from rpy import r, set_default_mode,NO_CONVERSION,BASIC_CONVERSION


class rpart_prediction:
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, fname1=None, \
		lm_bit1=None, fname2=None, lm_bit2=None, filter_type=1, is_correct_type=2, \
		debug=0, commit=0, report=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.fname1 = fname1
		self.lm_bit1 = lm_bit1
		self.fname2 = fname2
		self.lm_bit2 = lm_bit2
		self.filter_type = int(filter_type)
		self.is_correct_type = int(is_correct_type)
		self.debug = int(debug)
		self.commit = int(commit)
		self.report = int(report)
	
	def data_fetch(self, curs, p_gene_table, filter_type, is_correct_type):
		"""
		11-09-05
			1st get the data from p_gene_table and remove redundancy given filter_type
			2nd transform the data to three lists
		"""
		sys.stderr.write("Fetching data from old p_gene_table...\n")
		prediction_pair2instance = {}
		curs.execute("DECLARE crs CURSOR FOR SELECT p.p_gene_id, p.gene_no, p.go_no, p.is_correct, p.is_correct_l1, \
			p.is_correct_lca, p.avg_p_value, p.no_of_clusters, p.cluster_array, p.p_value_cut_off, p.recurrence_cut_off, \
			p.connectivity_cut_off, p.cluster_size_cut_off, p.unknown_cut_off, p.depth_cut_off, p.mcl_id, p.lca_list, p.vertex_gradient,\
			p.edge_gradient from %s p"%(p_gene_table))
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
			is_correct = p_attr_instance.is_correct_dict[is_correct_type]
			data_row = [p_attr_instance.p_value_cut_off, p_attr_instance.recurrence_cut_off, p_attr_instance.connectivity_cut_off,\
				p_attr_instance.cluster_size_cut_off, p_attr_instance.edge_gradient, is_correct]
			all_data.append(data_row)
			if is_correct!=-1:
				known_data.append(data_row)	#to do fitting
		all_data = array(all_data)
		known_data = array(known_data)
		
		sys.stderr.write("Done fetching data.\n")
		return prediction_ls, all_data, known_data
	
	def rpart_fit_and_predict(self, all_data, known_data, bit_string='11111'):
		"""
		11-09-05
			1st use known_data to get the fit model
			2nd use the fit model to do prediction on all_data, result is prob for each class
		"""
		sys.stderr.write("rpart fitting and predicting...\n")
		r.library("rpart")
		coeff_name_list = ['p_value', 'recurrence', 'connectivity', 'cluster_size', 'gradient']
		formula_list = []
		for i in range(len(bit_string)):
			if bit_string[i] == '1':
				formula_list.append(coeff_name_list[i])
		set_default_mode(NO_CONVERSION)
		data_frame = r.as_data_frame({"p_value":known_data[:,0], "recurrence":known_data[:,1], "connectivity":known_data[:,2], \
			"cluster_size":known_data[:,3], "gradient":known_data[:,4], "is_correct":known_data[:,-1]})
		fit = r.rpart(r("is_correct~%s"%'+'.join(formula_list)), data=data_frame, method="class")
		del data_frame
		all_data_frame = r.as_data_frame({"p_value":all_data[:,0], "recurrence":all_data[:,1], "connectivity":all_data[:,2], \
			"cluster_size":all_data[:,3], "gradient":all_data[:,4], "is_correct":all_data[:,-1]})
		set_default_mode(BASIC_CONVERSION) #04-07-05
		pred = r.predict(fit, all_data_frame)
		sys.stderr.write("Done rpart fitting and predicting.\n")
		return pred
	
	def record_data(self, curs, MpiPredictionFilter_instance, prediction_ls, pred, new_p_gene_tablle):
		"""
		11-09-05
		"""
		sys.stderr.write("Recoding prediction...\n")
		for i in range(len(prediction_ls)):
			p_attr_instance = prediction_ls[i]
			p_attr_instance.vertex_gradient = pred[i][1]
			MpiPredictionFilter_instance.submit_to_p_gene_table(curs, new_p_gene_tablle, p_attr_instance)
		sys.stderr.write("Done recoding prediction...\n")
	
	def run(self):
		"""
		11-09-05
		"""
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		old_schema_instance = form_schema_tables(self.fname1, self.lm_bit1)
		new_schema_instance = form_schema_tables(self.fname2, self.lm_bit2)
		
		prediction_ls, all_data, known_data = self.data_fetch(curs, old_schema_instance.p_gene_table, self.filter_type, self.is_correct_type)
		pred = self.rpart_fit_and_predict(all_data, known_data)
		MpiPredictionFilter_instance = MpiPredictionFilter()
		MpiPredictionFilter_instance.view_from_table(curs, old_schema_instance.splat_table, new_schema_instance.splat_table)
		MpiPredictionFilter_instance.view_from_table(curs, old_schema_instance.mcl_table, new_schema_instance.mcl_table)
		MpiPredictionFilter_instance.view_from_table(curs, old_schema_instance.pattern_table, new_schema_instance.pattern_table)
		MpiPredictionFilter_instance.createGeneTable(curs, new_schema_instance.p_gene_table)
		self.record_data(curs, MpiPredictionFilter_instance, prediction_ls, pred, new_schema_instance.p_gene_table)
		if self.commit:
			curs.execute("end")

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:l:j:m:f:y:bcr", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	fname1 = None
	lm_bit1 = '00001'
	fname2 = None
	lm_bit2 = '00001'
	filter_type = 1
	is_correct_type = 2
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
		elif opt in ("-l"):
			lm_bit1 = arg
		elif opt in ("-j"):
			fname2 = arg
		elif opt in ("-m"):
			lm_bit2 = arg
		elif opt in ("-f"):
			filter_type = int(arg)
		elif opt in ("-y"):
			is_correct_type = int(arg)
		elif opt in ("-b"):
			debug = 1
		elif opt in ("-c"):
			commit = 1
		elif opt in ("-r"):
			report = 1
	if schema and fname1 and fname2:
		instance = rpart_prediction(hostname, dbname, schema, fname1, lm_bit1, \
			fname2, lm_bit2, filter_type, is_correct_type, debug, commit, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
