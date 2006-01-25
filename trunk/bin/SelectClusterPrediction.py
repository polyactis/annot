#!/usr/bin/env python
"""
Usage: codense2db.py -k SCHEMA -i -s -j -m [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-i ...,	cluster filename, output of MpiBFSCluster.py
	-s ...,	stat filename, output of MpiStatCluster.py
	-j ...,	fname setting to create splat, mcl, pattern tables
	-f ...,	filter type (same as rpart_prediction.py) 1(default)
	-m ...,	gim(gene incidence matrix) inputfile(-y=4 only)
	-b,	debug version.
	-c,	commit this database transaction
	-r,	report the progress(a number)
	-h,	show this help
	
Examples:

Description:
	Program to select clusters from MpiBFSCluster.py and MpiStatCluster.py's
		output and submit clusters and predictions to database.
	filter_type:
		1. longest recurrence
		2. largest edge_gradient
		3. largest (recurrence + edge_gradient)
	
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import getopt, csv, re
from codense.common import db_connect, get_known_genes_dict, form_schema_tables,\
	get_gene_no2incidence_array, get_vertex_set_gim_array, get_gene_id2gene_no, \
	get_no_of_total_genes, cal_hg_p_value, get_go_no2gene_no_set
from sets import Set
from MpiPredictionFilter import prediction_attributes, MpiPredictionFilter
from codense.codense2db import codense2db
from rpy import r
if sys.version_info[:2] < (2, 3):       #python2.2 or lower needs some extra
	from python2_3 import *

class SelectClusterPrediction:
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, inputfile=None,\
		stat_fname=None, jnput_fname=None, filter_type=1, gim_inputfname=None, \
		debug=0, commit=0, report=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.inputfile = inputfile
		self.stat_fname = stat_fname
		self.jnput_fname = jnput_fname
		self.filter_type = int(filter_type)
		self.gim_inputfname = gim_inputfname
		self.debug = int(debug)
		self.commit = int(commit)
		self.report = int(report)
	
	def parse_stat_fname(self, stat_fname, filter_type):
		sys.stderr.write("Parsing stat_fname: %s ...\n"%os.path.basename(stat_fname))
		reader = csv.reader(open(stat_fname), delimiter='\t')
		prediction_pair2instance = {}
		
		#following is temp default values
		p_gene_id = -2	#-2 means no database submission for p_gene_id
		avg_p_value = -1
		no_of_clusters = 1
		p_value_cut_off = -1
		connectivity_cut_off = -1
		cluster_size_cut_off = -1
		unknown_cut_off = -1
		vertex_gradient = -1
		
		counter = 0
		real_counter = 0
		for row in reader:
			cluster_id, gene_no, go_no, go_no_depth, recurrence, gradient_score, edge_gradient, \
				is_correct, is_correct_L1, is_correct_lca, lca_list_string = row
			cluster_id = int(cluster_id)
			gene_no = int(gene_no)
			go_no = int(go_no)
			go_no_depth = int(go_no_depth)
			recurrence = float(recurrence)
			edge_gradient = float(edge_gradient)
			is_correct = int(is_correct)
			is_correct_L1 = int(is_correct_L1)
			is_correct_lca = int(is_correct_lca)
			cluster_array = '{%s}'%cluster_id
			new_row = [p_gene_id, gene_no, go_no, is_correct, is_correct_L1, \
				is_correct_lca, avg_p_value, no_of_clusters, cluster_array, p_value_cut_off, recurrence, \
				connectivity_cut_off, cluster_size_cut_off, unknown_cut_off, go_no_depth, cluster_id, \
				lca_list_string, vertex_gradient, edge_gradient]
			p_attr_instance = prediction_attributes(new_row, type=2)
			prediction_pair = (gene_no, go_no)
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
			if self.report and counter%10000==0:
				sys.stderr.write("%s%s/%s"%('\x08'*20, counter, real_counter))
		if self.report:
			sys.stderr.write("%s%s/%s"%('\x08'*20, counter, real_counter))
		del reader
		sys.stderr.write("Done.\n")
		return prediction_pair2instance
	
	def get_cluster_id_set(self, prediction_pair2instance):
		sys.stderr.write("Getting cluster_id_set from prediction_pair2instance...")
		cluster_id_set = Set()
		for prediction_pair,p_attr_instance in prediction_pair2instance.iteritems():
			cluster_id_set.add(p_attr_instance.mcl_id)
		sys.stderr.write("done with %s clusters.\n"%len(cluster_id_set))
		return cluster_id_set
	
	def parse_cluster_fname(self, curs, cluster_fname, gim_inputfname, cluster_id_set, schema_instance):
		"""
		01-24-06
			a lot of analogy to codense2db.py's run()
		"""
		sys.stderr.write("Parsing cluster_fname: %s ...\n"%os.path.basename(cluster_fname))
		codense2db_instance  = codense2db()
		codense2db_instance.create_tables(curs, schema_instance.splat_table, \
			schema_instance.mcl_table, schema_instance.pattern_table)
		gene_id2gene_no = get_gene_id2gene_no(curs)
		gene_no2incidence_array = get_gene_no2incidence_array(gim_inputfname, gene_id2gene_no)
		known_gene_no2go_no_set = get_known_genes_dict(curs)
		counter = 0
		real_counter = 0
		cluster_id2properties = {}	#additional properties for prediction_pair2instance
		reader = csv.reader(open(cluster_fname, 'r'), delimiter='\t')
		for row in reader:
			counter += 1
			#only those who are in cluster_id_set
			if counter in cluster_id_set:	#cluster_id starts from 1
				cluster_list = codense2db_instance.fimbfs_parser(row, gene_no2incidence_array, curs)
				for cluster in cluster_list:
					real_counter += 1
					cluster.unknown_gene_ratio = codense2db_instance.calculate_unknown_gene_ratio(cluster.vertex_set, \
						known_gene_no2go_no_set)
					cluster.cluster_id = counter	#line number is the cluster_id
					codense2db_instance.db_submit(curs, cluster, schema_instance.pattern_table)
					
					cluster_id2properties[cluster.cluster_id] = [cluster.connectivity, cluster.unknown_gene_ratio, cluster.vertex_set]
			if real_counter==len(cluster_id_set):
				#all relevant clusters have been got, ignore remaining clusters
				break
			if self.report and counter%2000==0:
				sys.stderr.write("%s%s/%s"%('\x08'*20, counter, real_counter))
		if self.report:
			sys.stderr.write("%s%s/%s"%('\x08'*20, counter, real_counter))
		del reader
		sys.stderr.write("Done.\n")
		return cluster_id2properties
		
	def submit_predictions(self, curs, schema_instance, prediction_pair2instance, cluster_id2properties):
		sys.stderr.write("Submitting predictions...\n")
		MpiPredictionFilter_instance = MpiPredictionFilter()
		MpiPredictionFilter_instance.createGeneTable(curs, schema_instance.p_gene_table)
		
		no_of_total_genes = get_no_of_total_genes(curs)
		go_no2gene_no_set = get_go_no2gene_no_set(curs)
		counter = 0
		for prediction_pair, p_attr_instance in prediction_pair2instance.iteritems():
			#1st fill those empty items
			properties = cluster_id2properties[p_attr_instance.mcl_id]
			vertex_set = properties[2]
			p_attr_instance.p_value_cut_off = cal_hg_p_value(p_attr_instance.gene_no, p_attr_instance.go_no,\
				vertex_set, no_of_total_genes, go_no2gene_no_set, r)
			p_attr_instance.avg_p_value = p_attr_instance.p_value_cut_off
			p_attr_instance.connectivity_cut_off = properties[0]
			p_attr_instance.cluster_size_cut_off = len(vertex_set)
			p_attr_instance.unknown_cut_off = properties[1]
			MpiPredictionFilter_instance.submit_to_p_gene_table(curs, schema_instance.p_gene_table, p_attr_instance)
			counter += 1
			if self.report and counter%2000==0:
				sys.stderr.write("%s%s"%('\x08'*20, counter))
		if self.report:
			sys.stderr.write("%s%s"%('\x08'*20, counter))
		sys.stderr.write("Done.\n")
		
	def run(self):
		"""
		01-24-06
			
		"""
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		schema_instance = form_schema_tables(self.jnput_fname)
		prediction_pair2instance = self.parse_stat_fname(self.stat_fname, self.filter_type)
		cluster_id_set = self.get_cluster_id_set(prediction_pair2instance)
		cluster_id2properties = self.parse_cluster_fname(curs, self.inputfile, self.gim_inputfname, cluster_id_set, schema_instance)
		self.submit_predictions(curs, schema_instance, prediction_pair2instance, cluster_id2properties)
		if self.commit:
			curs.execute("end")

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:s:j:f:m:bcr", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	inputfile = None
	stat_fname = None
	jnput_fname = None
	filter_type = 1
	gim_inputfname = None
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
			inputfile = arg
		elif opt in ("-s",):
			stat_fname = arg
		elif opt in ("-j",):
			jnput_fname = arg
		elif opt in ("-f",):
			filter_type = int(arg)
		elif opt in ("-m",):
			gim_inputfname = arg
		elif opt in ("-b",):
			debug = 1
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-r",):
			report = 1
	if schema and inputfile and stat_fname and jnput_fname:
		instance = SelectClusterPrediction(hostname, dbname, schema, inputfile, stat_fname, \
			jnput_fname, filter_type, gim_inputfname, debug, commit, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
