#!/usr/bin/env python
"""
Usage: context_specific.py -k SCHEMA -p P_VALUE_CUT_OFF [OPTION] [STAT_TABLE_FILE]

Option:
	STAT_TABLE_FILE is the file to store the function table of all predicted genes.
		If it's not given, no table will be outputed.
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	cluster_stat(default), IGNORE
	-m ..., --mcl_table=...	mcl_result(default), mcl_result table corresponding to above table. IGNORE
	-g ..., --gene_table=...	table storing the stat results, p_gene(default)
	-s ..., --contrast=...	the contrast differential ratio for two contexts to be different, 0.3(default)
	-r, --report	report the progress(a number) IGNORE
	-c, --commit	commit the database transaction, records in table gene. IGNORE
	-h, --help              show this help

Examples:
	context_specific.py -k sc_yh60 -g p_gene  context_table_file
	context_specific.py -k sc_yh60 -g p_gene -s 0.4

Description:
	
"""

import sys, os, psycopg, getopt, csv
from graphlib import Graph, GraphAlgo
from sets import Set
from rpy import r

class gene_prediction:
	# class holding prediction information of a gene
	def __init__(self):
		self.tp = {}
		self.tp1 = {}
		self.tn = 0
		self.fp = {}
		self.fn = {}
		self.p_functions_dict = {}
		self.p_functions_struc_dict = {}
		self.mcl_id_list = []

class function_struc:
	#data structure for p_functions_struc_dict in gene_prediction
	def __init__(self):
		self.is_correct = 0
		self.p_value_list = []
		self.cluster_array = []
		self.context_dict = {}

class accuracy_struc:
	#data structure for self.go_no2accuracy
	def __init__(self):
		self.all = 0
		self.correct = 0
		self.ratio = 0

class context_specific:
	def __init__(self, hostname, dbname, schema, table, mcl_table, contrast=0.3, report=0, \
		needcommit=0, gene_table='p_gene', stat_table_fname='null'):
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.table = table
		self.mcl_table = mcl_table
		self.contrast = float(contrast)
		self.report = int(report)
		self.needcommit = int(needcommit)
		self.gene_table = gene_table
		self.stat_table_fname = stat_table_fname
		
		#mapping between gene_no and go_no list
		self.known_genes_dict = {}
		#some counters
		self.no_of_predicted_genes = 0
		#only >1 numbers
		self.list_no_of_functions_each_gene = []
		self.list_no_of_contexts_each_gene = []
		self.no_of_records = 0
		self.log_file = open('/tmp/context_specific.log','w')
		self.gene_prediction_dict = {}
		self.no_of_p_known = 0
		self.go_no2go_id = {}
		#mapping between go_no and go's name
		self.go_no2go_name = {}
		#mapping between gene_no and gene_id
		self.gene_no2gene_id = {}
		#the overall prediction accuracy of each function
		self.go_no2accuracy = {}
		#GO DAG
		self.go_graph = Graph.Graph()
		#a list of the depths of all the predicted functions, for histogram
		self.depth_list=[]
		#a list storing the maximum distance between the predicted functions for one gene.
		self.max_distance_list = []
		
	def dstruc_loadin(self):
		sys.stderr.write("Loading Data STructure...")
		#setup self.known_genes_dict
		self.curs.execute("select gene_no,go_functions from gene where known=TRUE")
		rows = self.curs.fetchall()
		for row in rows:
			go_functions_list = row[1][1:-1].split(',')
			self.known_genes_dict[row[0]] = []
			for go_no in go_functions_list:
				self.known_genes_dict[row[0]].append(int(go_no))
		
		#setup self.go_no2go_id
		self.curs.execute("select go_no, go_id from go")
		rows = self.curs.fetchall()
		for row in rows:
			self.go_no2go_id[row[0]] = row[1]
		
		#setup self.go_no2go_name
		self.curs.execute("select g.go_no,t.name from go g, go.term t where g.go_id=t.acc")
		rows = self.curs.fetchall()
		for row in rows:
			self.go_no2go_name[row[0]] = row[1]
				
		#setup self.gene_no2gene_id
		self.curs.execute("select gene_no, gene_id from gene")
		rows = self.curs.fetchall()
		for row in rows:
			self.gene_no2gene_id[row[0]] = row[1]
		
		#setup self.gene_prediction_dict
		self.curs.execute("select gene_no, go_no, is_correct, cluster_context from %s"%self.gene_table)
		rows = self.curs.fetchall()
		for row in rows:
			gene_no = row[0]
			go_no = row[1]
			is_correct = row[2]
			cluster_context = row[3]
			context_dict = self.cluster_context2dict(cluster_context)
			item = function_struc()
			item.is_correct = is_correct
			item.context_dict = context_dict
			if gene_no not in self.gene_prediction_dict:
				self.gene_prediction_dict[gene_no] = gene_prediction()
				self.gene_prediction_dict[gene_no].p_functions_struc_dict[go_no] = item
			else:
				self.gene_prediction_dict[gene_no].p_functions_struc_dict[go_no] = item
		
		#setup self.no_of_predicted_genes
		self.no_of_predicted_genes = len(self.gene_prediction_dict)
		
		#get the non-obsolete biological_process GO DAG		
		self.curs.execute("select count(go_no) from go")
		rows = self.curs.fetchall()
		self.no_of_functions = rows[0][0]
		self.curs.execute("select t2t.term1_id, t2t.term2_id, t1.acc, t2.acc from \
			go.term2term t2t, go.term t1, go.term t2 where t2t.term1_id=t1.id and \
			t2t.term2_id=t2.id and t1.is_obsolete=0 and t2.is_obsolete=0 and \
			t1.term_type='biological_process' and t2.term_type='biological_process' ")
		rows = self.curs.fetchall()
		for row in rows:
		#setup the go_graph structure
			self.go_graph.add_edge(row[2], row[3])
			self.go_graph.add_edge(row[3], row[2])
	
		sys.stderr.write("Done\n")

	def cluster_context2dict(self, cluster_context):
		context_dict = {}
		complex_gene_list = cluster_context.split(';')
		for complex_gene in complex_gene_list:
			gene_no, support = map(int, complex_gene.split('/'))
			context_dict[gene_no] = support
		return context_dict
	
	def depth_of_one_node(self, go_no):
		go_id = self.go_no2go_id[go_no]
		return len(GraphAlgo.shortest_path(self.go_graph, 'GO:0008150', go_id))

	def run(self):
		if self.stat_table_fname != 'null':
			self.stat_table_f = csv.writer(open(self.stat_table_fname, 'w'), delimiter='\t')

		for gene_no in self.gene_prediction_dict:
			unit = self.gene_prediction_dict[gene_no]

			if len(unit.p_functions_struc_dict) > 1:
				self.list_no_of_functions_each_gene.append(len(unit.p_functions_struc_dict))
				#rearrange the context_dict into a go_no:Set structure
				self.go_merge_dict = {}
				go_no_list = []
				for go_no in unit.p_functions_struc_dict:
					#!!! MEANING CHANGE !!!#
					#function_struc.context_dict will be a Set data structure to store context genes
					#function_struc.cluster_array will store the go terms to be merged
					go_no_list.append(go_no)
					self.depth_list.append(self.depth_of_one_node(go_no))
					
					item = function_struc()
					item.context_dict = Set(unit.p_functions_struc_dict[go_no].context_dict.keys())
					self.go_merge_dict[go_no] = item
				max_distance = 0
			
				for i in range(len(go_no_list)):
					for j in range(i+1, len(go_no_list)):
						go_id1 = self.go_no2go_id[go_no_list[i]]
						go_id2 = self.go_no2go_id[go_no_list[j]]
						distance = len(GraphAlgo.shortest_path(self.go_graph, go_id1, go_id2))-1
						if distance > max_distance:
							max_distance = distance
				
				self.max_distance_list.append(max_distance)
				self.distinct_contexts()
				if len(self.go_merge_dict) > 1:
					self.list_no_of_contexts_each_gene.append(len(self.go_merge_dict))
				if self.stat_table_fname != 'null':
					self.table_output(gene_no, self.go_merge_dict)

		self.stat_output()
		self.plot()

	def distinct_contexts(self):
		go_no_list = self.go_merge_dict.keys()
		no_of_gos = len(go_no_list)
		need_recursive = 0
		diff_result = 1
		for i in range(no_of_gos):
			for j in range(i+1, no_of_gos):
				diff_result = self.is_context_diff(self.go_merge_dict[go_no_list[i]].context_dict, self.go_merge_dict[go_no_list[j]].context_dict)
				if diff_result == 0:
					#the merged go term is correct if and only if i and j are both correct
					self.go_merge_dict[go_no_list[i]].is_correct *= self.go_merge_dict[go_no_list[j]].is_correct
					#put j into the list to be merged by i
					self.go_merge_dict[go_no_list[i]].cluster_array.append(go_no_list[j])
					#union_update the context Set of i
					self.go_merge_dict[go_no_list[i]].context_dict |= self.go_merge_dict[go_no_list[j]].context_dict
					#recursively run the function
					need_recursive = 1
					#delete the j entry
					del self.go_merge_dict[go_no_list[j]]
					#stop the loop
					break
			if diff_result == 0:
				break
		if need_recursive:
			#recursively run on the modified go_merge_dict again.
			self.distinct_contexts()


	def is_context_diff(self, set1, set2):
		#defaultly, regard two sets as different
		diff_result = 1
		intersection_set = set1 & set2
		set1_extra = set1 - intersection_set
		set2_extra = set2 - intersection_set
		#two important ratios
		r1 = len(set1_extra)/float(len(set1))
		r2 = len(set2_extra)/float(len(set2))
		#if either of them are below the threshold, they are not different. Need discussion.
		if r1 < self.contrast or r2 < self.contrast:
			diff_result = 0
		return diff_result

	def table_output(self, gene_no, go_merge_dict):
		gene_id = self.gene_no2gene_id[gene_no]

		for go_no in go_merge_dict:
			#row is for output
			row = []
			row.append(gene_id)
			unit = go_merge_dict[go_no]
			go_no_list = unit.cluster_array
			go_no_list.append(go_no)
			
			go_name_list = []
			for term_no in go_no_list:
				go_name_list.append(self.go_no2go_name[term_no])
			
			#row.append('%s'%('/'.join(map(repr,go_no_list) )))
			row.append('%s'%('/'.join(go_name_list)))
			row.append('%d'%unit.is_correct)
			context_gene_id_list = []
			for gene_no in unit.context_dict:
				context_gene_id_list.append(self.gene_no2gene_id[gene_no])
			#context_gene_no_list = list(unit.context_dict)
			#context_gene_no_list.sort()
			#row.append('%s'%('/'.join( map(repr, context_gene_no_list) ) ))
			row.append('%s'%('/'.join(context_gene_id_list)))
			self.stat_table_f.writerow(row)
	
	def stat_output(self):
		no_of_multi_functions_genes = len(self.list_no_of_functions_each_gene)
		no_of_diff_contexts_genes = len(self.list_no_of_contexts_each_gene)
		avg_functions_per_gene = (sum(self.list_no_of_functions_each_gene) + self.no_of_predicted_genes - no_of_multi_functions_genes)/float(self.no_of_predicted_genes)
		avg_contexts_per_gene = (sum(self.list_no_of_contexts_each_gene) + self.no_of_predicted_genes - no_of_diff_contexts_genes)/float(self.no_of_predicted_genes)
		sys.stdout.write('\t%f of predicted genes with multiple functions.\n'%(no_of_multi_functions_genes/float(self.no_of_predicted_genes)))
		sys.stdout.write('\t%f of predicted genes with distinct contexts.\n'%(no_of_diff_contexts_genes/float(self.no_of_predicted_genes)))
		sys.stdout.write('\taverage functions per gene: %f.\n'%(avg_functions_per_gene))
		sys.stdout.write('\taverage contexts per gene: %f.\n'%(avg_contexts_per_gene))

	
	def plot(self):
		max_depth = max(self.depth_list)
		max_distance = max(self.max_distance_list)
		r.pdf('depth_hist.pdf')
		r.hist(self.depth_list, breaks=range(max_depth+1), las=1, main='depth histogram', xlab='depth of predicted function' )
		r.dev_off()
		
		r.pdf('distance_hist.pdf')
		r.hist(self.distance_hist, breaks=range(max_distance+1), las =1, main='max distance histogram', xlab='max distance of the predicted functions of one gene')
		r.dev_off()
		
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema=", "table=", "mcl_table=", \
		"contrast=", "report", "commit", "gene_table="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:m:s:rcg:", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = 'cluster_stat'
	mcl_table = 'mcl_result'
	contrast = 0.3
	report = 0
	commit = 0
	gene_table = 'p_gene'
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
		elif opt in ("-s", "--contrast"):
			contrast = float(arg)
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-g", "--gene_table"):
			gene_table = arg
	if len(args) == 1:
		stat_table_fname = args[0]
	else:
		stat_table_fname = 'null'
			
	if schema and gene_table:
		instance = context_specific(hostname, dbname, schema, table, mcl_table, contrast, report, \
			commit, gene_table, stat_table_fname)
		instance.dstruc_loadin()
		instance.run()
	else:
		print __doc__
		sys.exit(2)
