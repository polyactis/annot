#!/usr/bin/env python
"""
Usage: gene_stat.py -k SCHEMA -p P_VALUE_CUT_OFF [OPTION] [STAT_TABLE_FILE]

Option:
	STAT_TABLE_FILE is the file to store the function table of all predicted genes.
		If it's not given, no table will be outputed.
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	cluster_stat(default)
	-m ..., --mcl_table=...	mcl_result(default), mcl_result table corresponding to above table.
	-g ..., --gene_table=...	table to store the stat results, p_gene(default), needed if commit
	-p ..., --p_value_cut_off=...	p_value_cut_off
	-u ..., --unknown_cut_off=...	unknown_cut_off, 1(default), for wu's
	-n ..., --connectivity_cut_off=...	0.8(default), minimum connectivity of a mcl cluster
	-y ..., --recurrence_cut_off=...	5(default), minimum recurrences
	-x ..., --cluster_size_cut_off=...	1000(default), maximum cluster size
	-e ..., --depth_cut_off=...	the minimum depth for a go node to be valid, 3(default)
	-f ..., --dir_files=...	the directory containing all the files outputed by cluster_stat.py
	-j ..., --judger_type=...	how to judge predicted functions, 0(default), 1, 2
	-l, --leave_one_out	use the leave_one_out stat method, default is no leave_one_out
	-w, --wu	Wu's strategy(Default is Jasmine's strategy)
	-r, --report	report the progress(a number)
	-c, --commit	commit the database transaction, records in table gene.
	-v, --dominant	Only assign the dominant function(s) to a gene.
	-o, --log	enable logging
	-s ..., --plottype=...	0, 1, 2, 3(default). 0 is false positive. 1 is true positive. 2 is all. 3 is no plot.
	-q ..., --subgraph_cut_off=...	the cut_off for the subgraph to be valid in one dataset, 0(default)
		NOTICE: 0 means the binary conversion won't be used, just summing the floats.
	-a ..., --accuracy_cut_off=...	the accuracy_cut_off to be based for p_value adjusting, 0(default)
		NOTICE: 0 means using p_value_cut_off method
	-b, --debug	enable debugging, no debug by default
	-h, --help              show this help

Examples:
	gene_stat.py -k sc_yh60_splat_5 -t cluster_stat2 -m mcl_result2 -p 0.001 -l -w p_table_file
	gene_stat.py -k sc_yh60_splat_5 -t cluster_stat2 -m mcl_result2 -p 0.001 -l -w -c -j 2  -g p_gene_cluster_stat2 p_table_file

Description:
	This program is mainly for validation purpose. 
	leave_one_out method must be run after cluster_stat.py and mcl_result_stat.py.
	no leave_one_out method must be run after mcl_result_stat.py.
	So, the former approach requires both --table and --mcl_table.
	The latter only requires --mcl_table.
	
	The unknown_cut_off is the same as the one in gene_stat_on_mcl_result.
	which is unknown genes ratio.
	
	Difference from gene_stat.py, the program will generate some plots(plot 1-6 in Wed-09-22-04) as by-products.
	An important structure p_functions_struc_dict has been added to self.gene_prediction_dict.
	The program now outputs a detailed table about the function prediction for each gene.
	Functions, go_no_accuracy(), table_output() and _table_output() are added for this purpose.
	The column names of the table are stored in variable header of table_output().
	Most modifications only apply to the leave_one_out approach.
	
	If read cluster_stat.py results from files, mcl_table is still needed, but table(cluster_stat) is unnecessary.
	
	Judger-type
	0:	direct match
	1:	L1
	2:	common ancestor's depth is larger than the depth_cut_off

	There're two hidden parameters to control recurrence and connectivity gap_size.
	
	(02-16-05)
		This program can't handle non-leave_one_out.
	(02-16-05) L1 approach,
		predicted function is correct if its jasmine_distance with one of the known functions is 1.
		Jasmine distance is the shorter distance from the two functions to their lowest
		common ancestor.
		Furthermore, the depth of both the predicted and known function should  be >= depth_cut_off.
		Notice, this means, their lowest common ancestor could be on level, depth_cut_off-1.
	(02-16-05)
		the program adjusts p_value in each cell to meet the accuracy_cut_off
		gene_prediction_dict_setup() and p_value_outof_accuracy_cut_off() are added.
		Three judgers are seperated from final. final() becomes a real function, not a dictionary.
	(02-17-05)
		whether to use p_value_cut_off or accuracy_cut_off method depends on whether
		accuracy_cut_off=0 or not.
		In p_value_cut_off method, 
		First, prediction_tuple2accuracy is used to tell the accuracy
		in each cell of the recurrence-connectivity grid. 
		Second, gene_prediction_dict[gene_no].p_functions_dict is used to store a list of distinct
		functions of this gene_no.
	02-19-05
	In gene_distinct_functions_output()
		total predictions = len(self.prediction_tuple2list[tuple]) is totally wrong.
		prediction_tuple2list also contains some predictions that are below that p-value cutoff.
		So remove this column. To see total predictions for each tuple, look up the database
		table.
	02-19-05
	in function submit()
		Changes to table p_gene,
		1. one row means one gene, one cluster, one function. No merging of the clusters.
		2. avg_p_value is the real p-value.
		3. cluster_context and cluster_array loses its meaning. But I kept cluster_array
			because it's easy.
			context_specific.py and subgraph_visualize.py are going to be changed.
		4. p_value_cut_off is the real p_value(same as avg_p_value).
		5. recurrence_cut_off is the real recurrence of the cluster.
		6. connectivity_cut_off is the real connectivity of the cluster.
		7. cluster_size_cut_off is the real size of the cluster.
		8. all the predictions are kept, no any cutoff.
		9. add one field, mcl_id to the end of table p_gene to ease table linking.
"""

import sys, os, psycopg, getopt, csv, fileinput, math
from graphlib import Graph, GraphAlgo
from sets import Set
from rpy import r
from numarray import *

class gene_prediction:
	# class holding prediction information of a gene
	def __init__(self):
		self.tp = {}
		self.tp1 = {}
		self.tn = 0
		self.fp = {}
		self.fn = {}
		#from 02-17-05, p_functions_dict will be used as a list to store distinct functions
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

class gene_stat:
	"""
	dstruc_loadin()
	run()
		--core_from_files()
		or
		--core()
			--_gene_stat_leave_one_out()
				--index_tuple()
				--match()
					--direct_match()
					--L1_match()
					--common_ancestor_deep_enough()
				--add_item_to_dict()
		--gene_prediction_dict_setup()
			--p_value_dictionary_setup()
			--p_value_outof_accuracy_cut_off()
		--final()
			--return_distinct_functions()
		--stat_output()
		--go_no_accuracy()
		--table_output()
			--p_value_cut_off_output()
			--accuracy_cut_off_output()
			--gene_distinct_functions_output()
			--_table_output()
		--submit()
	02-19-05
		three old final_xxx() functions are deleted. also their sub-functions, is_L0 and is_L1.
	
	"""
	def __init__(self, hostname, dbname, schema, table, mcl_table, p_value_cut_off, unknown_cut_off, \
		connectivity_cut_off, recurrence_cut_off, cluster_size_cut_off, leave_one_out, wu, report=0, \
		log=0, judger_type=0, depth_cut_off =3, dir_files=None, needcommit=0, gene_table='p_gene', \
		dominant=0, plottype=3, stat_table_fname='null', subgraph_cut_off=0.8, debug=0, accuracy_cut_off=0):
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.schema = schema
		self.curs.execute("set search_path to %s"%schema)
		self.table = table
		self.mcl_table = mcl_table
		self.p_value_cut_off = float(p_value_cut_off)
		self.unknown_cut_off = float(unknown_cut_off)
		self.connectivity_cut_off = float(connectivity_cut_off)
		self.recurrence_cut_off = float(recurrence_cut_off)
		self.cluster_size_cut_off = int(cluster_size_cut_off)
		self.leave_one_out = int(leave_one_out)
		self.wu = int(wu)
		self.report = int(report)
		
		self.log = int(log)
		self.match_dict = {0: self.direct_match,
			1: self.L1_match,
			2: self.common_ancestor_deep_enough}
		self.match = self.match_dict[judger_type]
		self.judger_type = int(judger_type)
		self.depth_cut_off = int(depth_cut_off)
		self.dir_files = dir_files
		self.needcommit = int(needcommit)
		self.gene_table = gene_table
		self.dominant = int(dominant)
		self.plottype = int(plottype)
		self.stat_table_fname = stat_table_fname
		self.subgraph_cut_off = float(subgraph_cut_off)
		#debugging flag
		self.debug = int(debug)
		self.accuracy_cut_off = float(accuracy_cut_off)
		
		#debug flags in several functions
		self.debug_gene_prediction_dict_setup = 0
		self.debug_L1_match = 0
		self.debug_final = 0
		self.debug_return_distinct_functions = 0
		self.debug_common_ancestor_deep_enough = 0
		#the gap between two recurrences
		self.recurrence_gap_size = 2
		self.connectivity_gap_size = 2
		self.tp = 0.0
		self.tp_m = 0.0
		self.tp1 = 0.0
		self.tp1_m = 0.0
		self.tn = 0.0
		self.fp = 0.0
		self.fp_m =0.0
		self.fn = 0.0
		#mapping between gene_no and its directly assigned go functions
		#key is a gene_no, value is a set of go_id's.
		self.gene_no2direct_go = {}
		#mapping between gene_no and go_no list
		self.known_genes_dict = {}
		self.no_of_records = 0
		if self.log:
			self.log_file = open('/tmp/gene_stat_plot.log','w')
		#the dictionary having (recurrence, connectivity) as a key, [[p_value, cluster_id, gene_no, go_no, correct], [],...] as a value
		self.prediction_tuple2list ={}
		#the dictionary having (recurrence, connectivity) as a key, p_value_cut_off as a value
		self.prediction_tuple2p_value_cut_off = {}
		#the accuracy in each cell
		self.prediction_tuple2accuracy = {}
		
		self.gene_prediction_dict = {}
		self.no_of_p_known = 0
		#GO DAG
		self.go_graph = Graph.Graph()
		#mapping between go_no and go_id
		self.go_no2go_id = {}
		#mapping between go_no and go's name
		self.go_no2go_name = {}
		#mapping between gene_no and gene_id
		self.gene_no2gene_id = {}
		#the overall prediction accuracy of each function
		self.go_no2accuracy = {}
		#mapping between a pair of go_no's and its associated distances
		self.go_no2distance = {}
		#mapping go term id and go_no
		self.go_term_id2go_no = {}
		#mapping each go_no to its depth
		self.go_no2depth = {}
		#mapping each go term id to its depth
		self.go_term_id2depth = {}
		
		#the dictionary holding all clusters
		self.mcl_id2vertex_set = {}
		
		#data structures for plotting, key is go_no and value is a Set of genes of clusters
		self.go_no2gene = {}
		self.go_no2gene_fp = {}
		self.go_no2cluster = {}
		self.go_no2cluster_fp = {}
		self.cluster_size2cluster = {}
		self.cluster_size2cluster_fp = {}
		self.cluster_size2go_no = {}
		self.cluster_size2go_no_fp = {}
		self.dataset_no2cluster = {}
		self.dataset_no2cluster_fp = {}
		self.dataset_no2go_no = {}
		self.dataset_no2go_no_fp = {}
		self.gene_no2cluster = {}
		self.gene_no2cluster_fp = {}
	
	def parameter_reset_and_cleanup(self, unknown_cut_off, connectivity_cut_off, recurrence_cut_off, \
			cluster_size_cut_off, p_value_cut_off, depth_cut_off):
		'''
		this function is for batch_stat.py
		'''
		#parameter_reset
		self.unknown_cut_off = float(unknown_cut_off)
		self.connectivity_cut_off = float(connectivity_cut_off)
		self.recurrence_cut_off = float(recurrence_cut_off)
		self.cluster_size_cut_off = int(cluster_size_cut_off)
		self.p_value_cut_off = float(p_value_cut_off)	
		self.depth_cut_off = int(depth_cut_off)
		
		#cleanup
		self.gene_prediction_dict = {}
		self.no_of_p_known = 0
		self.no_of_records = 0
		self.tp = 0.0
		self.tp_m = 0.0
		self.tp1 = 0.0
		self.tp1_m = 0.0
		self.tn = 0.0
		self.fp = 0.0
		self.fp_m =0.0
		self.fn = 0.0
		#the overall prediction accuracy of each function
		self.go_no2accuracy = {}
		#for database reset
		self.conn.rollback()
		self.curs.execute("set search_path to %s"%self.schema)
		
	def dstruc_loadin(self):
		sys.stderr.write("Loading Data STructure...")
		
		#setup self.known_genes_dict
		self.curs.execute("select gene_no,go_functions from gene where known=TRUE")
		rows = self.curs.fetchall()
		for row in rows:
			go_functions_list = row[1][1:-1].split(',')
			self.known_genes_dict[row[0]] = Set()
			for go_no in go_functions_list:
				self.known_genes_dict[row[0]].add(int(go_no))
		
		#setup self.gene_no2direct_go
		self.curs.execute("select ge.gene_no, a.go_id from graph.association a, gene ge\
			where ge.gene_id=a.gene_id")
		rows = self.curs.fetchall()
		for row in rows:
			if row[0] in self.gene_no2direct_go:
				self.gene_no2direct_go[row[0]].add(row[1])
			else:
				self.gene_no2direct_go[row[0]] = Set([row[1]])
		
		#get the non-obsolete biological_process GO DAG
		self.curs.execute("select t2t.term1_id, t2t.term2_id, t1.acc, t2.acc from \
			go.term2term t2t, go.term t1, go.term t2 where t2t.term1_id=t1.id and \
			t2t.term2_id=t2.id and t1.is_obsolete=0 and t2.is_obsolete=0 and \
			t1.term_type='biological_process' and t2.term_type='biological_process' ")
		rows = self.curs.fetchall()
		for row in rows:
		#setup the go_graph structure
			self.go_graph.add_edge(row[2], row[3])
		
		#setup self.go_no2go_id and self.go_no2depth
		self.curs.execute("select go_no, go_id, depth from go")
		rows = self.curs.fetchall()
		for row in rows:
			self.go_no2go_id[row[0]] = row[1]
			self.go_no2depth[row[0]] = row[2]
		
		#setup self.go_no2go_name and self.go_term_id2go_no
		self.curs.execute("select g.go_no, t.name, t.id from go g, go.term t where g.go_id=t.acc")
		rows = self.curs.fetchall()
		for row in rows:
			self.go_no2go_name[row[0]] = row[1]
			self.go_term_id2go_no[row[2]] = row[0]
		
		#setup self.gene_no2gene_id
		self.curs.execute("select gene_no, gene_id from gene")
		rows = self.curs.fetchall()
		for row in rows:
			self.gene_no2gene_id[row[0]] = row[1]
		
		#setup self.go_no2distance
		if self.judger_type != 0:
			sys.stderr.write("loading distances ....")
			self.curs.execute("DECLARE dist_crs CURSOR FOR select go_id1, go_id2, raw_distance, lee_distance, jasmine_distance, \
				common_ancestor_list from go.node_dist")
			self.curs.execute("fetch 10000 from dist_crs")
			rows = self.curs.fetchall()
			while rows:
				for row in rows:
					go_no1 = self.go_term_id2go_no.get(row[0])
					go_no2 = self.go_term_id2go_no.get(row[1])
					common_ancestor_list = row[5][1:-1].split(',')
					common_ancestor_list = map(int, common_ancestor_list)
					common_ancestor_set = Set(common_ancestor_list)
					if go_no1 and go_no2:
						#key tuple in ascending order
						if go_no1<go_no2:
							self.go_no2distance[(go_no1, go_no2)] = (row[2], row[3], row[4], common_ancestor_set)
						else:
							self.go_no2distance[(go_no2, go_no1)] = (row[2], row[3], row[4], common_ancestor_set)
				self.curs.execute("fetch 5000 from dist_crs")
				rows = self.curs.fetchall()
			sys.stderr.write("done")

		#setup self.no_of_functions
		if self.wu:
			self.curs.execute("select count(go_no) from go")
		else:
			self.curs.execute("select count(go_no) from go where go_no!=0")
		rows = self.curs.fetchall()
		self.no_of_functions = rows[0][0]
		
		#setup self.go_term_id2depth
		self.curs.execute("select id, depth from go.term where depth NOTNULL")
		rows = self.curs.fetchall()
		for row in rows:
			self.go_term_id2depth[row[0]] = row[1]
		
		sys.stderr.write("Done\n")

	def core_from_files(self):
		#following codes are attaching directory path to each file in the list
		file_list = os.listdir(self.dir_files)
		file_path_list = []
		for filename in file_list:
			file_path_list.append(os.path.join(self.dir_files, filename))
		#multiple files constitute the source of data
		self.files = fileinput.input(file_path_list)
		#wrap it with a reader
		self.reader = csv.reader(self.files, delimiter='\t')
		for row in self.reader:
			row[0] = int(row[0])
			row[1] = int(row[1])
			row[3] = float(row[3])
			self.curs.execute("select recurrence_array, vertex_set from %s where mcl_id=%d"%(self.mcl_table, int(row[0])) )
			rows = self.curs.fetchall()
			#first append the recurrence_array
			row.append(rows[0][0])
			#second append the vertex_set
			row.append(rows[0][1])
			#only leave_one_out
			self._gene_stat_leave_one_out(row)

			if self.report and self.no_of_records%2000==0:
				sys.stderr.write('%s%s'%('\x08'*20, self.no_of_records))
		if self.report:
			sys.stderr.write('%s%s'%('\x08'*20, self.no_of_records))

	def core(self):
		#the central function of the class
		if self.leave_one_out:
			#leave_one_out method gets data from both cluster_stat-like and mcl_result-like table
			self.curs.execute("DECLARE crs CURSOR FOR select c.mcl_id, c.leave_one_out, c.p_value_vector, \
				 c.connectivity, m.recurrence_array, m.vertex_set from %s c, %s m where c.mcl_id=m.mcl_id"\
				%(self.table, self.mcl_table))
		else:
			#no leave_one_out method gets data only from mcl_result-like table
			self.curs.execute("DECLARE crs CURSOR FOR select mcl_id, vertex_set, p_value_min, go_no_vector, unknown_gene_ratio, \
				recurrence_array from %s where connectivity>=%f and p_value_min notnull and array_upper(recurrence_array, 1)>=%d\
				and array_upper(vertex_set, 1)<=%d"%(self.mcl_table, self.connectivity_cut_off, self.recurrence_cut_off, self.cluster_size_cut_off))
		
		self.curs.execute("fetch 5000 from crs")
		rows = self.curs.fetchall()
		while rows:
			for row in rows:
				if self.leave_one_out:
					#in leave_one_out, only one gene's function is predicted based on one row
					self._gene_stat_leave_one_out(row)
				else:
					#in no leave_one_out, function of all vertices in that cluster is predicted based on one row
					self._gene_stat_no_leave_one_out(row)
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, self.no_of_records))
			
			self.curs.execute("fetch 5000 from crs")
			rows = self.curs.fetchall()

	def run(self):
		if self.dir_files and self.leave_one_out==0:
			sys.stderr.write("working on files of cluster_stat.py results, it must be leave_one_out.\n")
			sys.exit(2)
		if self.dir_files:
			self.core_from_files()
		else:
			self.core()
		self.gene_prediction_dict_setup()
		#before self.final(), self.gene_prediction_dict has been filled.
		self.final()
		
		self.stat_output()
		#compute the expected accuracy for each specific function
		#needed by table_output() and submit()
		self.go_no_accuracy()
		
		if self.leave_one_out and self.stat_table_fname != 'null':
			self.table_output()

		if self.needcommit and self.leave_one_out:
			#Database updating is too slow. Do it only if needcommit.
			self.submit()

		if self.plottype != 3:
			self.hist_plot(self.go_no2cluster, 'go_no2cluster.png', 'go_no', 'number of clusters')
			self.hist_plot(self.go_no2gene, 'go_no2gene.png', 'go_no', 'number of genes')
			self.hist_plot(self.cluster_size2cluster, 'cluster_size2cluster.png', 'cluster_size', 'number of clusters')
			self.hist_plot(self.cluster_size2go_no, 'cluster_size2go_no.png', 'cluster_size', 'number of go_nos')
			self.hist_plot(self.dataset_no2cluster, 'dataset_no2cluster.png', 'dataset_no', 'number of clusters')
			self.hist_plot(self.dataset_no2go_no, 'dataset_no2go_no.png', 'dataset_no', 'number of go_nos')
			self.hist_plot(self.gene_no2cluster, 'gene_no2cluster.png', 'gene_no', 'number of clusters')
		
			if self.plottype == 1:
				self.hist_plot_ratio(self.go_no2cluster, self.go_no2cluster_fp, 'go_no2cluster_ratio.png', 'go_no', 'number of clusters(ratio)')
				self.hist_plot_ratio(self.go_no2gene, self.go_no2gene_fp, 'go_no2gene_ratio.png', 'go_no', 'number of genes(ratio)')
				self.hist_plot_ratio(self.cluster_size2cluster, self.cluster_size2cluster_fp, 'cluster_size2cluster_ratio.png', 'cluster_size', 'number of clusters(ratio)')
				self.hist_plot_ratio(self.cluster_size2go_no, self.cluster_size2go_no_fp, 'cluster_size2go_no_ratio.png', 'cluster_size', 'number of go_nos(ratio)')
				self.hist_plot_ratio(self.dataset_no2cluster, self.dataset_no2cluster_fp, 'dataset_no2cluster_ratio.png', 'dataset_no', 'number of clusters(ratio)')
				self.hist_plot_ratio(self.dataset_no2go_no, self.dataset_no2go_no_fp, 'dataset_no2go_no_ratio.png', 'dataset_no', 'number of go_nos(ratio)')
				self.hist_plot_ratio(self.gene_no2cluster, self.gene_no2cluster_fp, 'gene_no2cluster_ratio.png', 'gene_no', 'number of clusters(ratio)')

	def index_tuple(self, list):
		new_list = []
		for i in range(len(list)):
			#value is position 0, and index is position 1
			new_list.append((float(list[i]), i))
		#the sort is based on position 0
		new_list.sort()
		return new_list

	def _gene_stat_leave_one_out(self, row):
		"""
		03-08-05
			set a default(1.0) for min_p_value
			fix a bug, alter >self.depth_cut_off to >= self.depth_cut_off
		"""
		mcl_id = row[0]
		gene_no = row[1]
		p_value_vector = row[2][1:-1].split(',')
		connectivity = float(row[3])
		recurrence_array = row[4][1:-1].split(',')
		recurrence_array = map(float, recurrence_array)
		vertex_set = row[5][1:-1].split(',')
		vertex_set = map(int, vertex_set)
		
		#setup mcl_id2vertex_set
		if mcl_id not in self.mcl_id2vertex_set:
			self.mcl_id2vertex_set[mcl_id] = vertex_set
		
		#take the floor of the recurrence
		if self.subgraph_cut_off!=0:
			#0 means no cutoff
			recurrence_array = greater_equal(recurrence_array, self.subgraph_cut_off)
		recurrence = int(math.floor(sum(recurrence_array)/self.recurrence_gap_size)*self.recurrence_gap_size)
		#take the floor of the connectivity *10
		connectivity = int(math.floor(connectivity*10/self.connectivity_gap_size)*self.connectivity_gap_size)
		#setup in prediction_tuple2list
		prediction_tuple = (recurrence, connectivity)
		if prediction_tuple not in self.prediction_tuple2list:
			self.prediction_tuple2list[prediction_tuple] = []
		
		'''
		#cutoff cancelled
		if len(vertex_set) > self.cluster_size_cut_off:
			#the cluster is too big
			return
		'''
		"""
		#we don't do recurrence_cut_off and connectivity_cut_off anymore
		if recurrence < self.recurrence_cut_off or recurrence > (self.recurrence_cut_off+2):
			#the recurrence is not enough
			return
		if connectivity < self.connectivity_cut_off:
			#or connectivity > (self.connectivity_cut_off+0.1):
			#the cluster is not dense enough
			return
		"""
		#default min_p_value
		min_p_value = 1.0
		
		#transform into float type
		p_value_index_tuple_list = self.index_tuple(p_value_vector)
		for (p_value, index) in p_value_index_tuple_list:
			if self.wu:
				#index 0 corresponds to go_no 0.
				go_no = index
			else:
				#index 0 corresponds to go_no 1
				go_no = index+1
			if self.go_no2depth[go_no] >= self.depth_cut_off:
				min_p_value = p_value
				break

		#1.0 means the corresponding function has no associated genes in the cluster.
		if min_p_value >= 1.0:
			return
			
		if self.wu:
			if float(p_value_vector[0]) > self.unknown_cut_off:
				#too many unknown genes, and now cut_off is ratio.
				return
		"""
		#p_value_cut_off is flexible
		if min_p_value > self.p_value_cut_off:
			#none of the predicted functions in this cluster is significant
			return
		"""
		#The cluster is an eligible cluster. Passing all the cut_offs.
		#
		self.no_of_records += 1
		"""
		#gene_prediction_dict will be dealt with in gene_prediction_dict_setup()
		
		if gene_no not in self.gene_prediction_dict:
			item = gene_prediction()
			self.gene_prediction_dict[gene_no] = item
		self.gene_prediction_dict[gene_no].mcl_id_list.append(mcl_id)	
		"""
		for (p_value, index) in p_value_index_tuple_list:
			if p_value > min_p_value:
				break
			elif index == 0 or p_value==1.0:
				#0 is the unknown function, this is almost impossible because its depth = 2(see condition above)
				#1.0 is for function that has no associated genes
				continue
			elif p_value == min_p_value:
				if self.wu:
					#index 0 corresponds to go_no 0.
					go_no = index
				else:
					#index 0 corresponds to go_no 1
					go_no = index+1
				
				if gene_no in self.known_genes_dict:
					k_functions_set = self.known_genes_dict[gene_no]
					is_correct = self.match(go_no, k_functions_set)
				else:
					#unknown gene
					is_correct = -1
				
				prediction_list = [p_value, mcl_id, gene_no, go_no, is_correct]
				self.prediction_tuple2list[prediction_tuple].append(prediction_list)
				
				"""					
				#gene_prediction_dict will be dealt with in gene_prediction_dict_setup()
				
				#tp is indicator variable indicating the prediction is true positive or false.
				tp = 0
	
				if go_no not in self.gene_prediction_dict[gene_no].p_functions_dict:
					#value in p_functions_dict stores the number of associated clusters.
					self.gene_prediction_dict[gene_no].p_functions_dict[go_no] = 1
				else:
					self.gene_prediction_dict[gene_no].p_functions_dict[go_no] += 1
					
				#setup the p_functions_struc_dict
				#
				if go_no not in self.gene_prediction_dict[gene_no].p_functions_struc_dict:
					self.gene_prediction_dict[gene_no].p_functions_struc_dict[go_no] = function_struc()

				#push in the min_p_value
				self.gene_prediction_dict[gene_no].p_functions_struc_dict[go_no].p_value_list.append(min_p_value)
				#push in the mcl_id
				self.gene_prediction_dict[gene_no].p_functions_struc_dict[go_no].cluster_array.append(mcl_id)
				#push in the surrounding vertices of gene_no
				for vertex in vertex_set:
					if vertex != gene_no:
						if vertex not in self.gene_prediction_dict[gene_no].p_functions_struc_dict[go_no].context_dict:
							self.gene_prediction_dict[gene_no].p_functions_struc_dict[go_no].context_dict[vertex] = 1
						else:
							self.gene_prediction_dict[gene_no].p_functions_struc_dict[go_no].context_dict[vertex] += 1
				#
				#p_functions_struc_dict setup done
				"""
		
		#Fill in the dictionaries for plotting.
		#One possible error is the case that one cluster has >1 functions with min_p_value.
		#I didn't fix it because the possibility is low and the plotting dictionaries are not crucial.
		#
		if self.plottype != 3:
			if (self.plottype == 2) or (tp == self.plottype):
				#codes below fill in the data structures for plotting
				#variable for cluster_size2cluster and cluster_size2go_no
				cluster_size = len(vertex_set)
				self.add_item_to_dict(self.cluster_size2cluster, cluster_size, mcl_id)
				self.add_item_to_dict(self.cluster_size2go_no, cluster_size, go_no)
				self.add_item_to_dict(self.go_no2cluster, go_no, mcl_id)
				self.add_item_to_dict(self.go_no2gene, go_no, gene_no)
				for dataset_no in recurrence_array:
					self.add_item_to_dict(self.dataset_no2cluster, dataset_no, mcl_id)
					self.add_item_to_dict(self.dataset_no2go_no, dataset_no, go_no)
				self.add_item_to_dict(self.gene_no2cluster, gene_no, mcl_id)
			if tp == 0:
				cluster_size = len(vertex_set)
				self.add_item_to_dict(self.cluster_size2cluster_fp, cluster_size, mcl_id)
				self.add_item_to_dict(self.cluster_size2go_no_fp, cluster_size, go_no)
				self.add_item_to_dict(self.go_no2cluster_fp, go_no, mcl_id)
				self.add_item_to_dict(self.go_no2gene_fp, go_no, gene_no)
				for dataset_no in recurrence_array:
					self.add_item_to_dict(self.dataset_no2cluster_fp, dataset_no, mcl_id)
					self.add_item_to_dict(self.dataset_no2go_no_fp, dataset_no, go_no)
				self.add_item_to_dict(self.gene_no2cluster_fp, gene_no, mcl_id)
		#
		#Filling done.
		
	def _gene_stat_no_leave_one_out(self, row):
		mcl_id = row[0]
		vertex_set = row[1][1:-1].split(',')
		vertex_set = map(int, vertex_set)
		p_value_min = row[2]
		go_no_vector = row[3][1:-1].split(',')
		go_no_vector = map(int, go_no_vector)
		unknown_gene_ratio = row[4]
		if p_value_min>self.p_value_cut_off or unknown_gene_ratio>self.unknown_cut_off:
			return
		self.no_of_records += 1
		#variable for cluster_size2cluster and cluster_size2go_no
		cluster_size = len(vertex_set)
		if cluster_size not in self.cluster_size2cluster:
			self.cluster_size2cluster[cluster_size] = Set([mcl_id])
		else:
			self.cluster_size2cluster[cluster_size].add(mcl_id)
		if cluster_size not in self.cluster_size2go_no:
			self.cluster_size2go_no[cluster_size] = Set()
			#it will be expanded in the last block of this function
		for gene_no in vertex_set:
			if gene_no not in self.gene_prediction_dict:
				item = gene_prediction()
				self.gene_prediction_dict[gene_no] = item
			self.gene_prediction_dict[gene_no].mcl_id_list.append(mcl_id)
			for go_no in go_no_vector:
				#every go_no in go_no_vector is assigned to this gene_no
				if go_no not in self.gene_prediction_dict[gene_no].p_functions_dict:
					self.gene_prediction_dict[gene_no].p_functions_dict[go_no] = 1
				else:
					self.gene_prediction_dict[gene_no].p_functions_dict[go_no] += 1
				#every go_no's associated gene Set is expanded
				if go_no not in self.go_no2gene:
					self.go_no2gene[go_no] = Set([gene_no])
				else:
					self.go_no2gene[go_no].add(gene_no)
		for go_no in go_no_vector:
			#add each go_no into the corresponding cell in self.cluster_size2go_no
			self.cluster_size2go_no[cluster_size].add(go_no)
			#every go_no's associated cluster Set is expanded
			if go_no not in self.go_no2cluster:
				self.go_no2cluster[go_no] = Set([mcl_id])
			else:
				self.go_no2cluster[go_no].add(mcl_id)

	def gene_prediction_dict_setup(self):
		if self.debug_gene_prediction_dict_setup:
			print "\t\t### in function gene_prediction_dict_setup() ###"
		for tuple in self.prediction_tuple2list:
			#following two variables is used to count accuracy in each cell
			no_good_predictions = 0.0
			no_total_predictions = 0.0
			prediction_list = self.prediction_tuple2list[tuple]
			p_value_dictionary = self.p_value_dictionary_setup(prediction_list)
			if self.debug_gene_prediction_dict_setup:
				print "recurrence and connectivity:%s"%repr(tuple)
			if self.accuracy_cut_off>0:
				p_value_cut_off = self.p_value_outof_accuracy_cut_off(p_value_dictionary, self.accuracy_cut_off)
			else:
				#0 means using p_value_cut_off method.
				p_value_cut_off = self.p_value_cut_off
			self.prediction_tuple2p_value_cut_off[tuple] = p_value_cut_off
			for p_value in p_value_dictionary:
				prediction_list = p_value_dictionary[p_value]
				if p_value <=  p_value_cut_off:
					for unit in prediction_list:
						mcl_id = unit[0]
						gene_no = unit[1]
						go_no = unit[2]
						is_correct = unit[3]
						if is_correct != -1:
							#-1 is unknown gene
							no_good_predictions += is_correct
							no_total_predictions += 1
						if self.debug_gene_prediction_dict_setup:
							print "p_value: %s with cutoff :%f "%(p_value, p_value_cut_off)
							print "data is %s"%repr(unit)
							raw_input("pause:")

						if gene_no not in self.gene_prediction_dict:
							item = gene_prediction()
							self.gene_prediction_dict[gene_no] = item
						if go_no not in self.gene_prediction_dict[gene_no].p_functions_struc_dict:
							self.gene_prediction_dict[gene_no].p_functions_struc_dict[go_no] = function_struc()
							self.gene_prediction_dict[gene_no].p_functions_struc_dict[go_no].is_correct = is_correct
						#transfer to a function_struc to ease coding, it's like pointer transfer in C, copy() is used for value transfer.
						item = self.gene_prediction_dict[gene_no].p_functions_struc_dict[go_no]
						#push in the p_value
						item.p_value_list.append(p_value)
						#push in the mcl_id
						item.cluster_array.append(mcl_id)
						#push in the context
						for vertex in self.mcl_id2vertex_set[mcl_id]:
							if vertex not in item.context_dict:
								item.context_dict[vertex] = 1
							else:
								item.context_dict[vertex] += 1
			#push in the accuracy.
			if no_total_predictions != 0: 
				self.prediction_tuple2accuracy[tuple] = no_good_predictions/no_total_predictions
			else:
				
				self.prediction_tuple2accuracy[tuple] = -1
		if self.debug_gene_prediction_dict_setup:
			print "\t\t### leave gene_prediction_dict_setup() ###"
			
	def p_value_dictionary_setup(self, prediction_list):
		p_value_dictionary = {}
		for ls in prediction_list:
			p_value = ls[0]
			if p_value in p_value_dictionary:
				p_value_dictionary[p_value].append( ls[1:])
			else:
				p_value_dictionary[p_value] = [ls[1:]]
		return p_value_dictionary

	def p_value_outof_accuracy_cut_off(self, p_value_dictionary, accuracy_cut_off):
		if self.debug:
			print "\t\t###Enter p_value_outof_accuracy_cut_off()###"
		total_clusters = 0.0
		good_clusters = 0.0
		#if no good p_value_cut_off to get accuracy, it's to be 0
		p_value_cut_off = 0.0
		p_value2cumu_accuracy = {}
		p_value_list = p_value_dictionary.keys()
		p_value_list.sort()
		if self.debug:
			print 'p_value_list sorted:%s'%repr(p_value_list)
			
		for p_value in p_value_list:
			prediction_array = array(p_value_dictionary[p_value])
			#is_correct = 0 or 1 is a prediction of known genes
			no_of_clusters = sum(greater_equal(prediction_array[:,3], 0))
			#is_correct = 1 is a good prediction
			no_of_good_clusters = sum(greater(prediction_array[:,3],0))
			total_clusters += no_of_clusters
			good_clusters += no_of_good_clusters
			
			if total_clusters == 0:
				p_value2cumu_accuracy[p_value] = -1
				if good_clusters !=0:
					sys.stderr.write("ERROR: good_clusters are more than total_clusters \n")
			else:
				p_value2cumu_accuracy[p_value] = good_clusters/total_clusters
			if self.debug:
				print "p_value: %s"%p_value
				print "prediction_array: %s"%prediction_array
				print "no_of_clusters: %s"%no_of_clusters
				print "no_of_good_clusters: %s"%no_of_good_clusters
				print "total_clusters:%s"%total_clusters
				print "good_clusters: %s"%good_clusters
				print "accuracy: %s"%(good_clusters/total_clusters)
				raw_input("pause:")
		p_value_list.reverse()
		for p_value in p_value_list:
			if p_value2cumu_accuracy[p_value] >= accuracy_cut_off:
				p_value_cut_off = p_value
				if self.debug:
					print "accuracy: %s"%p_value2cumu_accuracy[p_value]
				break
		if self.debug:
			print "p_value_cut_off:  %s"%p_value_cut_off
			print "\t\t##leave p_value_outof_accuracy_cut_off()"
		return p_value_cut_off
			
	def direct_match(self, p_go_no, k_functions_set):
		if self.go_no2depth[p_go_no] < self.depth_cut_off:
			#first see if it's deep enough
			return 0
		else:
			if p_go_no in k_functions_set:
				return 1
			else:
				return 0
	
	def L1_match(self, p_go_no, k_functions_set):
		if self.debug_L1_match:
			print "\t\t ### In function L1_match() "
		#default not match
		flag = 0
		if self.go_no2depth[p_go_no] < self.depth_cut_off:
			#not good 
			return 0
		for k_go_no in k_functions_set:
			if self.go_no2depth[k_go_no] < self.depth_cut_off:
				#the known function is above the depth_cut_off, discard it
				continue
			elif k_go_no == p_go_no:
				flag = 1
				break
			elif k_go_no < p_go_no:
				key = (k_go_no, p_go_no)
			elif  k_go_no > p_go_no:
				key = (p_go_no, k_go_no)
			if key in self.go_no2distance:
				if self.go_no2distance[key][2] == 1:
					#jasmine distance = 1
					if self.debug_L1_match:
						print 'One of %s and %s are one step away from their lowest common ancestor, \
							with depth_cut_off, %d'%(p_go_no, k_go_no, self.depth_cut_off)
						raw_input("Pause:")
					flag = 1
					break
			else:
				if self.debug_L1_match:
					print "something wrong, %s's distance is unknown"%repr(key)
			
		if self.debug_L1_match:
			print "\t\t ###leave function L1_match()"
		return flag
	
	def common_ancestor_deep_enough(self, p_go_no, k_functions_set):
		if self.debug_common_ancestor_deep_enough:
			print "\t\t ### Enter common_ancestor_deep_enough() "
		ancestor_set = Set()
		for k_go_no in k_functions_set:
			if k_go_no == p_go_no:
				continue
			elif k_go_no < p_go_no:
				key = (k_go_no, p_go_no)
			elif k_go_no > p_go_no:
				key = (p_go_no, k_go_no)
			if key in self.go_no2distance:
				ancestor_set |= self.go_no2distance[key][3]
			elif self.debug_common_ancestor_deep_enough:
				print "distance for %s doesn't exist.\n"%(repr(key))
		#in case no ancestor at all
		depth = 0
		for ancestor in ancestor_set:
			depth = self.go_term_id2depth[ancestor]
			if depth >= self.depth_cut_off:
				if self.debug_common_ancestor_deep_enough:
					print "%s's common_ancestor %s\n"%(self.go_no2go_id[p_go_no], \
						self.go_no2go_id[self.go_term_id2go_no[ancestor]])
				#pre-stop the loop
				break
		if depth >= self.depth_cut_off:
			return 1
		else:
			return 0
		if self.debug_common_ancestor_deep_enough:
			print "\t\t ### Leave common_ancestor_deep_enough() "

	def final(self):
		'''
		get some global statistics, no_of_p_known, tp, fp, tp_m, fp_m
		'''
		if self.debug_final:
			print "\t\t##Enter final()"
		for gene_no in self.gene_prediction_dict:
			item = self.gene_prediction_dict[gene_no].p_functions_struc_dict
			if gene_no in self.known_genes_dict:
				self.no_of_p_known += 1
			
			#not a dictionary anymore, a list instead.
			#get distinct functions for this gene_no
			self.gene_prediction_dict[gene_no].p_functions_dict =\
					self.return_distinct_functions(item.keys())
			for go_no in item:
				is_correct = item[go_no].is_correct
				no_of_clusters = len(item[go_no].cluster_array)
				if is_correct!=-1:
					#-1 is unknown gene
					self.tp += is_correct
					self.fp += -(is_correct-1)
					self.tp_m += is_correct*no_of_clusters
					self.fp_m += -(is_correct-1)*no_of_clusters
		if self.debug_final:
			print "\t\t##leave final()"
	
	def return_distinct_functions(self, go_no_list):
		'''
		collapse those parent-child nodes, which one is kept is random
		'''
		if self.debug_return_distinct_functions:
			print "\t\t##Enter return_distinct_functions()"
			print "original go_no_list : %s"%repr(go_no_list)
		go_no_list_to_return = []
		#a dictionary to indicate the removal status of a go_no
		go_no_list2flag = {}
		for go_no in go_no_list:
			go_no_list2flag[go_no] = 0
		for i in range(len(go_no_list)):
			go_no = go_no_list[i]
			if go_no_list2flag[go_no] ==0 :
				#not flagged,
				go_no_list_to_return.append(go_no)
				if self.debug_return_distinct_functions:
					print "go_no %s added"%go_no
				for j in range(i+1, len(go_no_list)):
					go_no2 = go_no_list[j]
					if go_no < go_no2:
						key= (go_no, go_no2)
					else:
						key = (go_no2, go_no)
					if self.go_no2distance[key][2] == 0:
						#jasmine_distance=0 means they are parent-child
						go_no_list2flag[go_no2] = 1
						if self.debug_return_distinct_functions:
							print "go_no %s flagged to remove"%go_no2
		if self.debug_return_distinct_functions:
			print "\t\t##leave return_distinct_functions()"
		return go_no_list_to_return


	def list_stringlist(self, list):
		return '{' + repr(list)[1:-1] + '}'
	
	def submit(self):
		"""
		02-19-05
			Changes to table p_gene,
			1. one row means one gene, one cluster, one function. No merging of the clusters.
			2. avg_p_value is the real p-value.
			3. cluster_context and cluster_array loses its meaning. But I kept cluster_array
				because it's easy.
				context_specific.py and subgraph_visualize.py are going to be changed.
			4. p_value_cut_off is the real p_value(same as avg_p_value).
			5. recurrence_cut_off is the real recurrence of the cluster.
			6. connectivity_cut_off is the real connectivity of the cluster.
			7. cluster_size_cut_off is the real size of the cluster.
			8. all the predictions are kept, no any cutoff.
			9. add one field, mcl_id to the end of table p_gene to ease table linking.
			
		"""
		sys.stderr.write("Database transacting...")
		if self.gene_table!='p_gene':
			#create the table if it's not 'p_gene'
			self.curs.execute("create table %s(\
				p_gene_id       serial,\
				gene_no integer,\
				go_no   integer,\
				is_correct      integer,\
				avg_p_value     float,\
				e_accuracy      float,\
				no_of_clusters  integer,\
				cluster_context varchar,\
				cluster_array   integer[],\
				p_value_cut_off float,\
				recurrence_cut_off      float,\
				connectivity_cut_off    float,\
				cluster_size_cut_off    integer,\
				unknown_cut_off      float,\
				depth_cut_off integer,\
				mcl_id integer\
				)"%self.gene_table)
		
		"""the value of self.prediction_tuple2list, [[p_value, cluster_id, gene_no, go_no, correct], ... ] """
		for (tuple, prediction_list)  in self.prediction_tuple2list.iteritems():
			recurrence = tuple[0]
			connectivity = tuple[1]
			for unit in prediction_list:
				p_value = unit[0]
				mcl_id = unit[1]
				gene_no = unit[2]
				go_no = unit[3]
				is_correct = unit[4]
				if go_no not in self.go_no2accuracy:
					#There's possibility that one function can not be validated.
					#Because after masking a known gene to be unknown, the representative function of a cluster might be changed.
					#
					e_accuracy = -1.0
				else:
					e_accuracy = self.go_no2accuracy[go_no].ratio
				self.curs.execute("insert into %s(gene_no, go_no, is_correct, avg_p_value, e_accuracy,\
					no_of_clusters, cluster_array, p_value_cut_off, recurrence_cut_off,\
					connectivity_cut_off, cluster_size_cut_off, unknown_cut_off, depth_cut_off, mcl_id)\
					values(%d, %d, %d, %f, %f, %s, ARRAY%s, %f, %s, %s, %s, %s, %s, %s)"%\
					(self.gene_table, gene_no, go_no, is_correct, p_value, e_accuracy,\
					1, repr([mcl_id]), p_value, recurrence,\
					connectivity, len(self.mcl_id2vertex_set[mcl_id]), self.unknown_cut_off, self.depth_cut_off, mcl_id))
		"""
		#replaced by codes above, see function documentation
		for gene_no in self.gene_prediction_dict:
			unit = self.gene_prediction_dict[gene_no].p_functions_struc_dict
			for go_no in unit:
				no_of_clusters = len(unit[go_no].cluster_array)
				#convert the context_dict into a human-readable list.
				context_list = []
				context_dict_keys = unit[go_no].context_dict.keys()
				context_dict_keys.sort()
				for member in context_dict_keys:
					no_of_containing_clusters = unit[go_no].context_dict[member]
					context_list.append('%s/%d'%(member, no_of_containing_clusters ))
				if go_no not in self.go_no2accuracy:
					#There's possibility that one function can not be validated.
					#Because after masking a known gene to be unknown, the representative function of a cluster might be changed.
					#
					e_accuracy = -1.0
				else:
					e_accuracy = self.go_no2accuracy[go_no].ratio
				self.curs.execute("insert into %s(gene_no, go_no, is_correct, avg_p_value, e_accuracy,\
					no_of_clusters, cluster_context, cluster_array, p_value_cut_off, recurrence_cut_off,\
					connectivity_cut_off, cluster_size_cut_off, unknown_cut_off, depth_cut_off)\
					values(%d, %d, %d, %f, %f, %d, '%s', ARRAY%s, %f, %s, %f, %d, %f, %d)"%\
					(self.gene_table, gene_no, go_no, unit[go_no].is_correct, sum(unit[go_no].p_value_list)/no_of_clusters,\
					e_accuracy, no_of_clusters, ';'.join(context_list),\
					repr(unit[go_no].cluster_array), self.p_value_cut_off, self.recurrence_cut_off,\
					self.connectivity_cut_off, self.cluster_size_cut_off, self.unknown_cut_off, self.depth_cut_off))
		"""
		if self.needcommit:
			self.curs.execute("end")
		sys.stderr.write("done.\n")

	def stat_output(self):
		"""
		(02-17-05)
			redirect the output to self.stat_table_fname.
			table_output() open the same file with 'a' flag.
		"""
		outf = open(self.stat_table_fname, 'w')
		outf.write('\n\tp_value_cut_off:%f unknown_cut_off:%f connectivity_cut_off:%f\n'%(self.p_value_cut_off, self.unknown_cut_off, self.connectivity_cut_off))
		outf.write('\trecurrence_cut_off:%s cluster_size_cut_off:%d\n'%(self.recurrence_cut_off, self.cluster_size_cut_off))
		outf.write('\tdepth_cut_off:%d\n'%(self.depth_cut_off))
		outf.write('\tTotal genes: %d\n'%len(self.gene_prediction_dict))
		outf.write('\tTotal known genes: %d\n'%self.no_of_p_known)
		outf.write("\tBased on functions:\n")
		outf.write('\t\tTP0: %d  TP1: %d  TN: %d  FP: %d  FN: %d\n'%(self.tp, self.tp1, self.tn, self.fp, self.fn))
		if (self.tp+self.tp1+self.fn) == 0:
			outf.write('\t\tSensitvity: Null\n')
		else:
			outf.write('\t\tSensitvity: %f\n'%((self.tp+self.tp1)/(self.tp+self.tp1+self.fn)))
		if (self.fp+self.tn) == 0:
			outf.write('\t\tSpecificity: Null\n')
		else:
			outf.write('\t\tSpecificity: %f\n'%(self.tn/(self.fp+self.tn)))
		if (self.tp+self.tp1+self.fp) == 0:
			outf.write('\t\tFalse Positive Ratio: Null\n')
		else:
			outf.write('\t\tFalse Positive Ratio: %f\n'%(self.fp/(self.tp+self.tp1+self.fp)))
		outf.write("\tBased on clusters:\n")
		outf.write('\t\tTP0_M: %d  TP1_M: %d  FP_M: %d\n'%(self.tp_m, self.tp1_m, self.fp_m))
		if (self.tp_m+self.tp1_m+self.fp_m) == 0:
			outf.write('\t\tFalse Positive Ratio: Null\n')
		else:
			outf.write('\t\tFalse Positive Ratio: %f\n'%(self.fp_m/(self.tp_m+self.tp1_m+self.fp_m)))
		outf.close()

	def hist_plot(self, dict, filename, xlabel, ylabel):
		#convert self.go_no2cluster and self.go_no2gene into histograms
		r.png('%s'%filename)
		x_list = []
		y_list = []
		for (key, value) in dict.iteritems():
			x_list.append(key)
			y_list.append(len(value))
		r.plot(x_list, y_list, type='h', xlab=xlabel, ylab=ylabel, main='%s v.s. %s'%(ylabel, xlabel))
		r.dev_off()

	def hist_plot_ratio(self, dict1, dict2, filename, xlabel, ylabel):
		#convert self.go_no2cluster and self.go_no2gene into histograms
		r.png('%s'%filename)
		x_list = []
		y_list = []
		keys = Set(dict1.keys()).union( Set(dict2.keys()) )
		for key in keys:
			value1 = dict1.get(key, [])
			value2 = dict2.get(key, [])
			ratio = float(len(value1))/(len(value1)+len(value2))
			x_list.append(key)
			y_list.append(ratio)
		r.plot(x_list, y_list, type='h', xlab=xlabel, ylab=ylabel, main='%s v.s. %s'%(ylabel, xlabel))
		r.dev_off()
		
	def add_item_to_dict(self, dict, key, value):
		if key not in dict:
			dict[key] = Set([value])
		else:
			dict[key].add(value)
	
	def go_no_accuracy(self):
		#this function computes the prediction accuracy for each go function
		for gene_no in self.gene_prediction_dict:
			if gene_no in self.known_genes_dict:
				unit = self.gene_prediction_dict[gene_no].p_functions_struc_dict
				for go_no in unit:
					if go_no not in self.go_no2accuracy:
						self.go_no2accuracy[go_no] = accuracy_struc()
					#the counter 'all' always gets increased.
					self.go_no2accuracy[go_no].all += 1
					if unit[go_no].is_correct == 1:
						#the counter 'correct' only gets increased if the prediction is correct
						self.go_no2accuracy[go_no].correct += 1
		#finally, compute all the ratios.
		for go_no in self.go_no2accuracy:
			unit = self.go_no2accuracy[go_no]
			self.go_no2accuracy[go_no].ratio = unit.correct/float(unit.all)
		
	def table_output(self):
		'''
		--p_value_cut_off_output
		--accuracy_cut_off_output
		--gene_distinct_functions_output
		'''
		#only for leave_one_out approach and when stat_table_fname is not 'null'
		#mode is append, because stat_output() output first.(02-17-05)
		stat_table_f = csv.writer(open(self.stat_table_fname, 'a'), delimiter='\t')
		
		self.p_value_cut_off_output(stat_table_f)
		self.accuracy_cut_off_output(stat_table_f)
		self.gene_distinct_functions_output(stat_table_f)
		header = ['gene_id', 'function_known', 'function_predicted', 'is_correct', 'average p_value', 'expected accuracy',\
		'#supporting clusters', 'mcl_id_list', 'cluster_context']
		stat_table_f.writerow(header)
		for gene_no in self.gene_prediction_dict:
			#first output only the known genes
			if gene_no in self.known_genes_dict:
				unit = self.gene_prediction_dict[gene_no].p_functions_struc_dict
				self._table_output(stat_table_f, gene_no, unit)
		
		for gene_no in self.gene_prediction_dict:
			#second, output the unknown genes
			if gene_no not in self.known_genes_dict:
				unit = self.gene_prediction_dict[gene_no].p_functions_struc_dict
				self._table_output(stat_table_f, gene_no, unit)
		
		del stat_table_f

	def _table_output(self, f_handler, gene_no, unit):
		#indicator for the number of functions associated with this gene
		for go_no in unit:
			row = []
			row.append(self.gene_no2gene_id[gene_no])
			#function_known might not be available
			if gene_no in self.known_genes_dict:
				function_known_list = []
				for go_known in self.known_genes_dict[gene_no]:
					go_id = self.go_no2go_id[go_known]
					if go_id not in self.gene_no2direct_go[gene_no]:
						function_known_list.append('%s(parent)'%go_id)
					else:
						function_known_list.append('%s'%go_id)
				row.append(';'.join(function_known_list))
			else:
				row.append('')

			row.append('%s(%s)'%(self.go_no2go_name[go_no], self.go_no2go_id[go_no]))
			row.append(unit[go_no].is_correct)
			no_of_clusters = len(unit[go_no].p_value_list)
			row.append(sum(unit[go_no].p_value_list)/float(no_of_clusters))
			if go_no in self.go_no2accuracy:
				#There's possibility that one function can not be validated.
				#Because after masking a known gene to be unknown, the representative function of a cluster might be changed.
				#
				row.append('%2.2f%%'%(self.go_no2accuracy[go_no].ratio*100))
			else:
				row.append('NULL')
			row.append(no_of_clusters)
			row.append(repr(unit[go_no].cluster_array))
			#convert the context_dict into a human-readable list.
			context_list = []
			context_dict_keys = unit[go_no].context_dict.keys()
			context_dict_keys.sort()
			for member in context_dict_keys:
				no_of_containing_clusters = unit[go_no].context_dict[member]
				context_list.append('%s/%d'%(self.gene_no2gene_id[member], no_of_containing_clusters ))
			row.append(';'.join(context_list))
			#write the row into the stat_table_fname
			f_handler.writerow(row)

	def p_value_cut_off_output(self, f_handler):
		f_handler.writerow(["recurrence", "connectivity", "p_value"])
		for tuple in self.prediction_tuple2p_value_cut_off:
			f_handler.writerow([tuple[0],tuple[1],self.prediction_tuple2p_value_cut_off[tuple]])

	def accuracy_cut_off_output(self, f_handler):
		"""
		02-19-05
			total predictions = len(self.prediction_tuple2list[tuple]) is totally wrong.
			prediction_tuple2list also contains some predictions that are below that p-value cutoff.
			So remove this column. To see total predictions for each tuple, look up the database
			table.
		"""
		f_handler.writerow(['recurrence', 'connectivity', 'accuracy'])
		for tuple in self.prediction_tuple2accuracy:
			f_handler.writerow([tuple[0],tuple[1],self.prediction_tuple2accuracy[tuple]])
	
	def gene_distinct_functions_output(self, f_handler):
		f_handler.writerow(['gene_id', 'go_no_list', '#'])
		for gene_no in self.gene_prediction_dict:
			go_no_list = self.gene_prediction_dict[gene_no].p_functions_dict
			f_handler.writerow([self.gene_no2gene_id[gene_no], repr(go_no_list), len(go_no_list)])
	
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema=", "table=", "mcl_table=", "p_value_cut_off=",\
		"unknown_cut_off=", "connectivity_cut_off=", "recurrence_cut_off=", "cluster_size_cut_off=", "depth_cut_off=",\
		"log", "judger_type=", "dir_files=", "leave_one_out", "wu", "report", "commit", "gene_table=", "dominant", \
		"plottype=", "subgraph_cut_off=", "debug", "accuracy_cut_off="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:m:p:u:n:y:x:e:oj:f:lwrcg:vs:q:ba:", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = 'cluster_stat'
	mcl_table = 'mcl_result'
	p_value_cut_off = None
	connectivity_cut_off = 0.8
	recurrence_cut_off = 5
	cluster_size_cut_off = 1000
	depth_cut_off = 3
	dir_files = None
	log = 0
	judger_type = 0
	leave_one_out = 0
	wu = 0
	report = 0
	commit = 0
	unknown_cut_off = 1
	gene_table = 'p_gene'
	dominant = 0
	plottype = 3
	subgraph_cut_off = 0
	debug = 0
	accuracy_cut_off = 0
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
		elif opt in ("-u", "--unknown_cut_off"):
			unknown_cut_off = float(arg)
		elif opt in ("-n", "--connectivity_cut_off"):
			connectivity_cut_off = float(arg)
		elif opt in ("-y", "--recurrence_cut_off"):
			recurrence_cut_off = float(arg)
		elif opt in ("-x", "--cluster_size_cut_off"):
			cluster_size_cut_off = int(arg)
		elif opt in ("-e", "--depth_cut_off"):
			depth_cut_off = int(arg)
		elif opt in ("-f", "--dir_files"):
			dir_files = arg
		elif opt in ("-o", "--log"):
			log = 1
		elif opt in ("-j", "--judger_type"):
			judger_type = int(arg)
		elif opt in ("-l", "--leave_one_out"):
			leave_one_out = 1
		elif opt in ("-w", "--wu"):
			wu = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-g", "--gene_table"):
			gene_table = arg
		elif opt in ("-v", "--dominant"):
			dominant = 1
		elif opt in ("-s", "--plottype"):
			plottype = int(arg)
		elif opt in ("-q", "--subgraph_cut_off="):
			subgraph_cut_off = float(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-a", "--accuracy_cut_off="):
			accuracy_cut_off = float(arg)
	if len(args) == 1:
		stat_table_fname = args[0]
	else:
		stat_table_fname = 'null'
			
	if schema and p_value_cut_off:
		instance = gene_stat(hostname, dbname, schema, table, mcl_table, p_value_cut_off,\
			unknown_cut_off, connectivity_cut_off, recurrence_cut_off, cluster_size_cut_off,\
			leave_one_out, wu, report, \
			log, judger_type, depth_cut_off, dir_files, commit, gene_table, dominant, plottype, \
			stat_table_fname, subgraph_cut_off, debug, accuracy_cut_off)
		instance.dstruc_loadin()
		instance.run()
	else:
		print __doc__
		sys.exit(2)
