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
	-h, --help              show this help

Examples:
	gene_stat.py -k sc_yh60_splat_5 -t cluster_stat2 -m mcl_result2 -p 0.001 -l -w p_table_file
	gene_stat.py -k sc_yh60_splat_5 -t cluster_stat2 -m mcl_result2 -p 0.001 -l -w -c -g p_gene_cluster_stat2 p_table_file

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
	1:	L0 (predicted function is correct if it's one of  the parents of known functions
	2:	common ancestor's depth is larger than the depth_cut_off

"""

import sys, os, psycopg, getopt, csv, fileinput
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

class gene_stat:
	def __init__(self, hostname, dbname, schema, table, mcl_table, p_value_cut_off, unknown_cut_off, \
		connectivity_cut_off, recurrence_cut_off, cluster_size_cut_off, leave_one_out, wu, report=0, \
		log=0, judger_type=0, depth_cut_off =3, dir_files=None, needcommit=0, gene_table='p_gene', dominant=0, plottype=3, stat_table_fname='null'):
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.table = table
		self.mcl_table = mcl_table
		self.p_value_cut_off = float(p_value_cut_off)
		self.unknown_cut_off = float(unknown_cut_off)
		self.connectivity_cut_off = float(connectivity_cut_off)
		self.recurrence_cut_off = int(recurrence_cut_off)
		self.cluster_size_cut_off = int(cluster_size_cut_off)
		self.leave_one_out = int(leave_one_out)
		self.wu = int(wu)
		self.report = int(report)
		
		self.log = int(log)
		self.final_dict = {0: self.final_direct_match,
			1: self.final_L0,
			2: self.final_depth}
		self.final = self.final_dict[judger_type]
		self.depth_cut_off = int(depth_cut_off)
		self.dir_files = dir_files
		self.needcommit = int(needcommit)
		self.gene_table = gene_table
		self.dominant = int(dominant)
		self.plottype = int(plottype)
		self.stat_table_fname = stat_table_fname
		
		self.tp = 0.0
		self.tp_m = 0.0
		self.tp1 = 0.0
		self.tp1_m = 0.0
		self.tn = 0.0
		self.fp = 0.0
		self.fp_m =0.0
		self.fn = 0.0
		#mapping between gene_no and go_no list
		self.known_genes_dict = {}
		self.no_of_records = 0
		if self.log:
			self.log_file = open('/tmp/gene_stat_plot.log','w')
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
		
		#get the non-obsolete biological_process GO DAG
		self.curs.execute("select t2t.term1_id, t2t.term2_id, t1.acc, t2.acc from \
			go.term2term t2t, go.term t1, go.term t2 where t2t.term1_id=t1.id and \
			t2t.term2_id=t2.id and t1.is_obsolete=0 and t2.is_obsolete=0 and \
			t1.term_type='biological_process' and t2.term_type='biological_process' ")
		rows = self.curs.fetchall()
		for row in rows:
		#setup the go_graph structure
			self.go_graph.add_edge(row[2], row[3])
		
		#setup self.go_no2go_id
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
		#before self.final(), self.gene_prediction_dict has been filled.
		self.final()
		
		self.stat_output()
		
		if self.leave_one_out and self.stat_table_fname != 'null':
			#only for leave_one_out approach and when stat_table_fname is not 'null'
			self.stat_table_f = csv.writer(open(self.stat_table_fname, 'w'), delimiter='\t')
			self.go_no_accuracy()
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
		mcl_id = row[0]
		gene_no = row[1]
		p_value_vector = row[2][1:-1].split(',')
		connectivity = float(row[3])
		recurrence_array = row[4][1:-1].split(',')
		recurrence_array = map(int, recurrence_array)
		vertex_set = row[5][1:-1].split(',')
		vertex_set = map(int, vertex_set)
		if len(recurrence_array) < self.recurrence_cut_off:
			#the recurrence is not enough
			return
		if len(vertex_set) > self.cluster_size_cut_off:
			#the cluster is too big
			return
		if connectivity < self.connectivity_cut_off:
			#the cluster is not dense enough
			return
		#transform into float type
		p_value_index_tuple_list = self.index_tuple(p_value_vector)
		for (p_value, index) in p_value_index_tuple_list:
			if self.wu:
			#index 0 corresponds to go_no 0.
				go_no = index
			else:
			#index 0 corresponds to go_no 1
				go_no = index+1
			if self.go_no2depth[go_no] > self.depth_cut_off:
				min_p_value = p_value
				break

		if self.wu:
			if float(p_value_vector[0]) > self.unknown_cut_off:
			#too many unknown genes, and now cut_off is ratio.
				return
		if min_p_value > self.p_value_cut_off:
		#none of the predicted functions in this cluster is significant
			return
		#The cluster is an eligible cluster. Passing all the cut_offs.
		#
		self.no_of_records += 1
		if gene_no not in self.gene_prediction_dict:
			item = gene_prediction()
			self.gene_prediction_dict[gene_no] = item
		self.gene_prediction_dict[gene_no].mcl_id_list.append(mcl_id)	

		for (p_value, index) in p_value_index_tuple_list:
			if p_value > min_p_value:
				break
			elif p_value == min_p_value:
				if self.wu:
				#index 0 corresponds to go_no 0.
					go_no = index
				else:
				#index 0 corresponds to go_no 1
					go_no = index+1
					
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

	def final_L0(self):
		if self.dominant:
			#get the dominant function among all the functions predicted for one gene
			for gene_no in self.gene_prediction_dict:
				entry = self.gene_prediction_dict[gene_no]
				#the biggest support is what the dominant function needs
				dominant_support = max(entry.p_functions_dict.values())
				#delete those functions whose support is less than the dominant support. There's still possibility that >1 functions are kept.
				for go_no in entry.p_functions_dict.keys():
					if entry.p_functions_dict[go_no] < dominant_support:
						del entry.p_functions_dict[go_no]
		for gene_no in self.gene_prediction_dict:
			#initialize three data structures
			L0_dict = {}
			L1_dict = {}
			#L0 or L1 are all good known functions
			good_k_functions_dict = {}
			#each entry is a gene_prediction class
			entry = self.gene_prediction_dict[gene_no]
			p_functions_dict = entry.p_functions_dict.copy()
			p_functions = p_functions_dict.keys()
			if gene_no not in self.known_genes_dict:
				if self.log:
					self.log_file.write('unknown: %d %s %s\n'%(gene_no, repr(p_functions_dict),repr(entry.mcl_id_list)))
				continue
			k_functions_list = self.known_genes_dict[gene_no]
			self.no_of_p_known += 1
			#compare the known go_no with the predicted go_no to see if they match
			#L0 level or L1 level
			for p_go_no in p_functions_dict:
				for k_go_no in k_functions_list:
					if self.is_L0(p_go_no, k_go_no):
						L0_dict[p_go_no] = p_functions_dict[p_go_no]
						good_k_functions_dict[k_go_no] = 1
					elif self.is_L1(p_go_no, k_go_no):
						L1_dict[p_go_no] = p_functions_dict[p_go_no]
						good_k_functions_dict[k_go_no] = 1
			#if one L1 is also counted as L0, remove it from L1_dict.
			#this case happens when one known function has both L0 and L1 matches in the predicted functions
			for go_no in L0_dict:
				if go_no in L1_dict:
					del L1_dict[go_no]
			entry.tp = L0_dict
			entry.tp1 = L1_dict
			for k_go_no in k_functions_list:
				if k_go_no not in good_k_functions_dict:
					entry.fn[k_go_no] = 1
			for p_go_no in p_functions_dict:
				if p_go_no not in L0_dict and p_go_no not in L1_dict:
					entry.fp[p_go_no] = p_functions_dict[p_go_no]
			
			entry.tn = self.no_of_functions - (len(entry.tp)+len(entry.tp1)+len(entry.fp)+len(entry.fn))
			self.tp += len(entry.tp)
			self.tp_m += sum(entry.tp.values())
			self.tp1 += len(entry.tp1)
			self.tp1_m += sum(entry.tp1.values())
			self.tn += entry.tn
			self.fp += len(entry.fp)
			self.fp_m += sum(entry.fp.values())
			self.fn += sum(entry.fn.values())

	def final_depth(self):
		if self.dominant:
			#get the dominant function among all the functions predicted for one gene
			for gene_no in self.gene_prediction_dict:
				entry = self.gene_prediction_dict[gene_no]
				#the biggest support is what the dominant function needs
				dominant_support = max(entry.p_functions_dict.values())
				#delete those functions whose support is less than the dominant support. There's still possibility that >1 functions are kept.
				for go_no in entry.p_functions_dict.keys():
					if entry.p_functions_dict[go_no] < dominant_support:
						del entry.p_functions_dict[go_no]
		for gene_no in self.gene_prediction_dict:
			#initialize three data structures
			L0_dict = {}
			L1_dict = {}
			#each entry is a gene_prediction class
			entry = self.gene_prediction_dict[gene_no]
			p_functions_dict = entry.p_functions_dict.copy()
			if gene_no not in self.known_genes_dict:
				#unknown genes are not validatable. Ignore.
				#self.log_file.write('unknown: %d %s %s\n'%(gene_no, repr(p_functions_dict),repr(entry.mcl_id_list)))
				continue
			k_functions_set = self.known_genes_dict[gene_no]
			self.no_of_p_known += 1
			#compare the known go_no with the predicted go_no to see if they match

			for p_go_no in entry.p_functions_dict:
				if p_go_no in k_functions_set:
					L0_dict[p_go_no] = p_functions_dict[p_go_no]
					self.gene_prediction_dict[gene_no].p_functions_struc_dict[p_go_no].is_correct = 1
					del p_functions_dict[p_go_no]
					k_functions_set.remove(p_go_no)
				elif self.common_ancestor_deep_enough(p_go_no, k_functions_set):
					#Jasmine's advice.
					if self.log:
						self.log_file.write('%s\t%s\t'%(self.gene_no2gene_id[gene_no], repr(k_functions_set)))
					L0_dict[p_go_no] = 1
					self.gene_prediction_dict[gene_no].p_functions_struc_dict[p_go_no].is_correct = 1
					del p_functions_dict[p_go_no]
			
			#remaining function categories in k_functions_set are false negatives
			entry.fn = k_functions_set

			entry.tp = L0_dict
			entry.tp1 = L1_dict
			#the p_go_no's left in p_functions_dict are false positive's
			entry.fp = p_functions_dict
			
			entry.tn = self.no_of_functions - (len(entry.tp)+len(entry.tp1)+len(entry.fp)+len(entry.fn))
			self.tp += len(entry.tp)
			self.tp_m += sum(entry.tp.values())
			self.tp1 += len(entry.tp1)
			self.tp1_m += sum(entry.tp1.values())
			self.tn += entry.tn
			self.fp += len(entry.fp)
			self.fp_m += sum(entry.fp.values())
			self.fn += len(entry.fn)
			#self.log_file.write('known: %d %s %s %d %s %s %s %s\n'%(gene_no, repr(entry.tp),repr(entry.tp1),\
			#	entry.tn,repr(entry.fp),repr(entry.fn),repr(p_functions),repr(entry.mcl_id_list)))

	def common_ancestor_deep_enough(self, p_go_no, k_functions_set):
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
			elif self.log:
				self.log_file("distance for %s doesn't exist.\n"%(repr(key)))
		
		for ancestor in ancestor_set:
			depth = self.go_term_id2depth[ancestor]
			if depth >= self.depth_cut_off:
				if self.log:
					self.log_file.write("%s's common_ancestor %s\n"%(self.go_no2go_id[p_go_no], \
						self.go_no2go_id[self.go_term_id2go_no[ancestor]]))
				#pre-stop the loop
				break
		if depth >= self.depth_cut_off:
			return 1
		else:
			return 0
		

	def final_direct_match(self):
		if self.dominant:
			#get the dominant function among all the functions predicted for one gene
			for gene_no in self.gene_prediction_dict:
				entry = self.gene_prediction_dict[gene_no]
				#the biggest support is what the dominant function needs
				dominant_support = max(entry.p_functions_dict.values())
				#delete those functions whose support is less than the dominant support. There's still possibility that >1 functions are kept.
				for go_no in entry.p_functions_dict.keys():
					if entry.p_functions_dict[go_no] < dominant_support:
						del entry.p_functions_dict[go_no]
		for gene_no in self.gene_prediction_dict:
			#matched predictions
			L0_dict = {}
			#each entry is a gene_prediction class
			entry = self.gene_prediction_dict[gene_no]
			p_functions_dict = entry.p_functions_dict.copy()
			if gene_no not in self.known_genes_dict:
				#unknown genes are not validatable. Ignore.
				#self.log_file.write('unknown: %d %s %s\n'%(gene_no, repr(p_functions_dict),repr(entry.mcl_id_list)))
				continue
			k_functions_set = self.known_genes_dict[gene_no]
			self.no_of_p_known += 1
			#compare the known go_no with the predicted go_no to see if they match

			for k_go_no in k_functions_set:
				if k_go_no in p_functions_dict:
					L0_dict[k_go_no] = p_functions_dict[k_go_no]
					self.gene_prediction_dict[gene_no].p_functions_struc_dict[k_go_no].is_correct = 1
					del p_functions_dict[k_go_no]
				else:
					entry.fn[k_go_no] = 1


			entry.tp = L0_dict
			#the p_go_no's left in p_functions_dict are false positive's
			entry.fp = p_functions_dict
			
			entry.tn = self.no_of_functions - (len(entry.tp)+len(entry.fp)+len(entry.fn))
			self.tp += len(entry.tp)
			self.tp_m += sum(entry.tp.values())
			self.tn += entry.tn
			self.fp += len(entry.fp)
			self.fp_m += sum(entry.fp.values())
			self.fn += sum(entry.fn.values())
			#self.log_file.write('known: %d %s %s %d %s %s %s %s\n'%(gene_no, repr(entry.tp),repr(entry.tp1),\
			#	entry.tn,repr(entry.fp),repr(entry.fn),repr(p_functions),repr(entry.mcl_id_list)))

	def is_L0(self, p_go_no, k_go_no):
		k_go_id = self.go_no2go_id[k_go_no]
		p_go_id = self.go_no2go_id[p_go_no]
		k_go_family = Set(self.go_graph.forw_bfs(k_go_id))
		if p_go_id in k_go_family:
			#self.log_file.write('%d is L0 of %d:: %s %s\n'%(p_go_no, k_go_no, p_go_id, k_go_id))
			return 1
		else:
			return 0
		
	def is_L1(self, p_go_no, k_go_no):
		k_go_id = self.go_no2go_id[k_go_no]
		p_go_id = self.go_no2go_id[p_go_no]
		k_go_inc_nbrs = Set(self.go_graph.inc_nbrs(k_go_id))
		if p_go_id in k_go_inc_nbrs:
			#self.log_file.write("%d is direct parent of %d:: %s %s\n"%(p_go_no, k_go_no, p_go_id, k_go_id))	
			return 1
		for k_go_inc_nbr in k_go_inc_nbrs:
			k_go_inc_nbr_out_nbrs = Set(self.go_graph.out_nbrs(k_go_inc_nbr))
			if p_go_id in k_go_inc_nbr_out_nbrs:
				#self.log_file.write("%d and %d are siblings:: %s %s\n"%(p_go_no, k_go_no, p_go_id, k_go_id))			
				return 1
		return 0
	
	def list_stringlist(self, list):
		return '{' + repr(list)[1:-1] + '}'
	
	def submit(self):
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
				recurrence_cut_off      integer,\
				connectivity_cut_off    float,\
				cluster_size_cut_off    integer,\
				unknown_cut_off      float,\
				depth_cut_off integer\
				)"%self.gene_table)

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
					values(%d, %d, %d, %f, %f, %d, '%s', ARRAY%s, %f, %d, %f, %d, %f, %d)"%\
					(self.gene_table, gene_no, go_no, unit[go_no].is_correct, sum(unit[go_no].p_value_list)/no_of_clusters,\
					e_accuracy, no_of_clusters, ';'.join(context_list),\
					repr(unit[go_no].cluster_array), self.p_value_cut_off, self.recurrence_cut_off,\
					self.connectivity_cut_off, self.cluster_size_cut_off, self.unknown_cut_off, self.depth_cut_off))
		
		if self.needcommit:				
			self.curs.execute("end")	
		sys.stderr.write("done.\n")

	def stat_output(self):
		sys.stderr.write('\n\tp_value_cut_off:%f unknown_cut_off:%f connectivity_cut_off:%f\n'%(self.p_value_cut_off, self.unknown_cut_off, self.connectivity_cut_off))
		sys.stderr.write('\trecurrence_cut_off:%d cluster_size_cut_off:%d\n'%(self.recurrence_cut_off, self.cluster_size_cut_off))
		sys.stderr.write('\tdepth_cut_off:%d\n'%(self.depth_cut_off))
		sys.stderr.write('\tTotal genes: %d\n'%len(self.gene_prediction_dict))
		sys.stderr.write('\tTotal known genes: %d\n'%self.no_of_p_known)
		sys.stderr.write("\tBased on functions:\n")
		sys.stderr.write('\t\tTP0: %d  TP1: %d  TN: %d  FP: %d  FN: %d\n'%(self.tp, self.tp1, self.tn, self.fp, self.fn))
		if (self.tp+self.tp1+self.fn) == 0:
			sys.stderr.write('\t\tSensitvity: Null\n')
		else:
			sys.stderr.write('\t\tSensitvity: %f\n'%((self.tp+self.tp1)/(self.tp+self.tp1+self.fn)))
		if (self.fp+self.tn) == 0:
			sys.stderr.write('\t\tSpecificity: Null\n')
		else:
			sys.stderr.write('\t\tSpecificity: %f\n'%(self.tn/(self.fp+self.tn)))
		if (self.tp+self.tp1+self.fp) == 0:
			sys.stderr.write('\t\tFalse Positive Ratio: Null\n')
		else:
			sys.stderr.write('\t\tFalse Positive Ratio: %f\n'%(self.fp/(self.tp+self.tp1+self.fp)))
		sys.stderr.write("\tBased on clusters:\n")
		sys.stderr.write('\t\tTP0_M: %d  TP1_M: %d  FP_M: %d\n'%(self.tp_m, self.tp1_m, self.fp_m))
		if (self.tp_m+self.tp1_m+self.fp_m) == 0:
			sys.stderr.write('\t\tFalse Positive Ratio: Null\n')
		else:
			sys.stderr.write('\t\tFalse Positive Ratio: %f\n'%(self.fp_m/(self.tp_m+self.tp1_m+self.fp_m)))

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
		header = ['gene_id', 'function_known', 'function_predicted', 'is_correct', 'average p_value', 'expected accuracy', '#supporting clusters', 'cluster_context']
		self.stat_table_f.writerow(header)
		for gene_no in self.gene_prediction_dict:
			#first output only the known genes
			if gene_no in self.known_genes_dict:
				unit = self.gene_prediction_dict[gene_no].p_functions_struc_dict
				self._table_output(self.stat_table_f, gene_no, unit)
		
		for gene_no in self.gene_prediction_dict:
			#second, output the unknown genes
			if gene_no not in self.known_genes_dict:
				unit = self.gene_prediction_dict[gene_no].p_functions_struc_dict
				self._table_output(self.stat_table_f, gene_no, unit)			

	def _table_output(self, f_handler, gene_no, unit):
		#indicator for the number of functions associated with this gene
		for go_no in unit:
			row = []
			row.append(self.gene_no2gene_id[gene_no])
			#function_known might not be available
			if gene_no in self.known_genes_dict:
				function_known_list = []
				for go_known in self.known_genes_dict[gene_no]:
					function_known_list.append(self.go_no2go_name[go_known])
				row.append(';'.join(function_known_list))
			else:
				row.append('')

			row.append(self.go_no2go_name[go_no])
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

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema=", "table=", "mcl_table=", "p_value_cut_off=",\
		"unknown_cut_off=", "connectivity_cut_off=", "recurrence_cut_off=", "cluster_size_cut_off=", "depth_cut_off=",\
		"log", "judger_type=", "dir_files=", "leave_one_out", "wu", "report", "commit", "gene_table=", "dominant", "plottype"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:m:p:u:n:y:x:e:oj:f:lwrcg:vs:", long_options_list)
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
			recurrence_cut_off = int(arg)
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
	if len(args) == 1:
		stat_table_fname = args[0]
	else:
		stat_table_fname = 'null'
			
	if schema and p_value_cut_off:
		instance = gene_stat(hostname, dbname, schema, table, mcl_table, p_value_cut_off,\
			unknown_cut_off, connectivity_cut_off, recurrence_cut_off, cluster_size_cut_off,\
			leave_one_out, wu, report, \
			log, judger_type, depth_cut_off, dir_files, commit, gene_table, dominant, plottype, stat_table_fname)
		instance.dstruc_loadin()
		instance.run()
	else:
		print __doc__
		sys.exit(2)
