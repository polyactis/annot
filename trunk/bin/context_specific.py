#!/usr/bin/env python
"""
Usage: context_specific.py -k SCHEMA -t gene_p_table -m mcl_table -g p_gene_table [OPTION] OUTPUT_FILE

Option:
	OUTPUT_FILE is the file to store the function table of all predicted genes.
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	gene_p_table
	-m ..., --mcl_table=...	mcl_table corresponding to above table
	-g ..., --gene_table=...	table storing the stat results
	-n ..., --node_dist_table=...	the node_distance table, node_dist(default) (IGNORE)
	-s ..., --contrast=...	the contrast differential ratio for two contexts to be different, 0.3(default)
	-y ..., --type=...	1(all genes, default), 2(known genes), 3(unknown genes).
	-r, --report	report the progress(a number) IGNORE
	-c, --commit	commit the database transaction, records in table gene. IGNORE
	-l, --log	record down some stuff in the logfile(context_specific.log)
	-h, --help              show this help

Examples:
	context_specific.py -k sc_54 -g p_gene -s 0.4 -t gene_p -m mcl_result
		context_specific.out

Description:
	This program inspects the cluster contexts of the functions predicted for one gene.
	It only keeps those functions which only have distinct contexts. Functions with similar
	contexts will be merged.
	
	03-09-05
		Below in {} is dropped.
	{
	The depth of the predicted function on the GO tree is drawn as
	a histogram in 'depth_hist.pdf'.
	The max_raw_distance, max_lee_distance, max_jasmine_distance
	among the predicted functions for one gene
	are drawn as a histogram in files with those names respectively.
	}
	
"""

import sys, os, getopt, csv
from sets import Set
from rpy import r
from numarray import array, greater
from codense.common import db_connect

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
	"""
	data structure for p_functions_struc_dict in gene_prediction
	03-09-05
		add three is_correct's
	"""
	def __init__(self):
		self.is_correct = 0
		self.is_correct_l1 = 0
		self.is_correct_lca = 0
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
	"""

	run()
		--dstruc_loadin()
		(loop)
			--distinct_contexts()
				--is_context_diff()
			--table_output()
		--stat_output()
	
	03-09-05
		overhaul the whole program to meet the new p_gene_table and gene_p_table.
		related to another pipeline change on 02-21-05.
	"""
	def __init__(self, hostname='zhoub', dbname='graphdb', schema=None, table=None, \
		mcl_table=None, node_dist_table=None, contrast=0.3, type=1, \
		report=0, needcommit=0, log=0, gene_table='p_gene', stat_table_fname=None):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema	
		
		self.table = table
		self.mcl_table = mcl_table
		self.node_dist_table = node_dist_table
		self.contrast = float(contrast)
		self.type = int(type)
		self.report = int(report)
		self.needcommit = int(needcommit)
		self.gene_table = gene_table
		self.stat_table_fname = stat_table_fname
		self.log = int(log)
		if self.log:
			self.log_file = open('/tmp/context_specific.log','w')	
		#mapping between gene_no and go_no list
		self.known_genes_dict = {}
		#some counters
		self.no_of_records = 0
		#a dictionary recording how many predictions before and after the context merging
		#key is gene_no, value is [#before, #after]
		self.gene_no2no_of_predictions = {}

		self.gene_prediction_dict = {}
		self.no_of_p_known = 0
		self.go_no2go_id = {}	
		#mapping between go_no and go's name
		self.go_no2go_name = {}
		#mapping between gene_no and gene_id
		self.gene_no2gene_id = {}
		#the overall prediction accuracy of each function
		self.go_no2accuracy = {}
		#mapping each go_no to its depth
		self.go_no2depth = {}
		#a list of the depths of all the predicted functions, for histogram
		self.depth_list=[]
		#a list storing the three kinds of maximum distances between the predicted functions for one gene.
		self.max_raw_distance_list = []
		self.max_lee_distance_list = []
		self.max_jasmine_distance_list = []
		#mapping between a pair of go_no's and its associated distances
		self.go_no2distance = {}
		#mapping go term id and go_no
		self.go_term_id2go_no = {}
		
	def dstruc_loadin(self, curs):
		"""
		
		03-09-05
			get the context from mcl_table via linking through mcl_id of p_gene_table
			context_dict is set
		"""
		from codense.common import get_known_genes_dict, get_go_no2go_id, \
			get_go_no2name, get_gene_no2gene_id
		
		self.known_genes_dict = get_known_genes_dict(curs)
		self.go_no2go_id = get_go_no2go_id(curs)
		self.go_no2go_name = get_go_no2name(curs)
		self.gene_no2gene_id = get_gene_no2gene_id(curs)
		
		sys.stderr.write("Setting up gene_prediction_dict...")
		#setup self.gene_prediction_dict
		curs.execute("select p.gene_no, p.go_no, p.is_correct, p.is_correct_l1, p.is_correct_lca, m.vertex_set\
			from %s p, %s g, %s m where g.p_gene_id=p.p_gene_id and m.mcl_id=p.mcl_id"%(self.gene_table, self.table, self.mcl_table))
		rows = curs.fetchall()
		for row in rows:
			gene_no = row[0]
			if self.type==2 and gene_no not in self.known_genes_dict:
				#I only want the known genes, but this gene is unknown
				continue
			elif self.type==3 and gene_no in self.known_genes_dict:
				#i only want the unknown genes, but this gene is known
				continue
			go_no = row[1]
			is_correct = row[2]
			is_correct_l1 = row[3]
			is_correct_lca = row[4]
			vertex_set = row[5][1:-1].split(',')
			vertex_set = map(int, vertex_set)

			item = function_struc()
			item.is_correct = is_correct
			item.is_correct_l1 = is_correct_l1
			item.is_correct_lca = is_correct_lca
			#context_dict is a set
			item.context_dict = Set(vertex_set)
			if gene_no not in self.gene_prediction_dict:
				self.gene_prediction_dict[gene_no] = gene_prediction()
				self.gene_prediction_dict[gene_no].p_functions_struc_dict[go_no] = item
			else:
				self.gene_prediction_dict[gene_no].p_functions_struc_dict[go_no] = item
		

		sys.stderr.write("Done\n")
		
		"""
		#setup self.go_no2distance
		self.curs.execute("select go_id1, go_id2, raw_distance, lee_distance, jasmine_distance\
			from go.%s"%self.node_dist_table)
		rows = self.curs.fetchall()
		for row in rows:
			go_no1 = self.go_term_id2go_no.get(row[0])
			go_no2 = self.go_term_id2go_no.get(row[1])
			if go_no1 and go_no2:
				#key tuple in ascending order
				if go_no1<go_no2:
					self.go_no2distance[(go_no1, go_no2)] = (row[2], row[3], row[4])
				else:
					self.go_no2distance[(go_no2, go_no1)] = (row[2], row[3], row[4])
		"""
		

	def cluster_context2dict(self, cluster_context):
		"""
		convert a cluster_context string fetched from database to a context_dict
		
		03-09-05
			defunct, no use anymore
		"""
		context_dict = {}
		complex_gene_list = cluster_context.split(';')
		for complex_gene in complex_gene_list:
			gene_no, support = map(int, complex_gene.split('/'))
			context_dict[gene_no] = support
		return context_dict
	

	def run(self):
		"""
		03-09-05
			no conversion for context_dict from dict to set
			no distance stuff
		"""
		(conn, curs) = db_connect(self.hostname, self.dbname, self.schema)
		self.dstruc_loadin(curs)
		if self.stat_table_fname:
			self.stat_table_f = csv.writer(open(self.stat_table_fname, 'w'), delimiter='\t')

		for (gene_no,unit) in self.gene_prediction_dict.iteritems():
			self.gene_no2no_of_predictions[gene_no] = [len(unit.p_functions_struc_dict)]

			"""
			#rearrange the context_dict into a go_no:Set structure
			self.go_merge_dict = {}
			go_no_list = []
			for go_no in unit.p_functions_struc_dict:
				#!!! MEANING CHANGE !!!#
				#function_struc.context_dict will be a Set data structure to store context genes
				#function_struc.cluster_array will store the go terms to be merged
				go_no_list.append(go_no)
				depth = self.go_no2depth[go_no]
				self.depth_list.append(depth)
				
				item = function_struc()
				item.context_dict = Set(unit.p_functions_struc_dict[go_no].context_dict.keys())
				item.is_correct = unit.p_functions_struc_dict[go_no].is_correct
				self.go_merge_dict[go_no] = item
			
			"""
			
			"""
			
			max_raw_distance = 0
			max_lee_distance = 0
			max_jasmine_distance = 0
		
			for i in range(len(go_no_list)):
				for j in range(i+1, len(go_no_list)):
					go_no1 = go_no_list[i]
					go_no2 = go_no_list[j]
					if go_no1 < go_no2:
						key = (go_no1, go_no2)
					else:
						key = (go_no2, go_no1)
					raw_distance, lee_distance, jasmine_distance = self.go_no2distance[key]
					if raw_distance > max_raw_distance:
						max_raw_distance = raw_distance
					if lee_distance > max_lee_distance:
						max_lee_distance = lee_distance
					if jasmine_distance > max_jasmine_distance:
						max_jasmine_distance = jasmine_distance
			
			self.max_raw_distance_list.append(max_raw_distance)
			self.max_lee_distance_list.append(max_lee_distance)
			self.max_jasmine_distance_list.append(max_jasmine_distance)
			if self.log:
				self.log_file.write('%d\t%s\t%d\t%d\t%d\n'%(gene_no, repr(go_no_list), max_raw_distance,\
					max_lee_distance, max_jasmine_distance))
			"""
			#compute the distinct contexts
			function_struc_dict = self.distinct_contexts(unit.p_functions_struc_dict)
			self.gene_no2no_of_predictions[gene_no].append(len(function_struc_dict))
			
			if self.stat_table_fname:
				self.table_output(self.stat_table_f, gene_no, function_struc_dict)

		self.stat_output(self.stat_table_f, self.gene_no2no_of_predictions)
		"""
		self.plot('depth_hist.pdf', self.depth_list, 'depth histogram', 'depth of predicted function')
		self.plot('max_raw_distance_hist.pdf', self.max_raw_distance_list, 'max_raw_distance histogram', \
			'max_raw_distance of predicted functions of one gene')
		self.plot('max_lee_distance_hist.pdf', self.max_lee_distance_list, 'max_lee_distance histogram', \
			'max_lee_distance of predicted functions of one gene')
		self.plot('max_jasmine_distance_hist.pdf', self.max_jasmine_distance_list, 'max_jasmine_distance histogram', \
			'max_jasmine_distance of predicted functions of one gene')
		"""
		
	def distinct_contexts(self, go_merge_dict):
		"""
		shrink the go_nos based on their contexts,
		among those go_nos which have similar contexts(is_context_diff()), only keep one go_no
		and merge their contexts
		03-09-05
			new strategy, a similar one in gene_p_map_redundancy.py and p_gene_analysis.py
		"""
		#a mapping between go_no's based on the similarity of their contexts
		go_no_map = {}
		go_no_list = go_merge_dict.keys()
		no_of_gos = len(go_no_list)
		for i in range(no_of_gos):
			go_no1 = go_no_list[i]
			if go_no1 not in go_no_map:
				#initial encountering, not mapped, map to itself
				go_no_map[go_no1] = go_no1
				for j in range(i+1, no_of_gos):
					go_no2 = go_no_list[j]
					if go_no2 not in go_no_map:
						#otherwise, it's been mapped already.
						is_similar = self.is_context_diff(go_merge_dict[go_no1].context_dict, go_merge_dict[go_no2].context_dict)
						if is_similar:
							#map it to go_no1
							go_no_map[go_no2] = go_no1
							#the merged go term is correct if and only if i and j are both correct
							go_merge_dict[go_no1].is_correct *= go_merge_dict[go_no2].is_correct
							go_merge_dict[go_no1].is_correct_l1 *= go_merge_dict[go_no2].is_correct_l1
							go_merge_dict[go_no1].is_correct_lca *= go_merge_dict[go_no2].is_correct_lca
							#put j into the list to be merged by i
							go_merge_dict[go_no1].cluster_array.append(go_no2)
							#union_update the context Set of i
							go_merge_dict[go_no1].context_dict |= go_merge_dict[go_no2].context_dict

		#the data structure to return
		function_struc_dict = {}
		for go_no in go_no_map.values():
			#the one being mapped
			function_struc_dict[go_no] = go_merge_dict[go_no]
		return function_struc_dict

	def is_context_diff(self, set1, set2):
		"""
		to see whether two context sets are different or not.
		"""
		#defaultly, regard two sets as different
		is_similar = 0
		intersection_set = set1 & set2
		set1_extra = set1 - intersection_set
		set2_extra = set2 - intersection_set
		#two important ratios
		r1 = len(set1_extra)/float(len(set1))
		r2 = len(set2_extra)/float(len(set2))
		#if either of them are below the threshold, they are not different. Need discussion.
		if r1 < self.contrast or r2 < self.contrast:
			is_similar = 1
		return is_similar

	
	def table_output(self, stat_table_f, gene_no, go_merge_dict, max_raw_distance=0, max_lee_distance=0, max_jasmine_distance=0):
		"""
		03-09-05
			set the default values of those distances
			discard the distance stuff
		"""
		gene_id = self.gene_no2gene_id[gene_no]
		row = ['gene_id', 'function group', 'is_correct', 'is_correct_l1', \
			'is_correct_lca', 'context']
		stat_table_f.writerow(row)
		
		for go_no in go_merge_dict:
			#row is for output
			row = []
			row.append(gene_id)
			unit = go_merge_dict[go_no]
			#the go_no's merged + itself
			go_no_list = unit.cluster_array
			go_no_list.append(go_no)
			
			go_name_list = []
			for term_no in go_no_list:
				go_name_list.append('%s(%s)'%(self.go_no2go_name[term_no], self.go_no2go_id[term_no]) )
			
			#row.append('%s'%('/'.join(map(repr,go_no_list) )))
			row.append('%s'%('/'.join(go_name_list)))
			row.append(unit.is_correct)
			row.append(unit.is_correct_l1)
			row.append(unit.is_correct_lca)
			context_gene_id_list = []
			for gene_no in unit.context_dict:
				context_gene_id_list.append(self.gene_no2gene_id[gene_no])
			#context_gene_no_list = list(unit.context_dict)
			#context_gene_no_list.sort()
			#row.append('%s'%('/'.join( map(repr, context_gene_no_list) ) ))
			row.append('%s'%('/'.join(context_gene_id_list)))
			"""
			#append the three kinds of maximum distances
			row.append(max_raw_distance)
			row.append(max_lee_distance)
			row.append(max_jasmine_distance)
			"""

			stat_table_f.writerow(row)
	
	def stat_output(self, stat_table_f, gene_no2no_of_predictions):
		"""
		03-09-05
			
		"""
		sys.stderr.write("Outputting stats ... ")
		no_of_predictons_array = array(gene_no2no_of_predictions.values())
		before_array = no_of_predictons_array[:,0]
		after_array = no_of_predictons_array[:,1]
		no_of_predicted_genes = len(gene_no2no_of_predictions)
		genes_with_multiple_functions = sum(greater(before_array, 1))
		genes_with_multiple_contexts = sum(greater(after_array, 1))
		
		avg_functions_per_gene = sum(before_array)/float(no_of_predicted_genes)
		avg_contexts_per_gene = sum(after_array)/float(no_of_predicted_genes)
		
		stat_table_f.writerow(['Total genes predicted: %s'%no_of_predicted_genes])
		stat_table_f.writerow(['\t%s(%f) of predicted genes with multiple functions'%(genes_with_multiple_functions, \
			genes_with_multiple_functions/float(no_of_predicted_genes)) ])
		stat_table_f.writerow(['\t%s(%f) of predicted genes with multiple distinct contexts'%(genes_with_multiple_contexts,\
			genes_with_multiple_contexts/float(no_of_predicted_genes)) ])
		stat_table_f.writerow(['\taverage functions per gene: %f'%(avg_functions_per_gene)])
		stat_table_f.writerow(['\taverage distinct contexts per gene: %f'%(avg_contexts_per_gene)])
		sys.stderr.write("Done.\n")
	
	def plot(self, filename, list_to_plot, main_lab, xlab):
		max_length = max(list_to_plot)
		r.pdf(filename)
		r.hist(list_to_plot, breaks=range(max_length+1), las=1, main=main_lab, xlab=xlab)
		r.dev_off()
		
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema=", "table=", "mcl_table=", \
		"contrast=", "type=", "node_dist_table=", "report", "commit", "log", "gene_table="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:m:n:s:y:rclg:", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = None
	table = None
	mcl_table = None
	node_dist_table = 'node_dist'
	contrast = 0.3
	type = 1
	report = 0
	commit = 0
	log = 0
	gene_table = None
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
		elif opt in ("-n", "--node_dist_table"):
			node_dist_table = arg
		elif opt in ("-s", "--contrast"):
			contrast = float(arg)
		elif opt in ("-y", "--type"):
			type = int(arg)
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-l", "--log"):
			log = 1
		elif opt in ("-g", "--gene_table"):
			gene_table = arg
	if len(args) == 1:
		stat_table_fname = args[0]
	else:
		stat_table_fname = None
			
	if schema and gene_table and table and mcl_table and stat_table_fname:
		instance = context_specific(hostname, dbname, schema, table, mcl_table, \
			node_dist_table, contrast, type, report, \
			commit, log, gene_table, stat_table_fname)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
