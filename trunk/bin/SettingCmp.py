#!/usr/bin/env python
"""
Usage: SettingCmp.py -k -i -j -o 

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-i ...,	fname of setting1
	-l ...,	lm_bit of setting1, ('111', default)
	-a ...	acc_cutoff of setting1, (0.6 default)
	-j ...,	fname of setting2
	-m ...,	lm_bit of setting2, ('111' default)
	-e ...,	acc_cutoff of setting2, (0.6 default)
	-o ...,	output fname
	-n ...,	no_of_samples for each category, 30, (default)
	-b ...,	setting linear model bit, 11(default), 1=linear model present
	-s,	simple output, just the output fname table, no graph pictures.
	-c,	commit the database transaction
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	~/script/annot/bin/SettingCmp.py -k hs_fim_40 -i hs_fim_40m4x40rec0_8 -l 11100
		-j hs_fim_40m4x40e1s1__1_0l10 -m 11001 -o /tmp/cat.cmp

Description:
	A program to show the difference between two settings by sampling
	good predictions from 4 categories.
	Category 1: only in setting1. prediction from 1, score2 and is_accepted2 refer to 1
	Category 2: in both setting1 and setting2, prediction info from 2.
	Category 3: only in setting2, prediction info from 2.
	Category 4: neither(bad predictions of setting2, not setting1), prediction from 2.
	
	prediction information, score2 and is_accepted2 are all based on setting2.
	score1 and is_accepted1 are got by applying setting1's linear model coeffs
	to setting1's prediction infomation. But only setting2's prediction information is
	shown.(10-18-05)

"""

import sys, os, getopt, csv, math
sys.path += [os.path.expanduser('~/script/annot/bin')]
from codense.common import db_connect, p_gene_id_set_from_gene_p_table, form_schema_tables
from sets import Set
from codense.common import get_gene_no2gene_id, get_gene_no2go_no
from cluster_info import cluster_info	#transform vertex_set and edge_set into subgraph
from codense.common import system_call, graphDotOutput
from p_gene_analysis import p_gene_analysis	#related to linear model coeff
from MpiPredictionFilter import prediction_attributes	#to accept data from sql select

class SettingCmp:
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, fname1=None, \
		lm_bit1=None, acc_cutoff1=None, fname2=None, lm_bit2=None, acc_cutoff2=None, \
		ofname=None, no_of_samples=30, bit='11', simple=0, commit=0, report=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.fname1 = fname1
		self.lm_bit1 = lm_bit1
		self.acc_cutoff1 = float(acc_cutoff1)
		self.fname2 = fname2
		self.lm_bit2 = lm_bit2
		self.acc_cutoff2 = float(acc_cutoff2)
		self.ofname = ofname
		self.no_of_samples = int(no_of_samples)
		self.bit = bit
		self.simple = int(simple)
		self.commit = int(commit)
		self.report = int(report)
		
		self.category_no2information = {
			1:"Only in %s,%s,%s"%(self.fname1, self.lm_bit1, self.acc_cutoff1),
			2:"In both %s,%s,%s and %s,%s,%s"%(self.fname1, self.lm_bit1, self.acc_cutoff1, \
				self.fname2, self.lm_bit2, self.acc_cutoff2),
			3:"Only in %s,%s,%s"%(self.fname2, self.lm_bit2, self.acc_cutoff2),
			4:"Neither"}
	
	def sample_p_gene_id_set(self, p_gene_id_set, no_of_samples):
		"""
		10-15-05
			return a sample from p_gene_id_set
		10-18-05
			handle the situation that step=0(no_of_samples>len(p_gene_id_set)
		"""
		sys.stderr.write("Sampling...")
		step = len(p_gene_id_set)/no_of_samples
		if step==0:
			index_ls = range(0, len(p_gene_id_set))
		else:
			index_ls = range(0, len(p_gene_id_set), step)
		p_gene_id_list = list(p_gene_id_set)
		sys.stderr.write("End sampling.\n")
		return [p_gene_id_list[index] for index in index_ls]
	
	def output_lm_model(self, curs, schema_instance, writer):
		"""
		10-16-05
			output two linear model coefficients
		"""
		sys.stderr.write("Outputting linear model coefficients...")		
		curs.execute("select * from %s"%schema_instance.lm_table)
		rows = curs.fetchall()
		for row in rows:
			row = list(row)
			row[0] = schema_instance.lm_suffix
			writer.writerow(row)
		writer.writerow([])
		sys.stderr.write("End linear model output.\n")
	
	def output_p_gene_id_list(self, curs, schema_instance1, schema_instance2, p_gene_id_list, writer, pic_output_dir,\
		pga_instance1, pga_instance2, cluster_info_instance, simple):
		"""
		10-15-05
			add score1 and is_accepted1
		10-17-05
			score and is_accepted depend on whether pga_instance is None or not
		10-17-05 add simple to allow no graph pictures output
			also get prediction from schema_instance1 and calculate the score if prediction is available
		10-18-05
			sort the p_gene_id_list first
		"""
		#10-15-05 following sentence slightly different from PredictionFilterByClusterSize.py in the trailing edge_gradient
			#and d_matrix is a placeholder
		sql_sentence1 = "SELECT p.p_gene_id, p.gene_no, p.go_no, p.is_correct, p.is_correct_l1, \
			p.is_correct_lca, p.avg_p_value, p.no_of_clusters, p.cluster_array, p.p_value_cut_off, p.recurrence_cut_off, \
			p.connectivity_cut_off, p.cluster_size_cut_off, p.unknown_cut_off, p.depth_cut_off, p.mcl_id, p.lca_list, \
			p.vertex_gradient, p.edge_gradient, m.vertex_set, s.edge_set, 'd_matrix', 'r' from %s p, %s s, %s m where \
			p.mcl_id=s.splat_id and p.mcl_id=m.mcl_id"%(schema_instance1.p_gene_table, \
			schema_instance1.splat_table, schema_instance1.mcl_table)
		sql_sentence2 = "SELECT p.p_gene_id, p.gene_no, p.go_no, p.is_correct, p.is_correct_l1, \
			p.is_correct_lca, p.avg_p_value, p.no_of_clusters, p.cluster_array, p.p_value_cut_off, p.recurrence_cut_off, \
			p.connectivity_cut_off, p.cluster_size_cut_off, p.unknown_cut_off, p.depth_cut_off, p.mcl_id, p.lca_list, \
			p.vertex_gradient, p.edge_gradient, m.vertex_set, s.edge_set, 'd_matrix', 'r' from %s p, %s s, %s m where \
			p.mcl_id=s.splat_id and p.mcl_id=m.mcl_id"%(schema_instance2.p_gene_table, \
			schema_instance2.splat_table, schema_instance2.mcl_table)
		writer.writerow(['p_gene_id', 'gene_no', 'go_no', 'is_correct_lca', 'p_value', 'recurrence', 'connectivity',\
			'cluster_size', 'unknown_ratio', 'mcl_id', 'lca_list', 'edge_gradient', 'score1', 'is_accepted1', 'score2', 'is_accepted2'])
		p_gene_id_list.sort()
		for p_gene_id in p_gene_id_list:
			#sql_sentence1's prediction infomation is not gonna be displayed
			curs.execute("%s and p.p_gene_id=%s"%(sql_sentence1, p_gene_id))
			rows = curs.fetchall()
			if rows:
				p_attr_instance1 = prediction_attributes(rows[0], type=3)
			else:
				p_attr_instance1 = None
				
			#sql_sentence2's prediction infomation is going to be displayed
			curs.execute("%s and p.p_gene_id=%s"%(sql_sentence2, p_gene_id))
			rows = curs.fetchall()
			if rows:
				p_attr_instance2 = prediction_attributes(rows[0], type=3)
				if pga_instance1 and p_attr_instance1:
					(is_accepted1, score1) = pga_instance1.prediction_accepted(p_attr_instance1.go_no, \
						[-math.log(p_attr_instance1.p_value_cut_off), p_attr_instance1.recurrence_cut_off, \
						p_attr_instance1.connectivity_cut_off, p_attr_instance1.cluster_size_cut_off, \
						p_attr_instance1.edge_gradient])
				else:
					is_accepted1, score1 = None, None
				if pga_instance2:
					(is_accepted2, score2) = pga_instance2.prediction_accepted(p_attr_instance2.go_no, \
						[-math.log(p_attr_instance2.p_value_cut_off), p_attr_instance2.recurrence_cut_off, \
						p_attr_instance2.connectivity_cut_off, p_attr_instance2.cluster_size_cut_off, \
						p_attr_instance2.edge_gradient])
				else:
					is_accepted2, score2 = None, None
				writer.writerow([p_attr_instance2.p_gene_id, p_attr_instance2.gene_no, p_attr_instance2.go_no, \
					p_attr_instance2.is_correct_lca, p_attr_instance2.avg_p_value, p_attr_instance2.recurrence_cut_off,\
					p_attr_instance2.connectivity_cut_off, p_attr_instance2.cluster_size_cut_off, p_attr_instance2.unknown_cut_off,\
					p_attr_instance2.mcl_id, p_attr_instance2.lca_list, p_attr_instance2.edge_gradient, score1, is_accepted1, \
					score2, is_accepted2])
				if not simple:
					#prepare vertex_set and edge_set to draw graphs
					vertex_set = p_attr_instance2.vertex_set[1:-1].split(',')
					vertex_set = map(int, vertex_set)
					edge_set = p_attr_instance2.edge_set[2:-2].split('},{')
					for i in range(len(edge_set)):
						edge_set[i] = edge_set[i].split(',')
						edge_set[i] = map(int, edge_set[i])
					
					#following copied from GuiAnalyzer.py
					subgraph = cluster_info_instance.graph_from_node_edge_set(vertex_set, edge_set)
					graphSrcFname = '/tmp/GuiAnalyzer.dot'
					graphFname = os.path.join(pic_output_dir, '%s_%s_%s.png'%(p_attr_instance2.p_gene_id, \
						p_attr_instance2.gene_no, p_attr_instance2.go_no))
					graphSrcF = open(graphSrcFname, 'w')
					graphDotOutput(graphSrcF, subgraph, \
						self.gene_no2gene_id, self.gene_no2go_no, \
						centralnode=p_attr_instance2.gene_no, function=p_attr_instance2.go_no, weighted=0, )
					graphSrcF.close()
					plot_type_command='neato -Goverlap=false'
					commandline = '%s -Tpng %s -o %s'%(plot_type_command, graphSrcFname, graphFname)
					system_call(commandline)
	
	def run(self):
		"""
		10-17-05
			bit control whether that setting has linear model
		"""
		schema_instance1 = form_schema_tables(self.fname1, self.acc_cutoff1, self.lm_bit1)
		schema_instance2 = form_schema_tables(self.fname2, self.acc_cutoff2, self.lm_bit2)
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		p_gene_id_set1 = p_gene_id_set_from_gene_p_table(curs, schema_instance1.gene_p_table)
		p_gene_id_set2 = p_gene_id_set_from_gene_p_table(curs, schema_instance2.gene_p_table)
		p_gene_id_set_total = p_gene_id_set_from_gene_p_table(curs, schema_instance2.p_gene_table)
		
		catI_set = p_gene_id_set1 - p_gene_id_set2
		catII_set = p_gene_id_set1 & p_gene_id_set2
		catIII_set = p_gene_id_set2 - p_gene_id_set1
		catIV_set = p_gene_id_set_total-(p_gene_id_set1|p_gene_id_set2)
		
		sample_ls_ls = []
		for p_gene_id_set in [catI_set, catII_set, catIII_set, catIV_set]:
			sample_ls_ls.append(self.sample_p_gene_id_set(p_gene_id_set, self.no_of_samples))
		
		writer = csv.writer(open(self.ofname, 'w'), delimiter = '\t')
		writer.writerow(['linear model coeffs of two settings'])
		writer.writerow([])
		writer.writerow(['No.','intercept', 'coeff1', 'coeff2', 'coeff3', 'coeff4', 'coeff5', 'intercept_p_value',\
			'coeff1_p_value', 'coeff2_p_value', 'coeff3_p_value', 'coeff4_p_value', 'coeff5_p_value',\
			'score_cut_off'])
		
		#fetch linear model coefficients
		pga_instance_list = [None, None]	#10-17-05 default is nothing, none of them have linear model
		if self.bit[0] == '1':
			pga_instance1 = p_gene_analysis()
			pga_instance1.go_no2lm_results, lm_results_2d_list = pga_instance1.get_go_no2lm_results(curs, schema_instance1.lm_table)
			pga_instance1.general_lm_results = pga_instance1.get_general_lm_results(lm_results_2d_list)
			pga_instance_list[0] = pga_instance1
			self.output_lm_model(curs, schema_instance1, writer)
		if self.bit[1] == '1':
			pga_instance2 = p_gene_analysis()
			pga_instance2.go_no2lm_results, lm_results_2d_list = pga_instance2.get_go_no2lm_results(curs, schema_instance2.lm_table)
			pga_instance2.general_lm_results = pga_instance2.get_general_lm_results(lm_results_2d_list)
			pga_instance_list[1] = pga_instance2
			self.output_lm_model(curs, schema_instance2, writer)
		
		#following is for drawing graph in output_p_gene_id_list()
		self.gene_no2gene_id = get_gene_no2gene_id(curs)
		self.gene_no2go_no = get_gene_no2go_no(curs)

		cluster_info_instance = cluster_info()
		
		for i in range(len(sample_ls_ls)):
			cat_no = i+1
			sys.stderr.write("Category %s...\n"%cat_no)
			writer.writerow(['Category %s'%cat_no])
			writer.writerow([self.category_no2information[cat_no]])
			cat_dir = 'cat%s'%cat_no
			if not os.path.isdir(cat_dir):
				os.makedirs(cat_dir)
			if i==0:	#this is different, prediction only in schema_instance1, so swap it
				self.output_p_gene_id_list(curs, schema_instance2, schema_instance1, sample_ls_ls[i], writer, cat_dir, \
					pga_instance_list[1], pga_instance_list[0], cluster_info_instance, self.simple)
			else:
				self.output_p_gene_id_list(curs, schema_instance1, schema_instance2, sample_ls_ls[i], writer, cat_dir, \
					pga_instance_list[0], pga_instance_list[1], cluster_info_instance, self.simple)
			sys.stderr.write("End Category %s.\n"%cat_no)


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:l:a:j:m:e:o:n:b:scr", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	fname1 = None
	lm_bit1 = '111'
	acc_cutoff1 = 0.6
	fname2 = None
	lm_bit2 = '111'
	acc_cutoff2 = 0.6
	ofname = None
	no_of_samples = 30
	bit = '11'
	simple = 0
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
		elif opt in ("-a"):
			acc_cutoff1 = float(arg)
		elif opt in ("-j"):
			fname2 = arg
		elif opt in ("-m"):
			lm_bit2 = arg
		elif opt in ("-e"):
			acc_cutoff2 = float(arg)
		elif opt in ("-o"):
			ofname = arg
		elif opt in ("-n"):
			no_of_samples = int(arg)
		elif opt in ("-b"):
			bit = arg
		elif opt in ("-s"):
			simple = 1
		elif opt in ("-c"):
			commit = 1
		elif opt in ("-r"):
			report = 1
	if schema and fname1 and fname2 and ofname:
		instance = SettingCmp(hostname, dbname, schema, fname1, lm_bit1, acc_cutoff1,\
			fname2, lm_bit2, acc_cutoff2, ofname, no_of_samples, bit, simple, commit, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
