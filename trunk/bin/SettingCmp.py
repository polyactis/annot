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
	-n ...,	no_of_samples, 30, (default)
	-c,	commit the database transaction
	-b,	debug version.
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	~/script/annot/bin/SettingCmp.py -k hs_fim_40 -i hs_fim_40m4x40rec0_8 -l 11100
		-j hs_fim_40m4x40e1s1__1_0l10 -m 11001 -o /tmp/cat.cmp

Description:
	A program to show the difference between two settings by sampling
	predictions from 4 categories.
	Category 1: only in setting1.
	Category 2: in both setting1 and setting2
	Category 3: only in setting2
	Category 4: neither
	
	prediction information, score2 and is_accepted2 are all based on setting2.
	score1 and is_accepted1 are got by applying setting1's linear model coeffs
	to setting2's prediction infomation.

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
		ofname=None, no_of_samples=30, commit=0, debug=0, report=0):
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
		self.commit = int(commit)
		self.debug = int(debug)
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
		"""
		sys.stderr.write("Sampling...")
		step = len(p_gene_id_set)/no_of_samples
		index_ls = range(0, len(p_gene_id_set), step)
		p_gene_id_list = list(p_gene_id_set)
		sys.stderr.write("End sampling.\n")
		return [p_gene_id_list[index] for index in index_ls]
	
	def output_p_gene_id_list(self, curs, schema_instance, p_gene_id_list, writer, pic_output_dir,\
		pga_instance1, pga_instance2, cluster_info_instance):
		"""
		10-15-05
			add score1 and is_accepted1
		"""
		#10-15-05 following sentence slightly different from PredictionFilterByClusterSize.py in the trailing edge_gradient
			#and d_matrix is a placeholder
		sql_sentence = "SELECT p.p_gene_id, p.gene_no, p.go_no, p.is_correct, p.is_correct_l1, \
			p.is_correct_lca, p.avg_p_value, p.no_of_clusters, p.cluster_array, p.p_value_cut_off, p.recurrence_cut_off, \
			p.connectivity_cut_off, p.cluster_size_cut_off, p.unknown_cut_off, p.depth_cut_off, p.mcl_id, p.lca_list, \
			m.vertex_set, s.edge_set, 'd_matrix', p.edge_gradient from %s p, %s s, %s m where \
			p.mcl_id=s.splat_id and p.mcl_id=m.mcl_id"%(schema_instance.p_gene_table, \
			schema_instance.splat_table, schema_instance.mcl_table)
		writer.writerow(['p_gene_id', 'gene_no', 'go_no', 'is_correct_lca', 'p_value', 'recurrence', 'connectivity',\
			'cluster_size', 'unknown_ratio', 'mcl_id', 'lca_list', 'edge_gradient', 'score1', 'is_accepted1', 'score2', 'is_accepted2'])
		for p_gene_id in p_gene_id_list:
			curs.execute("%s and p.p_gene_id=%s"%(sql_sentence, p_gene_id))
			rows = curs.fetchall()
			for row in rows:
				p_attr_instance = prediction_attributes(row)
				edge_gradient = row[-1]
				(is_accepted1, score1) = pga_instance1.prediction_accepted(p_attr_instance.go_no, \
					[-math.log(p_attr_instance.p_value_cut_off), p_attr_instance.recurrence_cut_off, \
					p_attr_instance.connectivity_cut_off, p_attr_instance.cluster_size_cut_off, edge_gradient])
				(is_accepted2, score2) = pga_instance2.prediction_accepted(p_attr_instance.go_no, \
					[-math.log(p_attr_instance.p_value_cut_off), p_attr_instance.recurrence_cut_off, \
					p_attr_instance.connectivity_cut_off, p_attr_instance.cluster_size_cut_off, edge_gradient])
				writer.writerow([p_attr_instance.p_gene_id, p_attr_instance.gene_no, p_attr_instance.go_no, \
					p_attr_instance.is_correct_lca, p_attr_instance.avg_p_value, p_attr_instance.recurrence_cut_off,\
					p_attr_instance.connectivity_cut_off, p_attr_instance.cluster_size_cut_off, p_attr_instance.unknown_cut_off,\
					p_attr_instance.mcl_id, p_attr_instance.lca_list, edge_gradient, score1, is_accepted1, score2, is_accepted2])
				#prepare vertex_set and edge_set to draw graphs
				vertex_set = p_attr_instance.vertex_set[1:-1].split(',')
				vertex_set = map(int, vertex_set)
				edge_set = p_attr_instance.edge_set[2:-2].split('},{')
				for i in range(len(edge_set)):
					edge_set[i] = edge_set[i].split(',')
					edge_set[i] = map(int, edge_set[i])
				
				#following copied from GuiAnalyzer.py
				subgraph = cluster_info_instance.graph_from_node_edge_set(vertex_set, edge_set)
				graphSrcFname = '/tmp/GuiAnalyzer.dot'
				graphFname = os.path.join(pic_output_dir, '%s_%s_%s.png'%(p_attr_instance.p_gene_id, \
					p_attr_instance.gene_no, p_attr_instance.go_no))
				graphSrcF = open(graphSrcFname, 'w')
				graphDotOutput(graphSrcF, subgraph, \
					self.gene_no2gene_id, self.gene_no2go_no, \
					centralnode=p_attr_instance.gene_no, function=p_attr_instance.go_no, weighted=0, )
				graphSrcF.close()
				plot_type_command='neato -Goverlap=false'
				commandline = '%s -Tpng %s -o %s'%(plot_type_command, graphSrcFname, graphFname)
				system_call(commandline)
	
	def run(self):
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
		
		#fetch linear model coefficients
		pga_instance1 = p_gene_analysis()
		pga_instance2 = p_gene_analysis()
		pga_instance1.go_no2lm_results, lm_results_2d_list = pga_instance1.get_go_no2lm_results(curs, schema_instance1.lm_table)
		pga_instance1.general_lm_results = pga_instance1.get_general_lm_results(lm_results_2d_list)
		pga_instance2.go_no2lm_results, lm_results_2d_list = pga_instance2.get_go_no2lm_results(curs, schema_instance2.lm_table)
		pga_instance2.general_lm_results = pga_instance2.get_general_lm_results(lm_results_2d_list)
		
		#following is for drawing graph in output_p_gene_id_list()
		self.gene_no2gene_id = get_gene_no2gene_id(curs)
		self.gene_no2go_no = get_gene_no2go_no(curs)

		cluster_info_instance = cluster_info()
		writer = csv.writer(open(self.ofname, 'w'), delimiter = '\t')
		for i in range(len(sample_ls_ls)):
			cat_no = i+1
			sys.stderr.write("Category %s...\n"%cat_no)
			writer.writerow(['Category %s'%cat_no])
			writer.writerow([self.category_no2information[cat_no]])
			cat_dir = 'cat%s'%cat_no
			if not os.path.isdir(cat_dir):
				os.makedirs(cat_dir)
			self.output_p_gene_id_list(curs, schema_instance2, sample_ls_ls[i], writer, cat_dir, \
				pga_instance1, pga_instance2, cluster_info_instance)
			sys.stderr.write("End Category %s.\n"%cat_no)


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:l:a:j:m:e:o:n:cbr", ["help", "hostname=", \
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
	commit = 0
	debug = 0
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
		elif opt in ("-c"):
			commit = 1
		elif opt in ("-b"):
			debug = 1
		elif opt in ("-r"):
			report = 1
	if schema and fname1 and fname2 and ofname:
		instance = SettingCmp(hostname, dbname, schema, fname1, lm_bit1, acc_cutoff1,\
			fname2, lm_bit2, acc_cutoff2, ofname, no_of_samples, commit, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
