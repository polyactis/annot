#!/usr/bin/env python
""""
Usage: AnnotTest.py -y TestCaseType [OPTIONS]

Option:
	-y ..., --type=...	which test case should be invoked.
	-h, --help              show this help

Examples:
	AnnotTest.py -y 2

02-17-05
	1: unittest for gene_stat_plot.

02-20-05
	2: unittest for subgraph_visualize
"""
import unittest, os, sys, getopt
from gene_stat_plot import gene_stat
from visualize.subgraph_visualize import subgraph_visualize

class TestGeneStat(unittest.TestCase):
	def setUp(self):
		hostname = 'zhoudb'
		dbname = 'graphdb'
		schema = 'sc_54'
		table = 'cluster_stat_tmp'
		mcl_table = 'mcl_result_p3g5e6d4q5n80'
		p_value_cut_off = 0.01
		connectivity_cut_off = 0.8
		recurrence_cut_off = 5
		cluster_size_cut_off = 1000
		depth_cut_off = 5
		dir_files = None
		log = 0
		judger_type = 1
		leave_one_out = 1
		wu = 1
		report = 1
		commit = 0
		unknown_cut_off = 1
		gene_table = 'p_gene'
		dominant = 0
		plottype = 3
		subgraph_cut_off = 0
		debug = 0
		accuracy_cut_off = 0
		stat_table_fname = 'null'
		self.instance = gene_stat(hostname, dbname, schema, table, mcl_table, p_value_cut_off,\
			unknown_cut_off, connectivity_cut_off, recurrence_cut_off, cluster_size_cut_off,\
			leave_one_out, wu, report, \
			log, judger_type, depth_cut_off, dir_files, commit, gene_table, dominant, plottype, \
			stat_table_fname, subgraph_cut_off, debug, accuracy_cut_off)
		
	def test_return_distinct_functions(self):
		"""
		02-17-05
		"""
		#loadin data structures used in return_distinct_functions()
		self.instance.dstruc_loadin()
		self.instance.debug_return_distinct_functions = 1
		#1029 is father of 824, 824 is father of 1030
		go_no_list = [824, 1029, 1030]
		new_go_no_list = self.instance.return_distinct_functions(go_no_list)
		print new_go_no_list
		self.assertEqual(new_go_no_list, [824])
	
	def test_L1_match(self):
		"""
		02-17-05
		"""
		#loadin data structures used in L1_match()
		self.instance.dstruc_loadin()
		self.instance.debug_L1_match =1
		#it's GO:0006512
		go_no = 449
		#it's YFR053C
		gene_no = 2048
		k_functions_set = self.instance.known_genes_dict[gene_no]
		is_correct = self.instance.L1_match(go_no, k_functions_set)
		print "L1_match result for %s and gene %s: %s"%(go_no, gene_no, is_correct)
		#it's GO:0046165, a child of GO:0006166, should be correct
		go_no = 1531
		is_correct = self.instance.L1_match(go_no, k_functions_set)
		print "L1_match result for %s and gene %s: %s"%(go_no, gene_no, is_correct)
		#it's GO:0046173, a child of GO:0046165, should be wrong
		go_no = 1532
		is_correct = self.instance.L1_match(go_no, k_functions_set)
		print "L1_match result for %s and gene %s: %s"%(go_no, gene_no, is_correct)
		
	def test_common_ancestor_deep_enough(self):
		"""
		02-17-05
		"""
		#loadin data structures used in common_ancestor_deep_enough()
		self.instance.dstruc_loadin()
		self.instance.debug_common_ancestor_deep_enough =1
		#it's GO:0006512
		go_no = 449
		#it's YFR053C
		gene_no = 2048
		k_functions_set = self.instance.known_genes_dict[gene_no]
		is_correct = self.instance.common_ancestor_deep_enough(go_no, k_functions_set)
		print "common_ancestor_deep_enough result for %s and gene %s: %s"%(go_no, gene_no, is_correct)
	
	def test_submit(self):
		"""
		02-20-05
			testing the new submit()
		"""
		#setting the testing gene_table and flag the needcommit
		self.instance.gene_table = 'p_gene_test'
		self.instance.needcommit = 1
		#construct a testing prediction_tuple2list
		tuple = (3,  4)	#(recurrence, connectivity)
		unit = [0.001, 3, 5000, 50, 0]	#[p-value, mcl_id, gene_no, go_no, is_correct]
		prediction_list = [unit]
		self.instance.prediction_tuple2list[tuple] = prediction_list
		#to get the cluster size.
		self.instance.mcl_id2vertex_set[3] = [1,2,3,4,5,6,7]
		self.instance.submit()

class TestSubgraphVisualize(unittest.TestCase):
	def setUp(self):
		hostname = 'zhoudb'
		dbname = 'graphdb'
		schema = 'sc_54'
		table = 'splat_result_1'
		mcl_table = 'mcl_result_1'
		gene_table = 'p_gene_1'
		function = 393
		functioncolor = 'green'
		centralnode = 3
		mcl_id = 0
		#context subgraph
		type = 2
		r_fname = '/tmp/AnnotTest.R'
		self.instance = subgraph_visualize(hostname, dbname, schema, table, mcl_table,\
			gene_table, function, functioncolor, centralnode, mcl_id, type, r_fname)
		self.instance.dstruc_loadin()
	
	def test_type2(self):
		"""
		02-20-05
			test the capability to draw context_subgraph
		"""
		self.instance.centralnode = 6166
		self.instance.function = 704
		subgraph = self.instance.context_subgraph(self.instance.centralnode, self.instance.function)
		self.instance.weighted_subgraph(subgraph)
		self.instance.r_f.close()
		
		
		
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "type="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hy:", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	TestCaseDict = {1:TestGeneStat,
		2: TestSubgraphVisualize}
	type = 0
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-y", "--type"):
			type = int(arg)
			
	if type:
		suite = unittest.TestSuite()
		"""
		add one by one
		"""
		#suite.addTest(TestGeneStat("test_return_distinct_functions"))
		#suite.addTest(TestGeneStat("test_L1_match"))
		#suite.addTest(TestGeneStat("test_submit"))
		"""
		add all
		"""
		suite.addTest(unittest.makeSuite(TestCaseDict[type]))
		unittest.TextTestRunner(verbosity=2).run(suite)

	else:
		print __doc__
		sys.exit(2)
