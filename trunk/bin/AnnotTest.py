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

02-20-05
	3: unittest for gene_stat

02-25-05
	4: unittest for CrackSplat

"""
import unittest, os, sys, getopt
from gene_stat_plot import gene_stat
from gene_stat import gene_stat as gene_stat_slim
from visualize.subgraph_visualize import subgraph_visualize
from CrackSplat import CrackSplat

class TestGeneStatPlot(unittest.TestCase):
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
		

class TestGeneStat(unittest.TestCase):
	"""
	02-21-05
	"""
	def setUp(self):
		hostname = 'zhoudb'
		dbname = 'graphdb'
		schema = 'sc_54'
		table = 'cluster_stat_tmp'
		mcl_table = 'mcl_result_p3g5e6d4q5n80'
		depth_cut_off = 5
		dir_files = None
		leave_one_out = 1
		wu = 1
		report = 1
		commit = 0
		gene_table = 'p_gene'
		subgraph_cut_off = 0
		debug = 0
		self.instance = gene_stat_slim(hostname, dbname, schema, table, mcl_table, \
			leave_one_out, wu, report, depth_cut_off, dir_files, commit, gene_table,\
			subgraph_cut_off, debug)

	
	def test_L1_match(self):
		"""
		02-21-05
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
		02-21-05
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
		02-21-05
			testing the new submit()
		"""
		#setting the testing gene_table and flag the needcommit
		self.instance.gene_table = 'p_gene_test'
		self.instance.needcommit = 1
		#construct a testing prediction_tuple2list
		tuple = (3,  4)	#(recurrence, connectivity)
		#[p-value, mcl_id, gene_no, go_no, is_correct, is_correct_L1, \
		#is_correct_lca, cluster_size, unknown_gene_ratio]
		unit = [0.001, 3, 5000, 50, 0, 1, 1, 10, 0.25]
		prediction_list = [unit]
		self.instance.prediction_tuple2list[tuple] = prediction_list
		self.instance.submit()
	
	def test__gene_stat_leave_one_out(self):
		"""
		02-21-05
			new _gene_stat_leave_one_out() adds more fields to a prediction_list
			
		"""
		"""fake a row like [mcl_id, c.leave_one_out, c.p_value_vector, c.connectivity, \
			m.recurrence_array, m.vertex_set]"""
		row = [3, 6166, '{0.25,0.1,0.01,0.001}', 0.34, '{0.7,0.8,0.9}', '{1,2,3,4,6166}']
		#setup a go_no2depth
		self.instance.go_no2depth[3] = 10
		self.instance._gene_stat_leave_one_out(row)
		for (tuple, prediction_list) in self.instance.prediction_tuple2list.iteritems():
			self.assertEqual(tuple, (2,2))
			self.assertEqual(prediction_list, [[0.001, 3, 6166, 3, -1, -1, -1, 5, 0.25]] )

class TestCrackSplat(unittest.TestCase):
	"""
	02-25-05
	"""
	def setUp(self):
		hostname = 'zhoudb'
		dbname = 'graphdb'
		schema = 'sc_54'
		table = 'splat_result_p3g5e6d4q5n80'
		mcl_table = 'mcl_result_p3g5e6d4q5n80'
		self.instance = CrackSplat(hostname, dbname, schema, table, mcl_table, \
			report=1,needcommit=0)
	
	def test_splat2graph_dict(self):
		edge_set='{{1,2},{3,2},{3,4}}'
		real_graph_dict = {(1,2):1,
			(2,3):1,
			(3,4):1}
		graph_dict = self.instance.splat2graph_dict(edge_set)
		self.assertEqual(graph_dict, real_graph_dict)
	
	def test_graph_dict2graph(self):
		graph_dict = {(1,2):1,
			(2,3):1,
			(3,4):1}
		
		(index2no, graph) = self.instance.graph_dict2graph(graph_dict)
		print index2no
		print graph
		return (index2no, graph)
	
	def test_call_modes(self):
		infname = '/tmp/g1.matrix'
		outfname = '/tmp/g1.output'
		no_of_genes = 10
		result_code = self.instance.call_modes(infname, outfname, no_of_genes)
		print "result_code: %s"%result_code
	
	def test_parse_modes_results(self):
		"""
		in alphabetical order, this should be run after test_call_modes
		"""
		from codense.codense2db import codense2db
		from graphlib import Graph
		codense2db_instance = codense2db()
		(conn, curs) = self.instance.db_connect(self.instance.hostname, self.instance.dbname, self.instance.schema)
		#really hard to pick the mapping, must make sure the graph edges exist in the database
		index2no = {0:898, 1:993, 2:761,3:915,4:3784,5:3808,6:3971,7:5500,8:5517,9:2621}
		graph = Graph.Graph()
		graph.add_edge(0,1)
		graph.add_edge(1,2)
		graph.add_edge(2,3)
		graph.add_edge(2,5)
		graph.add_edge(4,5)
		graph.add_edge(4,6)
		graph.add_edge(7,8)
		graph.add_edge(2,9)
		ls = self.instance.parse_modes_results(1, '/tmp/g1.output',index2no, graph, codense2db_instance, curs)
		for mclResult in ls:
			print "splat_id:%s"%mclResult.splat_id
			print "connectivity: %s"%mclResult.connectivity
			print "vertex_set: %s"%(repr(mclResult.vertex_set))
			print "recurrence_array: %s"%(repr(mclResult.recurrence_array))

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
	
	TestCaseDict = {1:TestGeneStatPlot,
		2: TestSubgraphVisualize,
		3: TestGeneStat,
		4: TestCrackSplat}
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
		#suite.addTest(TestGeneStatPlot("test_return_distinct_functions"))
		#suite.addTest(TestGeneStatPlot("test_L1_match"))
		#suite.addTest(TestGeneStat("test__gene_stat_leave_one_out"))
		#suite.addTest(TestGeneStat("test_submit"))
		#suite.addTest(TestGeneStat("test_common_ancestor_deep_enough"))
		"""
		add all
		"""
		suite.addTest(unittest.makeSuite(TestCaseDict[type]))
		unittest.TextTestRunner(verbosity=2).run(suite)

	else:
		print __doc__
		sys.exit(2)
