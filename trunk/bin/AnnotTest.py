#!/usr/bin/env python
"""
(02-17-05)
	unittest for gene_stat_plot.
"""
import unittest
from gene_stat_plot import gene_stat

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
		self.instance.dstruc_loadin()
		
	def test_return_distinct_functions(self):
		"""
		02-17-05
		"""
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
		self.instance.debug_common_ancestor_deep_enough =1
		#it's GO:0006512
		go_no = 449
		#it's YFR053C
		gene_no = 2048
		k_functions_set = self.instance.known_genes_dict[gene_no]
		is_correct = self.instance.common_ancestor_deep_enough(go_no, k_functions_set)
		print "common_ancestor_deep_enough result for %s and gene %s: %s"%(go_no, gene_no, is_correct)
		
if __name__ == '__main__':
	suite = unittest.TestSuite()
	"""
	add one by one
	"""
	#suite.addTest(TestGeneStat("test_return_distinct_functions"))
	#suite.addTest(TestGeneStat("test_L1_match"))
	"""
	add all
	"""
	suite.addTest(unittest.makeSuite(TestGeneStat))
	unittest.TextTestRunner(verbosity=2).run(suite)
