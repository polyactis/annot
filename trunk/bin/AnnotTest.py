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
02-28-05
	5: unittest for module_cc.linear_model
02-28-05
	6: unittest for p_gene_lm
03-01-05
	7: unittest for p_gene_analysis
03-01-05
	8: unittest for gene_p_map_redundancy
03-02-05
	9: graph_modeling
03-04-05
	10: codense2db
03-19-05
	11: crack_by_modes
03-19-05
	12: crack_by_splat
03-31-05
	13: p_gene_factor
04-03-05
	14: connectivity2homogeneity
04-03-05
	15: TestGraphModelingGraphCC
09-18-05
	16: Test_attr_of_mt_no
"""
import unittest, os, sys, getopt, csv

class TestGeneStatPlot(unittest.TestCase):
	def setUp(self):
		from gene_stat_plot import gene_stat
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
		from visualize.subgraph_visualize import subgraph_visualize
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
	
	03-14-05
		modify because gene_stat.py changed
	03-15-05
		further change (submit()) because lca_list also submitted to database by gene_stat
	"""
	def setUp(self):
		"""
		03-14-05
			use schema sc_54_6661 to replace sc_54
		"""
		from gene_stat import gene_stat as gene_stat_slim
		hostname = 'zhoudb'
		dbname = 'graphdb'
		schema = 'sc_54_6661'
		table = 'cluster_stat_hhu_2nd'
		mcl_table = 'mcl_result_hhu_2nd'
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
			leave_one_out, wu, report, depth_cut_off, dir_files, commit, gene_table)
		
	
	def test_L1_match(self):
		"""
		02-21-05
		
		03-14-05
			
		"""
		#call dstruc_loadin()
		from codense.common import db_connect, get_go_id2go_no
		(conn, curs) = db_connect(self.instance.hostname, self.instance.dbname, self.instance.schema)
		self.instance.dstruc_loadin(curs)
		#get go_id2go_no
		go_id2go_no = get_go_id2go_no(curs)
		from gene_p_map_redundancy import gene_p_map_redundancy
		node_distance_class = gene_p_map_redundancy()
		self.instance.debug_L1_match =1
		#it's GO:0006512
		go_no = go_id2go_no['GO:0006512']
		#it's YFR053C
		gene_no = 2048
		k_functions_set = self.instance.known_genes_dict[gene_no]
		is_correct = self.instance.L1_match(go_no, k_functions_set, node_distance_class, curs)
		print "L1_match result for %s and gene %s: %s"%(go_no, gene_no, is_correct)
		#it's GO:0046165, a child of GO:0006166, should be correct
		go_no = go_id2go_no['GO:0046165']
		is_correct = self.instance.L1_match(go_no, k_functions_set, node_distance_class, curs)
		print "L1_match result for %s and gene %s: %s"%(go_no, gene_no, is_correct)
		#it's GO:0046173, a child of GO:0046165, should be wrong, after 03-14-05, using sc_54_6661, it's also correct.
		go_no = go_id2go_no['GO:0046173']
		is_correct = self.instance.L1_match(go_no, k_functions_set, node_distance_class, curs)
		print "L1_match result for %s and gene %s: %s"%(go_no, gene_no, is_correct)
		
	def test_common_ancestor_deep_enough(self):
		"""
		02-21-05
		
		03-14-05
		
		"""
		#call dstruc_loadin()
		from codense.common import db_connect
		(conn, curs) = db_connect(self.instance.hostname, self.instance.dbname, self.instance.schema)
		self.instance.dstruc_loadin(curs)
		
		from gene_p_map_redundancy import gene_p_map_redundancy
		node_distance_class = gene_p_map_redundancy()
		
		self.instance.debug_common_ancestor_deep_enough =1
		#it's GO:0006512
		go_no = 449
		#it's YFR053C
		gene_no = 2048
		k_functions_set = self.instance.known_genes_dict[gene_no]
		is_correct = self.instance.common_ancestor_deep_enough(go_no, k_functions_set, node_distance_class, curs)
		print "common_ancestor_deep_enough result for %s and gene %s: %s"%(go_no, gene_no, is_correct)
	
	def test_submit(self):
		"""
		02-21-05
			testing the new submit()
		03-15-05
			add lca_list to submit
		"""
		from codense.common import db_connect
		(conn, curs) = db_connect(self.instance.hostname, self.instance.dbname, self.instance.schema)
		#setting the testing gene_table and flag the needcommit
		self.instance.gene_table = 'p_gene_test'
		#construct a testing prediction_tuple2list
		tuple = (3,  4)	#(recurrence, connectivity)
		#[p-value, mcl_id, gene_no, go_no, is_correct, is_correct_L1, \
		#is_correct_lca, cluster_size, unknown_gene_ratio, lca_list]
		unit = [0.001, 3, 5000, 50, 0, 1, 1, 10, 0.25, []]
		unit1 = [0.001, 3, 5000, 50, 0, 1, 1, 10, 0.25, [12,23]]
		prediction_list = [unit, unit1]
		self.instance.prediction_tuple2list[tuple] = prediction_list
		self.instance.submit(curs, self.instance.gene_table)
	
	def test__gene_stat_leave_one_out(self):
		"""
		02-21-05
			new _gene_stat_leave_one_out() adds more fields to a prediction_list
		
		03-14-05
		"""
		#call dstruc_loadin()
		from codense.common import db_connect
		(conn, curs) = db_connect(self.instance.hostname, self.instance.dbname, self.instance.schema)
		self.instance.dstruc_loadin(curs)
		from gene_p_map_redundancy import gene_p_map_redundancy
		node_distance_class = gene_p_map_redundancy()
		"""fake a row like [mcl_id, c.leave_one_out, c.p_value_vector, c.connectivity, \
			m.recurrence_array, m.vertex_set]"""
		row = [3, 6166, '{0.25,0.1,0.01,0.001}', 0.34, '{0.7,0.8,0.9}', '{1,2,3,4,6166}']
		#setup a go_no2depth
		self.instance.go_no2depth[3] = 10
		self.instance._gene_stat_leave_one_out(row, node_distance_class, curs)
		for (tuple, prediction_list) in self.instance.prediction_tuple2list.iteritems():
			print "prediction_space is %s"%repr(tuple)
			print "prediction_list is %s"%repr(prediction_list)

class TestCrackSplat(unittest.TestCase):
	"""
	02-25-05
	"""
	def setUp(self):
		from CrackSplat import CrackSplat
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
	

class TestLinearModel(unittest.TestCase):
	"""
	02-28-05
	"""
	def setUp(self):
		from module_cc.linear_model import linear_model
		self.instance = linear_model()
		
	def test_splat2graph_dict(self):
		x_2d_list = [[1, 14, 2],
					[1, 16, 0],
					[1, 16, 2],
					[1, 18, 0],
					[1, 18, 2],
					[1, 18, 4],
					[1, 20, 0],
					[1, 20, 2],
					[1, 20, 4],
					[1, 20, 6],
					[1, 22, 0],
					[1, 22, 2],
					[1, 22, 4],
					[1, 22, 6],
					[1, 22, 8],
					[1, 24, 2],
					[1, 26, 0],
					[1, 26, 2],
					[1, 26, 4]]
		y_list = [0.5500000, 0.6842105, 0.7647059, 0.5030303, 0.4000000, 0.7647059, \
			0.5439739, 0.5727412, 0.5603448, 0.5625000, 0.5781991, 0.6218905, \
			0.561538, 0.5573770, 0.5208333, 0.6666667, 0.8666667, 0.6862745, 0.6551724]

		self.instance.prepare_data(y_list, x_2d_list);
		self.instance.run()
		coeff_list = self.instance.coefficients()
		chisq_tuple = self.instance.chisq_return()
		self.assertEqual
		real_coeff_list = [0.44622090360304795, 0.0096039858657607172, -0.01244257855672223]
		real_chisq_tuple = (0.18409639558598592,)
		self.assertEqual(coeff_list, real_coeff_list)
		self.assertEqual(chisq_tuple, real_chisq_tuple)
		"""
		print
		print "Coefficient list: %s"%repr(self.instance.coefficients())
		print "chisq: %s"%repr(self.instance.chisq_return())
		"""
	
	def tearDown(self):
		self.instance.cleanup()

class TestPGeneLm(unittest.TestCase):
	"""
	02-28-05
	
	03-27-05
		p_gene_lm changed
	"""
	def setUp(self):
		hostname = 'zhoudb'
		dbname = 'graphdb'
		schema = ''
		table = None
		lm_table = None
		accuracy_cut_off = 0.5
		judger_type = 0
		min_data_points = 20
		commit = 0
		report = 0
		debug = 0
		from p_gene_lm import p_gene_lm
		self.instance = p_gene_lm(hostname, dbname, schema, table, lm_table, accuracy_cut_off,\
			judger_type, min_data_points, commit, report, debug)
		self.instance.init()
	"""
	def test_p_value_outof_accuracy_cut_off(self):
		self.instance.debug = 0
		tuple_list = [[0.001, 1], [0.001, -1], [0.001, 1], [0.001,1], [0.01, 1], [0.01, 0], [0.01, 1], [0.01, -1]]
		accuracy_cut_off = 0.84
		p_value_cut_off = self.instance.p_value_outof_accuracy_cut_off(tuple_list, accuracy_cut_off)
		self.assertEqual(p_value_cut_off, 0.001)
	"""
	
	def test_data_fetch(self):
		from codense.common import db_connect
		self.instance.schema = 'sc_54'
		self.instance.table = 'p_gene_repos_2_e5'
		(conn, curs) = db_connect(self.instance.hostname, self.instance.dbname, self.instance.schema)
		self.instance.data_fetch(curs, self.instance.table)
		#print self.instance.go_no2prediction_space
	
	"""
	def test_data_prepare(self):
		from p_gene_analysis import prediction_space_attr
		prediction_space2attr = {}
		connectivity_list = [0,2,4]*3
		recurrence_list = range(2,20,2)
		for i in range(len(connectivity_list)):
			prediction_space2attr[(recurrence_list[i], connectivity_list[i])] = prediction_space_attr()
			prediction_space2attr[(recurrence_list[i], connectivity_list[i])].known_predictions = 5*i
			prediction_space2attr[(recurrence_list[i], connectivity_list[i])].p_value_cut_off = 0.001*i
		
		y_list, x_2d_list = self.instance.data_prepare(prediction_space2attr)
		print y_list
		print x_2d_list
	"""
	
	def test_lm_fit(self):
		"""
		03-27-05
			the testing data is a piece from a real p_gene table
			format: [-lg(p_value), recurrence, connectivity, is_correct]
		"""
		go_no2prediction_space = {-1:[]}
		go_no2prediction_space[-1].append([0.3904164, 9.466667, 0.034783, 0])
		go_no2prediction_space[-1].append([8.044069435, 7.820000, 0.017834, 1])
		go_no2prediction_space[-1].append([0.005514175, 8.095238, 0.019373, 0])
		go_no2prediction_space[-1].append([3.734169750, 8.843137, 0.023794, 0])
		go_no2prediction_space[-1].append([7.195437351, 7.843750, 0.013958, 1])
		go_no2prediction_space[-1].append([1.497909955, 7.222222, 0.009553, 0])
		go_no2prediction_space[-1].append([4.216512196, 7.693878, 0.021872, 0])
		go_no2lm_results = self.instance.lm_fit(None, go_no2prediction_space)
		print "\nlinear model fitting results: %s"%repr(go_no2lm_results)
		
	def test_return_score_cut_off(self):
		"""
		03-27-05
			The lm model is got from the real p_gene_table, same as above.(via R)
		"""
		self.instance.debug=1
		score_list = [[6,1],[3,0],[3,0],[4,1],[3,0],[4,0],[4,0],[5,1],[6,1],[5,1],[5,0],[6,1]]
		score_cut_off = self.instance.return_score_cut_off(score_list, 0.7, -1)
		print "\nthe score cutoff is %s"%repr(score_cut_off)

class TestPGeneAnalysis(unittest.TestCase):
	"""
	02-28-05
	
	03-06-05
		add more testing functions
	"""
	def setUp(self):
		from p_gene_analysis import p_gene_analysis
		hostname = 'zhoudb'
		dbname = 'graphdb'
		schema = 'sc_54'
		table = 'cluster_stat_repos_2'
		mcl_table = 'mcl_result_repos_2'
		p_value_cut_off = 0
		judger_type = 1
		report = 0
		commit = 0
		gene_table = 'p_gene_repos_2_e5'
		lm_table = 'lm_p_gene_repos_2_e5_v40'
		self.instance = p_gene_analysis(hostname, dbname, schema, table, mcl_table, p_value_cut_off,\
			report, judger_type, commit, gene_table, lm_table)
			
	def test_lm(self):
		"""
		03-27-05
			test the new prediction_accepted()
		"""
		import math
		from codense.common import db_connect
		(conn, curs) = db_connect(self.instance.hostname, self.instance.dbname, 'sc_54_6661')
		self.instance.go_no2lm_results, lm_results_2d_list = self.instance.get_go_no2lm_results(curs,\
			'p_gene_lm_tmp_2')
		self.instance.general_lm_results = self.instance.get_general_lm_results(lm_results_2d_list)
		print self.instance.go_no2lm_results
		print self.instance.general_lm_results
		print self.instance.prediction_accepted(56, [-math.log(0.00001),15,0.03])
		
	def test_return_go_no_map(self):
		"""
		03-06-05
		"""
		from codense.common import db_connect
		(conn, curs) = db_connect(self.instance.hostname, self.instance.dbname, self.instance.schema)
		#1029 is father of 824, 824 is father of 1030
		go_no_list = [824, 1029, 1030]
		go_no_map = self.instance.return_go_no_map(go_no_list, curs, 'go.node_dist_sc_54')
		print go_no_map
	
	def test_dict_map2group(self):
		"""
		03-06-05
		"""
		from codense.common import db_connect
		(conn, curs) = db_connect(self.instance.hostname, self.instance.dbname, self.instance.schema)
		#1029 is father of 824, 824 is father of 1030
		go_no_list = [824, 1029, 1030]
		go_no_map = self.instance.return_go_no_map(go_no_list, curs, 'go.node_dist_sc_54')
		go_no_groups = self.instance.dict_map2group(go_no_map)
		print go_no_groups
	
	def test_return_data_of_the_group(self):
		"""
		03-06-05
		"""
		from p_gene_analysis import prediction_space_attr
		#setup a go_no2prediction_space
		go_no2prediction_space = {}
		go_no2prediction_space[1] = {}
		go_no2prediction_space[1][(3,2)] = prediction_space_attr()
		unit = go_no2prediction_space[1][(3,2)]
		unit.correct_predictions = 2
		unit.known_predictions = 3
		unit.unknown_predictions = 5
		unit.tuple_list.append([0.001, 1])
		
		go_no2prediction_space[2] = {}
		go_no2prediction_space[2][(4,3)] = prediction_space_attr()
		unit = go_no2prediction_space[2][(4,3)]
		unit.correct_predictions = 3
		unit.known_predictions = 4
		unit.unknown_predictions = 6
		unit.tuple_list.append([0.05, 0])
		
		go_no_list = [1,2]
		new = self.instance.return_data_of_the_group(go_no_list, go_no2prediction_space)
		
		#output
		writer = csv.writer(sys.stdout, delimiter='\t')
		self.instance.table_output(writer, new)
		#self.output_prediction_space2attr(new)
		
	def test_return_cumulative_prediction_space2attr(self):
		"""
		03-06-05
		"""
		from p_gene_analysis import prediction_space_attr
		prediction_space2attr = {}
		prediction_space2attr[(3,2)] = prediction_space_attr()
		unit = prediction_space2attr[(3,2)]
		unit.correct_predictions = 2
		unit.known_predictions = 3
		unit.unknown_predictions = 5
		unit.tuple_list.append([0.001, 1])
		
		prediction_space2attr[(4,3)] = prediction_space_attr()
		unit = prediction_space2attr[(4,3)]
		unit.correct_predictions = 3
		unit.known_predictions = 4
		unit.unknown_predictions = 6
		unit.tuple_list.append([0.05, 0])
		
		new = self.instance.return_cumulative_prediction_space2attr(prediction_space2attr)
		#output
		writer = csv.writer(sys.stdout, delimiter='\t')
		self.instance.table_output(writer, new)
		#self.output_prediction_space2attr(new)
		
	
	def output_prediction_space2attr(self, prediction_space2attr):
		"""
		03-06-05
			output a prediction_space2attr
		"""
		for (prediction_space,unit) in prediction_space2attr.iteritems():
			print "prediction_space: %s"%repr(prediction_space)
			print "\tcorrect_predictions: %s"%unit.correct_predictions
			print "\tknown_predictions: %s"%unit.known_predictions
			print "\tunknown_predictions: %s"%unit.unknown_predictions
			print "\ttuple_list: %s"%repr(unit.tuple_list)
		
		
		
class TestGenePMapRedundancy(unittest.TestCase):
	"""
	02-28-05

	03-03-05:
		modify the example in test__p_gene_map to show the bug
	03-04-05
		the source module changed a lot, change accordingly.
	"""
	def setUp(self):
		from gene_p_map_redundancy import gene_p_map_redundancy
		hostname = 'zhoudb'
		dbname = 'graphdb'
		schema = 'sc_54'
		p_gene_table = 'p_gene_repos_2_e5'
		gene_p_table = 'gene_p_repos_2_e5'
		
		self.instance = gene_p_map_redundancy(hostname, dbname, schema,\
			p_gene_table, gene_p_table)
		#prepare a cursor
		from codense.common import db_connect, get_go_no2term_id
		(conn, curs) = db_connect(hostname, dbname, schema)
		self.instance.go_no2term_id = get_go_no2term_id(curs, self.instance.schema, self.instance.term_table)
	
	def test__p_gene_map(self):
		from codense.common import db_connect
		(conn, curs) = db_connect(self.instance.hostname, self.instance.dbname, self.instance.schema)
		#set the debug flag first
		self.instance.debug = 1
		#1029 is father of 824, 824 is father of 1030
		p_gene_id2go_no = {1:[824,0.01], 2:[1029,0.001], 3:[1030,0.05], 4:[824,0.002], 5:[1029,0.003], 6:[1029,0.01]}
		self.instance._p_gene_map(p_gene_id2go_no, self.instance.p_gene_id_map, curs,\
			self.instance.distance_table, self.instance.go_no2distance, self.instance.go_no2term_id)
		print self.instance.p_gene_id_map
	
	def test_get_distance(self):
		from codense.common import db_connect
		(conn, curs) = db_connect(self.instance.hostname, self.instance.dbname, self.instance.schema)
		#1029 is parent of 824
		jasmine_distance = self.instance.get_distance(curs, 824, 1029, self.instance.distance_table,\
			self.instance.go_no2distance, self.instance.go_no2term_id)
		self.assertEqual(jasmine_distance, 0)
	
	def test_submit(self):
		from codense.common import db_connect
		(conn, curs) = db_connect(self.instance.hostname, self.instance.dbname, self.instance.schema)
		p_gene_id_map = {33:[32,0.01]}
		self.instance.submit(curs, self.instance.gene_p_table, p_gene_id_map)
	
	def test_gene_no2p_gene_setup(self):
		row = [2, 6661, 365, 0.01]
		self.instance.gene_no2p_gene_setup(row)
		print self.instance.gene_no2p_gene
	
class TestGraphModeling(unittest.TestCase):
	"""
	03-02-05
	"""
	def setUp(self):
		from graph import graph_modeling
		self.instance = graph_modeling
		self.instance.cor_cut_off_vector_construct(0.01, 0.6)
		
	def test_ind_min_cor(self):
		"""
		03-02-05
		"""
		v1 = [0.1,0.2,0.3,0.4,0.6, 0.7, 0.8, 0.9]
		v2 = [0.2, 0.3, 0.35, 0.4, 0.6, 0.7, 0.8, 0.9]
		data = self.instance.ind_min_cor(v1, v2)
		self.assertEqual(data.value, 0.99447262287139893)
		self.assertEqual(data.degree, 5)
		self.assertEqual(data.significance, 1)
		print
		print "min jackknife correlation: %s"%data.value
		print "degree: %s"%data.degree
		print "significance: %s"%data.significance
	
	def test_ind_cor(self):
		"""
		03-02-05
		"""
		v1 = [0.1,0.2,0.3,0.4,0.6, 0.7, 0.8, 0.9]
		v2 = [0.2, 0.3, 0.35, 0.5, 0.7, 0.8, 0.9, 1.0]
		data = self.instance.ind_cor(v1,v2, 2)
		self.assertEqual(data.value, 1.0)
		self.assertEqual(data.degree, 5)
		self.assertEqual(data.significance, 0)
		print 
		print "correlation: %s"%data.value
		print "degree: %s"%data.degree
		print "significance: %s"%data.significance
	
class TestCodense2db(unittest.TestCase):
	"""
	03-04-05
	"""
	def setUp(self):
		from codense.codense2db import codense2db
		hostname = 'zhoudb'
		dbname = 'graphdb'
		schema = 'sc_54'
		table = 'splat_result'
		mcl_table = 'mcl_result'
		mapping_file = "whatever"
		input_file = 'whatever'
		self.instance = codense2db(input_file, hostname, dbname, schema, table, mcl_table, mapping_file)

	
	def test_codense_parser(self):
		#get the cursor
		from codense.common import db_connect, get_gene_id2gene_no
		(conn, curs) = db_connect(self.instance.hostname, self.instance.dbname, self.instance.schema)
		
		gene_id2gene_no = get_gene_id2gene_no(curs)
		
		input_file = raw_input("the path to the input_file to test codense_parser:")
		inf = csv.reader(open(input_file, 'r'), delimiter='\t')
		row = inf.next()
		cluster = self.instance.codense_parser(row, gene_id2gene_no, curs)
		self.cluster_output(cluster)
		
	def test_copath_parser(self):
		#get the cursor
		from codense.common import db_connect, get_haiyan_no2gene_no
		(conn, curs) = db_connect(self.instance.hostname, self.instance.dbname, self.instance.schema)
		
		mapping_file = raw_input("the path to the mapping file:")
		haiyan_no2gene_no = get_haiyan_no2gene_no(mapping_file)
		
		input_file = raw_input("the path to the input_file to test copath_parser:")
		inf = csv.reader(open(input_file, 'r'), delimiter='\t')
		row = inf.next()
		cluster = self.instance.copath_parser(row, haiyan_no2gene_no, curs)
		self.cluster_output(cluster)
	
	def cluster_output(self, cluster):
		print "cluster_id: %s"%cluster.cluster_id
		print "splat_connectivity: %s"%cluster.splat_connectivity
		print "connectivity: %s"%cluster.connectivity
		print "no_of_edges: %s"%cluster.no_of_edges
		print "no_of_nodes: %s"%len(cluster.vertex_set)
		print "vertex_set: %s"%repr(cluster.vertex_set)
		print "edge_set: %s"%repr(cluster.edge_set)
		print "recurrence_array: %s"%repr(cluster.recurrence_array)
				
class TestCrackByModes(unittest.TestCase):
	"""
	03-19-05
		
	"""
	def setUp(self):
		from CrackSplat import crack_by_modes
		self.instance = crack_by_modes()
	
	def test_graph_dict2graph(self):
		graph_dict = {(1,2):1,
			(2,3):1,
			(3,4):1}
		
		(index2no, graph) = self.instance.graph_dict2graph(graph_dict)
		print index2no
		print graph
		return (index2no, graph)
	
	def test_call_modes(self):
		infname = os.path.join(os.path.expanduser('~'),'script/hhu_clustering/data/input/g1.matrix') 
		outfname = os.path.join(os.path.expanduser('~'),'script/hhu_clustering/data/output/g1.output') 
		no_of_genes = 10
		result_code = self.instance.call_modes(infname, outfname, no_of_genes)
		print "result_code: %s"%result_code
	
	def test_parse_modes_results(self):
		"""
		in alphabetical order, this should be run after test_call_modes
		"""
		from codense.codense2db import codense2db
		from graphlib import Graph
		from codense.common import db_connect
		hostname = 'zhoudb'
		dbname = 'graphdb'
		schema = 'sc_54'
		(conn, curs) = db_connect(hostname, dbname, schema)
		codense2db_instance = codense2db()
		codense2db_instance.debug = 1
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
		infname = os.path.join(os.path.expanduser('~'),'script/hhu_clustering/data/output/g1.output') 
		ls = self.instance.parse_modes_results(1, infname,index2no, graph, codense2db_instance, curs)
		for mclResult in ls:
			print "splat_id:%s"%mclResult.splat_id
			print "connectivity: %s"%mclResult.connectivity
			print "vertex_set: %s"%(repr(mclResult.vertex_set))
			print "edge_set: %s"%(repr(mclResult.edge_set))
			print "recurrence_array: %s"%(repr(mclResult.recurrence_array))

class TestCrackBySplat(unittest.TestCase):
	"""
	03-19-05
		
	"""
	def setUp(self):
		from CrackSplat import crack_by_splat
		self.instance = crack_by_splat(debug=1)
		#go into following directory first
		dir_files = '/tmp/yh'
		if not os.path.isdir(dir_files):
			os.makedirs(dir_files)
		else:
			sys.stderr.write("Warning, directory %s already exists.\n"%(dir_files))
		os.chdir(dir_files)
	
	def test_write_splat_input_files(self):
		edge_set = [(898,993), (761,993), (761,915), (761,3808), (3784,3808), (3784,3971), (5500,5517), (761,2621)]
		from codense.common import db_connect
		hostname = 'zhoudb'
		dbname = 'graphdb'
		schema = 'sc_54'
		(conn, curs) = db_connect(hostname, dbname, schema)
		from codense.codense2db import codense2db
		codense2db_instance = codense2db()
		combined_cor_vector, combined_sig_vector = codense2db_instance.get_combined_cor_vector(curs, edge_set)
		no_of_datasets = self.instance.write_splat_input_files(edge_set, combined_sig_vector)
		print "no_of_datasets: %s"%no_of_datasets
		return no_of_datasets

	
	def test_call_splat(self):
		sys.stderr.write("call write_splat_input_files() first...")
		no_of_datasets = self.test_write_splat_input_files()
		sys.stderr.write("done\n")
		result_code = self.instance.call_splat('gph', no_of_datasets, 3, min_cut=1)
		print "result_code: %s"%result_code
	
	def test_parse_splat_results(self):
		"""
		in alphabetical order, this should be run after test_call_splat
		"""
		from splat_to_db import splat_to_db
		from codense.codense2db import codense2db
		from graphlib import Graph
		from codense.common import db_connect
		splat_to_db_instance = splat_to_db()
		hostname = 'zhoudb'
		dbname = 'graphdb'
		schema = 'sc_54'
		codense2db_instance = codense2db()
		(conn, curs) = db_connect(hostname, dbname, schema)
		ls = self.instance.parse_splat_results(1, 'patterns-splat', splat_to_db_instance, codense2db_instance, curs)
		for mclResult in ls:
			print "splat_id:%s"%mclResult.splat_id
			print "connectivity: %s"%mclResult.connectivity
			print "vertex_set: %s"%(repr(mclResult.vertex_set))
			print "edge_set: %s"%(repr(mclResult.edge_set))
			print "recurrence_array: %s"%(repr(mclResult.recurrence_array))
	
class TestPGeneFactor(unittest.TestCase):
	"""
	03-31-05
		
	"""
	def setUp(self):
		from p_gene_factor import p_gene_factor
		self.instance = p_gene_factor(debug=1)
		
	def test_group_data(self):
		list_2d = [[2,3],[1.2,1],[1.3,0],[2,5],[3,6],[3,4],[4,2],[5,7]]
		key2data = self.instance.group_data(list_2d, key_column=0, no_of_groups=3)
		print key2data

class TestConnectivity2Homogeneity(unittest.TestCase):
	"""
	04-03-05
		
	"""
	def setUp(self):
		"""
		04-03-05
		"""
		from connectivity2homogeneity import connectivity2homogeneity
		schema = 'sc_54_6661'
		homedir = os.path.expanduser('~')
		input_file = os.path.join(homedir, 'bin/hhu_clustering/data/input/sc_54_6661_7')
		self.instance = connectivity2homogeneity(input_file=input_file, schema=schema, \
			output_fname='/tmp/data', debug=1, random_subgraph_size=6)

	def test_get_summary_graph(self):
		summary_graph = self.instance.get_summary_graph(self.instance.input_file)
		print "The summary_graph is %s"%repr(summary_graph)

	def test_get_gene_no2go_no_list(self):
		"""
		04-03-05
		"""
		from codense.common import db_connect
		(conn, curs) = db_connect(self.instance.hostname, self.instance.dbname, self.instance.schema)
		self.instance.debug = 0
		gene_no2go_no_list = self.instance.get_gene_no2go_no_list(curs)
		print "The length of gene_no2go_no_list is %s"%len(gene_no2go_no_list)
		return gene_no2go_no_list
	
	def test__connectivity2homogeneity(self):
		"""
		04-03-05
		"""
		from graphlib import Graph
		subgraph = Graph.Graph()
		#testing example from mcl_result_copathsc_54_6661_merge_6g5e6d40q20s80c50so of sc_54_6661
		subgraph.add_edge(1215,5616,5)
		subgraph.add_edge(1215,6325,5)
		subgraph.add_edge(4240,6325,5)
		subgraph.add_edge(5135,6325,5)
		subgraph.add_edge(1428,5135,5)
		subgraph_id =1
		node_list = subgraph.node_list()
		edge_list = subgraph.edge_list()
		edge_set = map(subgraph.edge_by_id, edge_list)
		gene_no2go_no_list = self.test_get_gene_no2go_no_list()
		writer = csv.writer(sys.stdout, delimiter='\t')
		self.instance.debug = 1
		self.instance._connectivity2homogeneity(subgraph_id, node_list, edge_set, gene_no2go_no_list,writer)
		

class TestGraphModelingGraphCC(unittest.TestCase):
	"""
	04-03-05
		this is not a unittest for other classes. Rather, it's a class
		to check whether summary_graph generated by graph.cc and
		the sig_vector generated by complete_cor_vector.py(backed
		by graph_modeling.cc) are have the same edges given a 
		minimum weight.
	"""
	def setUp(self):
		"""
		04-03-05
		"""
		self.hostname = 'zhoudb'
		self.dbname = 'graphdb'
		self.schema = 'sc_54_6661'
		homedir = os.path.expanduser('~')
		self.input_file = os.path.join(homedir, 'bin/hhu_clustering/data/input/sc_54_6661_5')
	
	def test_compare_edge_weights(self):
		"""
		04-03-05
		"""
		from codense.common import db_connect
		(conn, curs) = db_connect(self.hostname, self.dbname, self.schema)
		from connectivity_original import connectivity_original
		conn_ori_instance = connectivity_original()
		
		reader = csv.reader(open(self.input_file,'r'), delimiter=' ')
		no_of_edges = 0
		for row in reader:
			if row[0] == 'e':
				no_of_edges += 1
				gene_no1 = int(row[1])
				gene_no2 = int(row[2])
				edge = [gene_no1, gene_no2]
				edge.sort()
				weight = int(row[3])
				edge_weight = conn_ori_instance._get_edge(curs, edge, edge_table='edge_cor_vector')
				if weight!=edge_weight:
					sys.stderr.write("edge: %s, summary weight: %s, database weight: %s.\n"\
						%(repr(edge), weight, edge_weight))
				sys.stderr.write('%s%s'%('\x08'*10, no_of_edges))

class Test_attr_of_mt_no(unittest.TestCase):
	"""
	09-18-05
	"""
	
	def test_attr_of_mt_no(self):
		"""
		04-03-05
		"""
		ls = [
		[0.847, 'BC073913', 1],
		[0.858, 'BC073912', 2],
		[0.883, 'BC073911', 3],
		[0.925, 'BC073910', 4],
		[0.925, 'BC073913', 5],
		[0.847, 'BC073912', 6],
		[0.858, 'BC073911', 7],
		[0.883, 'BC073910', 8],
		[0.858, 'BC073913', 9],
		[0.883, 'BC073912', 10],
		[0.925, 'BC073911', 11],
		[0.847, 'BC073910', 12],
		[0.883, 'BC073913', 13],
		[0.925, 'BC073912', 14],
		[0.847, 'BC073911', 15],
		[0.858, 'BC073910', 16],
		]
		from binding_site2gene_id2mt_no import attr_of_mt_no
		unit = attr_of_mt_no(4, 1)	#top_number is 3, and debug is enabled.
		for matrix_similarity_score, prom_acc, id in ls:
			unit.consume_new_row(matrix_similarity_score, prom_acc, id)
		from heapq import heappush, heappop
		
		while len(unit.hq_list)>0:
			row = heappop(unit.hq_list)
			print row
		print unit.acc2id_set
		
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
		4: TestCrackSplat,
		5: TestLinearModel,
		6: TestPGeneLm,
		7: TestPGeneAnalysis,
		8: TestGenePMapRedundancy,
		9: TestGraphModeling,
		10: TestCodense2db,
		11: TestCrackByModes,
		12: TestCrackBySplat,
		13: TestPGeneFactor,
		14: TestConnectivity2Homogeneity,
		15: TestGraphModelingGraphCC,
		16: Test_attr_of_mt_no}
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
