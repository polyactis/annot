#!/usr/bin/env python

import unittest
from go_node_distance import go_node_distance

class TestGoNodeDistance(unittest.TestCase):
	def setUp(self):
		hostname = 'zhoudb'
		dbname = 'graphdb'
		schema = 'sc_54'
		table = 'node_distance'
		branch = 1
		depth = 0
		all = 0
		new_table = 0
		report = 0
		commit = 0
		log = 0
		
		self.instance = go_node_distance(hostname, dbname, schema, table, branch, depth, all, new_table, report, \
			commit, log)
		self.instance.dstruc_loadin()
	
	def testprocess_2indices(self):
		self.instance.debug_process_2indices = 1
		#a plain example
		index_set1 = self.instance.go_id2index[14240]
		index_set2 = self.instance.go_id2index[14265]
		for index1 in index_set1:
			for index2 in index_set2:
				self.instance.process_2indices(index1, index2)
		#one is parent of another, and merging their path creates a loop structure
		index_set1 = self.instance.go_id2index[9613]
		index_set2 = self.instance.go_id2index[9614]
		for index1 in index_set1:
			for index2 in index_set2:
				self.instance.process_2indices(index1, index2)

if __name__ == '__main__':
	suite = unittest.TestSuite()
	suite.addTest(unittest.makeSuite(TestGoNodeDistance))
	unittest.TextTestRunner(verbosity=2).run(suite)
