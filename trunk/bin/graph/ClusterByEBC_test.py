#!/usr/bin/env python
"""
Usage:	ClusterByEBC_test.py size_cutoff connectivity_cutoff

A program to test the ClusterByEBC algorithm from cc_from_edge_list.cc

"""


import sys,os

if len(sys.argv)==3:
	edge_list=[(0, 5), (0, 1), (0, 6), (1, 2), (1, 3), (1, 4), (2, 3), (4, 5), (6, 8), (6, 7), (7, 8)]
	from cc_from_edge_list import ClusterByEBC
	cf_instance =ClusterByEBC(edge_list,int(sys.argv[1]),float(sys.argv[2]))
	cf_instance.run()
	print "final clusters is"
	print cf_instance.cc_vertex_list
	print cf_instance.cc_list
else:
	print __doc__
	sys.exit(2)
