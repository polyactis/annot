#!/usr/bin/env python
"""
Usage: DrawPredGOExptTFCompTF_Patterns.py -k SCHEMA -i xxx -g xx -c xxx -e xxx -p xxx
	-o xxx [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-i ..., p_gene_table
	-g ...,	gene_p_table
	-n ...,	prot_interaction_table, 'mrinal_pi.intact_interaction'(default)
	-c ...,	comp_cluster_bs_table
	-e ...,	expt_cluster_bs_table
	-p ...,	pattern_table = 
	-o ...,	output_dir = 
	-m ...,	comp_tf_mapping_table = 'graph.gene_id2mt_no'
	-t ...,	expt_tf_mapping_table = 'graph.tf_mapping'
	-x ...,	tax_id, 9606(default)
	-b,	debug version.
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	python2.3 ./DrawPredGOExptTFCompTF_Patterns.py -k hs_fim_65 -c cluster_bs_hs_fim_65_m5x65s4l5e0p001geneid \
		-e cluster_bs_hs_fim_65_m5x65s4l5e0p001expt2 -p pattern_hs_fim_65_m5x65s4l5 -o /tmp/tf_go_patterns/ -i p_gene_hs_fim_65_n2s175_m5x65s4l5_ft2_e5
		-g gene_p_hs_fim_65_n2s175_m5x65s4l5_ft2_e5_000001a60
	
Description:
	Draw GO function graph, Experimental TF graph, Computational TF graph.
		(depends on whether they exist or not)
	color notation:
	 	green: standout but not associated
		yellow: standout and associated
		red: associated
		blue: not associated and not associated
	Edges colored red in GO function graph denote protein interaction.
		They might not be part of the pattern (compare it with TF graph).
	
	For the augmented PI graph:
		nodes and edges from the coexpr pattern are marked in 'green' and 'red' color, respectively
		nodes and edges added by protein interaction are marked in 'red' and 'gray' color, respectively
		overlapping edges are widened with 'red' color
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import sys, os, getopt, csv, math
from codense.common import db_connect, get_gene_id2gene_symbol, get_go_id2name, pg_1d_array2python_ls
from sets import Set
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import pylab
import Numeric

"""
2006-11-19
	three functions below used to investigate the extent of difference among the condition's and vertex_set's based on which
		a gene is projected into different function
	ex:
		compare_gene_condition_vertex_set(curs, 'p_gene_hs_fim_65_n2s200_m5x65s4l5_ft2_e5',\
			'gene_p_hs_fim_65_n2s200_m5x65s4l5_ft2_e5_000001a60', 'good_cl_hs_fim_65_n2s200_m5x65s4l5_ft2_e5_000001a60',\
			'tmp/gene_condition_vertex_set.hs_fim_65_n2s200_m5x65s4l5_ft2.cmp')
	
"""


def get_fisher_p_value_of_two_lists(recurrence_array1, recurrence_array2):
	import rpy
	contigency_table = [[0,0], [0,0]]
	for i in range(len(recurrence_array1)):
		contigency_table[recurrence_array1[i]][recurrence_array2[i]] += 1
	flat_list = [contigency_table[0][0], contigency_table[0][1], contigency_table[1][0], contigency_table[1][1]]
	contigency_matrix = rpy.r.matrix(flat_list, 2,2, byrow=rpy.r.TRUE)
	p_value = rpy.r.fisher_test(contigency_matrix)['p.value']
	return [contigency_table, p_value]


def get_mcl_id_sharing(mcl_id_go_no_list, mcl_id2recurrence_array_vertex_set):
	no_of_mcl_ids = len(mcl_id_go_no_list)
	result = []
	if no_of_mcl_ids>1:
		for i in range(no_of_mcl_ids):
			for j in range(i+1, no_of_mcl_ids):
				mcl_id1, go_no1 = mcl_id_go_no_list[i]
				mcl_id2, go_no2 = mcl_id_go_no_list[j]
				recurrence_array1, vertex_set1 = mcl_id2recurrence_array_vertex_set[mcl_id1]
				recurrence_array2, vertex_set2 = mcl_id2recurrence_array_vertex_set[mcl_id2]
				contigency_table, p_value = get_fisher_p_value_of_two_lists(recurrence_array1, recurrence_array2)
				vertex_sharing_perc = float(len(vertex_set1&vertex_set2))/len(vertex_set1|vertex_set2)
				result.append([mcl_id1, go_no1, mcl_id2, go_no2, contigency_table, p_value, vertex_sharing_perc])
	return result


def compare_gene_condition_vertex_set(curs, p_gene_table, gene_p_table, good_cluster_table, output_fname):
	import os, sys, csv
	from sets import Set
	from codense.common import pg_1d_array2python_ls
	sys.stderr.write("Getting gene_no2mcl_id_go_no_list ...\n")
	gene_no2mcl_id_go_no_list = {}
	curs.execute("DECLARE crs1 CURSOR for select p.gene_no, p.mcl_id, p.go_no from %s p, %s g\
		where p.p_gene_id=g.p_gene_id"%(p_gene_table, gene_p_table))
	counter = 0
	curs.execute("fetch 1000 from crs1")
	rows = curs.fetchall()
	mcl_id_set = Set()
	while rows:
		for row in rows:
			gene_no, mcl_id, go_no = row
			mcl_id_set.add(mcl_id)
			if gene_no not in gene_no2mcl_id_go_no_list:
				gene_no2mcl_id_go_no_list[gene_no] = []
			gene_no2mcl_id_go_no_list[gene_no].append([mcl_id, go_no])
			counter += 1
		sys.stderr.write("%s%s"%('\x08'*30, counter))
		curs.execute("fetch 1000 from crs1")
		rows = curs.fetchall()
	curs.execute("close crs1")
	sys.stderr.write("Done.\n")
	
	sys.stderr.write("Getting mcl_id2recurrence_array_vertex_set ...\n")
	mcl_id2recurrence_array_vertex_set = {}
	curs.execute("DECLARE crs0 CURSOR for select mcl_id, vertex_set, recurrence_array from %s"%good_cluster_table)
	counter = 0
	curs.execute("fetch 1000 from crs0")
	rows = curs.fetchall()
	while rows:
		for row in rows:
			mcl_id, vertex_set, recurrence_array = row
			vertex_set = pg_1d_array2python_ls(vertex_set)
			recurrence_array = pg_1d_array2python_ls(recurrence_array)
			mcl_id2recurrence_array_vertex_set[mcl_id] =	[recurrence_array, Set(vertex_set)]
			counter += 1
		sys.stderr.write("%s%s"%('\x08'*30, counter))
		curs.execute("fetch 1000 from crs0")
		rows = curs.fetchall()
	curs.execute("close crs0")
	sys.stderr.write("Done.\n")
	
	sys.stderr.write("Comparing gene condition and vertex_set ...\n")
	writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
	for gene_no, mcl_id_go_no_list in gene_no2mcl_id_go_no_list.iteritems():
		cmp_condition_vertex_result = get_mcl_id_sharing(mcl_id_go_no_list, mcl_id2recurrence_array_vertex_set)
		for row in cmp_condition_vertex_result:
			writer.writerow([gene_no]+row)
	sys.stderr.write("Done.\n")

