#!/usr/bin/env python
"""
Usage: batch_stat.py -k SCHEMA [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	which table contains the cluster stat, cluster_stat(default)
	-m ..., --mcl_table=...	which table is the corresponding mcl_table of above table, mcl_result(default)
	-g ..., --tag=...		a tag to identify all the results, required if commit
	-p ..., --p_value_list=...,	a list of p_value's, "0.1,0.01,0.001,0.0001,0.00001"(default)
	-u ..., --unknown_list=...	a list of unknown_gene ratio, "0.25"(default)
	-n ..., --connectivity_list=...	a list of connectivity cut_off, "0.8"(default)
	-r ..., --recurrence_list=...	a list of recurrence cut_off, "5"(default)
	-s ..., --size_list=...	a list of cluster_size cut_off, "1000"(default)
	-e ..., --depth_list=...	the minimum depth for a go node to be valid, 3(default)
	-f ..., --dir_files=...	the directory containing all the files outputed by cluster_stat.py
	-j ..., --judger_type=...	how to judge predicted functions, 0(default), 1, 2
	-l, --leave_one_out	use the leave_one_out stat method(gene_stat).
			default is gene_stat_on_mcl_result.
	-w, --wu	Ignore it, default is also wu's strategy.
	-v, --dominant	Only assign the dominant function(s) to a gene.
	-c, --commit	commit the database transaction, records in table stat_plot_data
	-h, --help              show this help
	
Examples:
	batch_stat.py -k shu
	batch_stat.py -k shu -t sc_known -c
	batch_stat.py -k sc_yh60_splat_5 -l -t cluster_stat2 -m mcl_result2
		-g sc_yh60_splat_5_cluster_stat2 -s "4,6,8,10,14,18,22" -c
	batch_stat.py -k sc_yh60_splat_5 -m mcl_result2 -g sc_yh60_splat_5_mcl_result2
		-s "4,6,8,10,12" -c

Description:
	Parameter leave_one_out controls which stat class to call. For gene_stat,
	both the 'table' and 'mcl_table' must be given. For gene_stat_on_mcl_result,
	mcl_table is enough.
"""

import sys, os, psycopg, getopt
from gene_stat_plot import gene_stat
#from gene_stat_on_mcl_result import gene_stat_on_mcl_result

class batch_stat:
	def __init__(self, hostname, dbname, schema, table, mcl_table, tag, p_value_list, unknown_list, connectivity_list,\
			recurrence_list, size_list, leave_one_out, wu, judger_type, depth_list, dir_files=None, dominant=0, needcommit=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.table = table
		self.mcl_table = mcl_table
		self.tag = tag
		self.p_value_list = p_value_list
		self.unknown_list = unknown_list
		self.connectivity_list = connectivity_list
		self.recurrence_list = recurrence_list
		self.cluster_size_list = size_list
		self.leave_one_out = int(leave_one_out)
		self.wu = int(wu)
		self.judger_type = int(judger_type)
		self.depth_list = depth_list
		self.dir_files = dir_files
		self.dominant = int(dominant)
		self.needcommit = int(needcommit)
		#a dictionary mapping the user choices to stat classes
		#self.stat_class_dict = {1: gene_stat,
		#	0: gene_stat_on_mcl_result}
		self.stat_data = []


		#default parameters to invoke gene_stat
		self.report = 0
		self.log = 0
		self.needcommit_of_gene_stat = 0
		self.gene_table = 'p_gene'
		
	def run(self):
		instance = gene_stat(self.hostname, self.dbname, self.schema, self.table, self.mcl_table, self.p_value_list[0],\
			self.unknown_list[0], self.connectivity_list[0], self.recurrence_list[0], self.cluster_size_list[0], \
			self.leave_one_out, self.wu, self.report,\
			self.log, self.judger_type, self.depth_list[0], self.dir_files, self.needcommit_of_gene_stat,\
			self.gene_table, self.dominant)
		instance.dstruc_loadin()
		for unknown in self.unknown_list:
			for connectivity in self.connectivity_list:
				for recurrence in self.recurrence_list:
					for size in self.cluster_size_list:
						for p_value in self.p_value_list:
							for depth in self.depth_list:
								instance.parameter_reset_and_cleanup(unknown, connectivity, recurrence, size, p_value, depth)
								instance.run()
								if self.needcommit:
									self.curs.execute("insert into graph.stat_plot_data(tag, p_value_cut_off, unknown_cut_off,\
									connectivity_cut_off, recurrence_cut_off, cluster_size_cut_off, depth_cut_off, tp, tp_m, tp1, tp1_m,\
									tn, fp, fp_m, fn) values ('%s', %1.5f, %1.2f, %1.2f,  %d, %d, %d, %d,  %d, %d, %d, %d,  %d, %d, %d)"%\
									(self.tag, p_value, unknown, connectivity, recurrence, size, depth, instance.tp,   instance.tp_m,\
									instance.tp1, instance.tp1_m, instance.tn,  instance.fp, instance.fp_m, instance.fn))
									self.conn.commit()

		
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:m:g:p:u:n:r:s:e:j:f:lwvc", ["help", "hostname=", "dbname=", "schema=", "table=",\
			"mcl_table=", "tag=", "p_value_list=", "unknown_list=", "connectivity_list=", "recurrence_list=",\
			"size_list=", "depth_list=", "judger_type=", "dir_files=", "leave_one_out", "wu", "dominant", "commit"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = 'cluster_stat'
	mcl_table = 'mcl_result'
	tag = ''
	p_value_list = [0.001]
	unknown_list = [1.0]
	connectivity_list = [0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]
	recurrence_list = [5]
	size_list = [10000]
	depth_list = [4]
	dir_files = None
	judger_type = 0
	leave_one_out = 0
	wu = 1 
	dominant = 0
	commit = 0
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
		elif opt in ("-g", "--tag"):
			tag = arg
		elif opt in ("-p", "--p_value_list"):
			p_value_list = arg.split(',')
			p_value_list = map(float, p_value_list)
		elif opt in ("-u", "--unknown_list"):
			unknown_list = arg.split(',')
			unknown_list = map(float, unknown_list)
		elif opt in ("-n", "--connectivity_list"):
			connectivity_list = arg.split(',')
			connectivity_list = map(float, connectivity_list)
		elif opt in ("-r", "--recurrence_list"):
			recurrence_list = arg.split(',')
			recurrence_list = map(int, recurrence_list)
		elif opt in ("-s", "--size_list"):
			size_list = arg.split(',')
			size_list = map(int, size_list)
		elif opt in ("-e", "--depth_list"):
			depth_list = arg.split(',')
			depth_list = map(int, depth_list)
		elif opt in ("-f", "--dir_files"):
			dir_files = arg
		elif opt in ("-o", "--log"):
			log = 1
		elif opt in ("-j", "--judger_type"):
			judger_type = int(arg)
		elif opt in ("-l", "--leave_one_out"):
			leave_one_out = 1
		elif opt in ("-w", "--wu"):
			wu = 1
		elif opt in ("-v", "--dominant"):
			dominant = 1
		elif opt in ("-c", "--commit"):
			commit = 1
	
	if schema:
		if commit == 1 and tag == '':
			print __doc__
			sys.exit(2)
		instance = batch_stat(hostname, dbname, schema, table, mcl_table, tag, p_value_list, unknown_list,\
			connectivity_list, recurrence_list, size_list, leave_one_out, wu, judger_type, depth_list, dir_files, dominant, commit)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
