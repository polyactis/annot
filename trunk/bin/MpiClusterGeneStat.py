#!/usr/bin/env mpipython
"""
Usage: MpiClusterGeneStat.py -k SCHEMA -s MCLTABLE -o OUTPUTFILE [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-s ..., --source_table=...	which source table, mcl_result(default) or fim_result
	-p ..., --output=...	specifiy the filename to output the cluster stat results
	-g ..., --gene_table=...	gene_table for gene_stat.py
	-t ..., --times_nodes=...	how many times of the number of original nodes(which is size-1),
		2.0(default), this is used to balance the amount of work per node
	-c, --commit	commit this database transaction
	-b, --debug	debug version.
	-r, --report	enable report flag
	-h, --help	Display the usage infomation.
	
Examples:
	mpirun.mpich -np 30 -machinefile ~/hostfile /usr/bin/mpipython ~/script/annot/bin/MpiClusterGeneStat.py
	-k mm_fim_97 -s mcl_test5 -p /scratch/00/yuhuang/cluster_stat/cluster_test5 -g p_gene_test5_e5
	
Description:
	This program parallels cluster_stat.py, other default parameters are
	uniformity = 0, min_node_size = 0, min_node_depth = 5, bonferroni = 0
	commit = 0, report = 0, log = 0, wu = 1.
	and also gene_stat.py, depth_cut_off = 5, leave_one_out = 1, wu = 1,
	report = 0, subgraph_cut_off = 0, debug = 0, new_table = 0
"""

import sys, os, getopt, csv, math, Numeric
sys.path += [os.path.expanduser('~/script/annot/bin')]
from Scientific import MPI
from codense.common import system_call, mpi_schedule_jobs, mpi_synchronize, db_connect
from netmine_wrapper import netmine_wrapper
from codense.codense2db import codense2db
from sets import Set
from gene_stat import gene_stat

def node_cluster_stat(index, parameter_list):
	"""
	08-14-05
	"""
	index  = int(index)	#08-14-05	convert it to integer form, cause it's transferred as a string
	from cluster_stat import cluster_stat
	hostname, dbname, schema, source_table, output, gene_table, \
	gene_table_commit, OffsetLimitList, debug = parameter_list
	
	target_table = 'cluster_stat'
	offset, limit = OffsetLimitList[index]
	output = '%s.%s'%(output, index)
	uniformity = 0
	min_node_size = 0
	min_node_depth = 5
	bonferroni = 0
	commit = 0
	report = 0
	log = 0
	wu = 1
	if debug:
		sys.stderr.write("Index: %s; Offset: %s; Limit: %s\n"%(index, offset, limit))
	instance = cluster_stat(hostname, dbname, schema, source_table, target_table, \
		offset, limit, output, bonferroni, report, log, wu, commit, uniformity, \
		min_node_size, min_node_depth)
	instance.run()
	del instance
	
	depth_cut_off = 5
	dir_files = output
	leave_one_out = 1
	wu = 1
	report = 0
	commit = gene_table_commit
	subgraph_cut_off = 0
	debug = 0
	new_table = 0
	recurrence_gap_size = 2
	connectivity_gap_size = 2
	instance = gene_stat(hostname, dbname, schema, target_table, source_table, \
		leave_one_out, wu, report, depth_cut_off, dir_files, commit, gene_table, \
		subgraph_cut_off, debug, new_table, recurrence_gap_size, connectivity_gap_size)
	instance.run()
	del instance
	return output

class MpiClusterGeneStat:
	"""
	08-14-05
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, source_table=None,\
		output=None, gene_table=None, times_nodes=2.0, commit=0, report=0, debug=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.source_table = source_table
		self.output = output
		self.gene_table = gene_table
		self.times_nodes = float(times_nodes)
		self.commit = int(commit)
		self.report = int(report)
		self.debug = int(debug)
	
	def createOffsetLimitList(self, curs, source_table, no_of_nodes):
		"""
		08-14-05
			create [offset,limit] list
			no_of_nodes is the number of jobs to schedule, not the total number of nodes
		"""
		sys.stderr.write("Creating OffsetLimitList...")
		curs.execute("select count(mcl_id) from %s"%source_table)
		rows = curs.fetchall()
		no_of_clusters = rows[0][0]
		OffsetLimitList = []
		step, remainder = divmod(no_of_clusters, no_of_nodes-1)
		for i in range(no_of_nodes):
			OffsetLimitList.append([i*step, step])	#the last step might be too long, but it's ok.
		sys.stderr.write("Done.\n")
		return OffsetLimitList
	
	def run(self):
		"""
		08-14-05
		"""
		communicator = MPI.world.duplicate()
		fake_no_of_nodes = int((communicator.size-1)*times_nodes)	#NOTICE: fake_no_of_nodes is used to enlarge(or shrink) the actual number of nodes,
			#to balance the amount of work on each node
		OffsetLimitList = Numeric.zeros((fake_no_of_nodes,2), Numeric.Int)
		if communicator.rank == 0:
			(conn, curs) = db_connect(self.hostname, self.dbname, self.schema)
			OffsetLimitList = self.createOffsetLimitList(curs, self.source_table, fake_no_of_nodes)
			OffsetLimitList = Numeric.array(OffsetLimitList, Numeric.Int)	#transform it into Numeric array to broadcast()
			if self.commit:	#08-14-05	create the gene_table
				instance = gene_stat()
				instance.createGeneTable(curs, self.gene_table)
				curs.execute('end')
			if self.debug:
				sys.stderr.write("OffsetLimitList: %s"%repr(OffsetLimitList))
			del conn, curs
		
		communicator.broadcast(OffsetLimitList, 0)	#share the OffsetLimitList
		
		mpi_synchronize(communicator)
		job_list = range(len(OffsetLimitList))	#corresponding to the indices in the OffsetLimitList
		parameter_list =[self.hostname, self.dbname, self.schema, self.source_table, self.output, \
			self.gene_table, self.commit, OffsetLimitList, self.debug]
		if self.debug:
			sys.stderr.write("The common parameter_list is %s.\n"%repr(parameter_list))
		of_name_list = mpi_schedule_jobs(communicator, job_list, node_cluster_stat, parameter_list, self.debug)
		
		mpi_synchronize(communicator)
		
		#collecting 08-14-05 not really necessary, but just to make the number of files small
		if communicator.rank==0:
			netmine_wrapper_instance = netmine_wrapper()
			netmine_wrapper_instance.collect_and_merge_output(of_name_list, self.output)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:s:p:g:t:cbr", ["help", "hostname=", \
			"dbname=", "schema=", "source_table=", "output=", "gene_table=", "times_nodes=",\
			"commit", "debug", "report"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	source_table = None
	output = None
	gene_table = None
	times_nodes = 2.0
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
		elif opt in ("-s", "--source_table"):
			source_table = arg
		elif opt in ("-p", "--output"):
			output = arg
		elif opt in ("-g", "--gene_table"):
			gene_table = arg
		elif opt in ("-t", "--times_nodes"):
			times_nodes = float(arg)
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	if schema and source_table and output and gene_table:
		instance = MpiClusterGeneStat(hostname, dbname, schema, \
			source_table, output, gene_table, times_nodes, commit, report, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
