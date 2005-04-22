#!/usr/bin/env python
"""
Usage: cluster_stat.py -k SCHEMA [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-s ..., --source_table=...	which source table, mcl_result(default) or fim_result
	-t ..., --target_table=...	which target table, cluster_stat(default)
	-o ..., --offset=...	the number of rows to skip before returning rows, 0 (default)
	-m ..., --limit=...	the maximum number of rows to return, all (default)
	-p ..., --output=...	specifiy the filename to output the cluster stat results
	-u ..., --uniformity=...	the percentage of associated-genes over total known genes, 0.5(default)
	-b, --bonferroni	bonferroni correction
	-w, --wu	apply Wu's strategy(Default is Jasmine's strategy)
	-c, --commit	commit the database transaction
	-r, --report	report the progress(a number)
	-l, --log	record down some stuff in the logfile(cluster_stat.log)
	-h, --help	show this help
	
Examples:
	cluster_stat.py -k shu -c
	cluster_stat.py -k ming1
		:jasmine's strategy
	cluster_stat.py -k sc_yh60_splat_5 -b -c
		:jasmine's strategy + bonferroni correction
	cluster_stat.py -k sc_yh_60_fp -s fim_result -c -w -b
		:wu's strategy + bonferroni correction
	cluster_stat.py -k sc_yh60_splat_5 -s mcl_result2 -t cluster_stat2 -c -b -r
	cluster_stat.py -k sc_38_no_informative -o 10*2 -m 10*1 -p /tmp/tmp -l
		-w -s mcl_result_merge_sup_6 -r -b

Description:
	Program for computing cluster p-value vectors(leave-one-out).
	The length of p_value_vector is subject to the number of functional categories.
	The difference between two strategies:
	validatoin:	Jasmine's only works on known-gene clusters.
		Wu's works on whole-gene clusters.
	prediciton:  Jasmine's only works on clusters with only one unknown gene.
		Wu's works on all clusters with any number of unknown genes and needs
		a cut_off for unknown functional class.
	The above difference is never true(08/23/04). Modified Jasmine's strategy is
	quite similar to Wu's strategy. If you call gene_stat.py after this version of
	cluster_stat.py, add the -w parameter to gene_stat.py
"""

import sys,os,psycopg,pickle,getopt
from rpy import r
from sets import Set

class cluster_stat:
	def __init__(self, hostname, dbname, schema, source_table, target_table, offset, limit, \
		output, bonferroni=0, report=0, log=0, wu=0, needcommit=0, uniformity=0.5):
		"""
		04-18-05
			add parameter uniformity
		"""
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.source_table = source_table
		self.target_table = target_table
		self.offset = offset
		self.limit = limit
		self.output = output
		if self.output:
			#filename exists
			self.outf = open(self.output, 'w')

		self.bonferroni = int(bonferroni)
		self.report = int(report)
		self.log = int(log)
		self.wu = int(wu)
		self.needcommit = int(needcommit)
		self.uniformity = float(uniformity)
		self.global_go_id_to_no_dict = {}
		self.global_go_no_to_size_dict = {}
		self.global_gene_to_go_dict = {}
		self.no_of_records = 0
		if self.log:
			self.logfile = open('/tmp/cluster_stat.log','w')
		self.cluster_memory = {}
		#mapping between gene_no and go_no list
		self.known_genes_dict = {}

	def dstruc_loadin(self):
		self.curs.execute("select go_id,go_no,array_upper(gene_array,1) from go")
		rows = self.curs.fetchall()
		for row in rows:
			self.global_go_id_to_no_dict[row[0]] = row[1]
			self.global_go_no_to_size_dict[row[1]] = row[2]
		self.no_of_functions = len(self.global_go_no_to_size_dict)
		
		self.curs.execute("select gene_no,go_functions from gene")	
		rows = self.curs.fetchall()
		for row in rows:
			self.global_gene_to_go_dict[row[0]] = []
			go_functions_list = row[1][1:-1].split(',')
			for go_no in go_functions_list:
				self.global_gene_to_go_dict[row[0]].append(int(go_no))
		self.no_of_genes = len(self.global_gene_to_go_dict)
		
		#setup self.known_genes_dict
		self.curs.execute("select gene_no,go_functions from gene where known=TRUE")
		rows = self.curs.fetchall()
		for row in rows:
			go_functions_list = row[1][1:-1].split(',')
			self.known_genes_dict[row[0]] = Set()
			for go_no in go_functions_list:
				self.known_genes_dict[row[0]].add(int(go_no))
		
	def run(self):
		if self.output and self.needcommit:
			sys.stderr.write("output and needcommit are two incompatible options.\n")
			sys.exit(2)
	
		if self.target_table != 'cluster_stat' and self.output == None:
			try:
				#03-18-05 cluster_stat_id is changed to be an integer, because the the name of xxx_cluster_stat_id_seq has limited length and got collisions.
				self.curs.execute("create table %s(\
					cluster_stat_id	integer,\
					mcl_id	integer,\
					leave_one_out	integer,\
					p_value_vector	float[],\
					connectivity	float)"%self.target_table)
			except:
				sys.stderr.write("Error occurred when creating table %s\n"%self.target_table)
		self.curs.execute("begin")
		self.curs.execute("DECLARE crs CURSOR FOR select mcl_id,vertex_set,connectivity from %s order by mcl_id offset %s limit %s"%\
			(self.source_table, self.offset, self.limit))
		self.curs.execute("fetch 5000 from crs")
		rows = self.curs.fetchall()
		while rows:
			for row in rows:
				mcl_id = row[0]
				vertex_set = row[1]
				connectivity = row[2]
				self._cluster_stat(mcl_id, vertex_set, connectivity)
				
			if self.report:
				sys.stderr.write("%s%s"%("\x08"*20,self.no_of_records))
			self.curs.execute("fetch 5000 from crs")
			rows = self.curs.fetchall()
		if self.needcommit:
			self.curs.execute("create index %s_connectivity_idx on %s(connectivity)"%(self.target_table, self.target_table))
			self.curs.execute("create index %s_mcl_id_idx on %s(mcl_id)"%(self.target_table, self.target_table))
			self.curs.execute("end")
			sys.stderr.write('\n\tTotal %d records.\n'%self.no_of_records)
		else:
			self.conn.rollback()
			sys.stderr.write('\n\tNo real updates\n')
		
	def _cluster_stat(self, mcl_id, vertex_set, connectivity):
		"""
		04-18-05
			add two important criteria to avoid the situation that hypergeometric test
			is powerless (the population size of the go-no is too small).
			1. percentage of associated-genes over total known genes >= uniformity (0.5 default)
			2. apart from the percentage, the absolute number is also needed in case the cluster is too small.
		04-19-05
			fix a bug. self._no_of_known_genes_of_the_cluster could be 0.
		"""	
		vertex_list_all = vertex_set[1:-1].split(',')
		vertex_list = []
		for i in range(len(vertex_list_all)):
			vertex_list_all[i] = int(vertex_list_all[i])
			if vertex_list_all[i] in self.global_gene_to_go_dict:
			#this filter will only be useful when Jasmine's strategy is applied to whole gene-set(unknown included)
				vertex_list.append(vertex_list_all[i])
		cluster_size = len(vertex_list)
		p_value_vector = [1] * self.no_of_functions
		self.local_go_no_dict_construct(vertex_list)
		for gene_no in vertex_list_all:
			self.go_no_dict_adjust(gene_no)
			for go_no in self._local_go_no_dict:
				if self.wu or (gene_no not in self.global_gene_to_go_dict):
				# code after 'or' deals with the situation that Jasmine's strategy is applied to whole gene-set(unknown included)
					x = self._local_go_no_dict[go_no]
					m = self._global_go_no_dict[go_no]
					n = self.no_of_genes - m
					k = cluster_size
				else:
					x = self._local_go_no_dict[go_no]
					m = self._global_go_no_dict[go_no]
					n = self.no_of_genes -1 - m
					k = cluster_size-1
				if self._no_of_known_genes_of_the_cluster == 0:
					go_no_ratio = 0
				else:
					go_no_ratio = float(x)/self._no_of_known_genes_of_the_cluster	#NOTE: it's different from no_of_known_genes_of_the_cluster
				if  go_no_ratio < self.uniformity and go_no!=0:	#It doesn't apply to the 0(unknown) category.
					#ignore the function category if its percentage is < uniformity
					continue
				if x < 3 and go_no!=0:	#apart from the percentage, the absolute number is also needed in case the cluster is too small.
					continue
				if self.bonferroni:
					p_value = r.phyper(x-1,m,n,k,lower_tail = r.FALSE)*len(self._local_go_no_dict)
				else:
					p_value = r.phyper(x-1,m,n,k,lower_tail = r.FALSE)
				if self.log:
					self.logfile.write('%d %d %d %d %d %d %d %f\n'%\
						(mcl_id,gene_no,go_no,x,m,n,k,p_value))
				p_value_vector[go_no] = p_value
			#for the unknown class, use the ratio instead of p_value, in accordance with mcl_result_stat.py
			if self.wu:
				p_value_vector[0] = self._local_go_no_dict[0]/float(cluster_size)
			else:
				#not wu's strategy, throw away the gene, the cluster_size is down by 1.
				if self._local_go_no_dict.has_key(0):
					#after leave_one_out, still unknown genes present
					p_value_vector[0] = self._local_go_no_dict[0]/float(cluster_size-1)
				else:
					#no unknown genes
					p_value_vector[0] = 1
			#03-18-05increment before inserted into table, cluster_stat_id starting from 1
			self.no_of_records += 1
			if self.output:
				self.outf.write('%d\t%d\t%s\t%f\n'%(mcl_id, gene_no, repr(p_value_vector), connectivity))
			elif self.needcommit:
				self.curs.execute("insert into %s(cluster_stat_id, mcl_id, leave_one_out, p_value_vector, connectivity)\
				values(%d, %d, %d, ARRAY%s, %8.6f)"%(self.target_table, self.no_of_records, mcl_id, gene_no, repr(p_value_vector), connectivity))

	def local_go_no_dict_construct(self, vertex_list):
		'''
		construct a local go_no:size dictionary for a specific cluster.
		
		04-18-05
			set the no_of_known_genes_of_the_cluster
		'''
		self.local_go_no_dict = {}
		self.no_of_known_genes_of_the_cluster = 0
		for gene_no in vertex_list:
			if gene_no in self.known_genes_dict:
				self.no_of_known_genes_of_the_cluster += 1
			go_no_list = self.global_gene_to_go_dict[gene_no]
			for go_no in go_no_list:
				if go_no in self.local_go_no_dict:
					self.local_go_no_dict[go_no] += 1
				else:
					self.local_go_no_dict[go_no] = 1

		
	def go_no_dict_adjust(self, gene_no):
		'''
		After one gene is left out, both local and global go_no:size dictionaries
		need to be adjusted.
		04-18-05:
			set _no_of_known_genes_of_the_cluster
		'''
		self._local_go_no_dict = self.local_go_no_dict.copy()
		self._global_go_no_dict = self.global_go_no_to_size_dict.copy()
		if gene_no not in self.global_gene_to_go_dict:
		#this filter will only be useful when Jasmine's strategy is applied to whole gene-set(unknown included)
			return
		if gene_no in self.known_genes_dict:
			self._no_of_known_genes_of_the_cluster = self.no_of_known_genes_of_the_cluster-1
		else:
			self._no_of_known_genes_of_the_cluster = self.no_of_known_genes_of_the_cluster
		
		go_no_list = self.global_gene_to_go_dict[gene_no]
		for go_no in go_no_list:
			self._local_go_no_dict[go_no] -= 1
			self._global_go_no_dict[go_no] -= 1
			if self._local_go_no_dict[go_no] == 0:
				del self._local_go_no_dict[go_no]
			if self._global_go_no_dict[go_no] == 0:
				del self._global_go_no_dict[go_no]
		if self.wu:
		# for Wu's strategy
			if self._local_go_no_dict.has_key(0):
				self._local_go_no_dict[0] += 1
			else:
				self._local_go_no_dict[0] = 1
			self._global_go_no_dict[0] += 1


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hrlz:d:k:s:t:o:m:p:u:bcw", \
			["help", "report", "log", "hostname=", "dbname=", "schema=", "source_table=", "target_table=", \
			"offset=", "limit=", "output=", "uniformity=", "bonferroni", "commit", "wu"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	source_table = 'mcl_result'
	target_table = 'cluster_stat'
	offset = 0
	limit = 'all'
	output = None
	uniformity = 0.5
	bonferroni = 0
	commit = 0
	report = 0
	log = 0
	wu = 0
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
		elif opt in ("-t", "--target_table"):
			target_table = arg
		elif opt in ("-o", "--offset"):
			offset = arg
		elif opt in ("-m", "--limit"):
			limit = arg
		elif opt in ("-p", "--output"):
			output = arg
		elif opt in ("-u", "--uniformity"):
			uniformity = float(arg)
		elif opt in ("-b", "--bonferroni"):
			bonferroni = 1
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-l", "--log"):
			log = 1
		elif opt in ("-w", "--wu"):
			wu = 1

	if schema:
		instance = cluster_stat(hostname, dbname, schema, source_table, target_table, \
			offset, limit, output, bonferroni, report, log, wu, commit, uniformity)
		instance.dstruc_loadin()
		instance.run()

	else:
		print __doc__
		sys.exit(2)
