#!/usr/bin/env python
"""
Usage: mcl_result_stat2file.py -k SCHEMA [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	mcl_result(default) or fim_result
	-g ..., --organism=...	specify the organism whose transcription
		relationship will be used, hs(default). Two-letter abbreviation.
	-b, --bonferroni	bonferroni correction
	-w, --wu	apply Wu's strategy(default, ignore it)
	-c, --commit	commit the database transaction(ignore it)
	-r, --report	report the progress(a number)
	-h, --help              show this help
	
Examples:
	mcl_result_stat2file.py -k shu
	mcl_result_stat2file.py -k shu_whole -t fim_result

Description:
	Program for computing cluster p-value vectors.
"""

import sys,os,psycopg,pickle,getopt
from rpy import r

class mcl_result_stat2file:
	def __init__(self, hostname, dbname, schema, table, organism, bonferroni=0, report=0, wu=1, needcommit=0):
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.table = table
		self.organism = organism
		self.bonferroni = int(bonferroni)
		self.report = int(report)
		self.wu = int(wu)
		self.needcommit = int(needcommit)
		self.global_go_no2go_name = {}
		self.global_go_no_to_size_dict = {}
		self.global_gene_to_go_dict = {}
		self.gene_no2gene_id = {}
		self.gene_no2tf = {}
		self.no_of_records = 0
		self.logfile = open('/tmp/mcl_result_stat2file.log','w')
		self.org_short2long = {'at':'Arabidopsis thaliana',
							'ce':'Caenorhabditis elegans',
							'dm':'Drosophila melanogaster',
							'hs':'Homo sapiens',
							'mm':'Mus musculus',
							'sc':'Saccharomyces cerevisiae',
							'rn':'Rattus norvegicus'}
		
	def dstruc_loadin(self):
		sys.stderr.write("Loading Data STructure...")
		#setup global_go_no2go_name and global_go_no_to_size_dict
		if self.wu:
			self.curs.execute("select name,go_no,array_upper(gene_array,1) from go")
		else:
			self.curs.execute("select name,go_no,array_upper(gene_array,1) from go where go_no!=0")
		rows = self.curs.fetchall()
		for row in rows:
			self.global_go_no2go_name[row[1]] = row[0]
			self.global_go_no_to_size_dict[row[1]] = row[2]
		self.no_of_functions = len(self.global_go_no_to_size_dict)
		#setup global_gene_to_go_dict and gene_no2gene_id
		if self.wu:
			self.curs.execute("select gene_no,go_functions,gene_id from gene")
		else:
			self.curs.execute("select gene_no,go_functions,gene_id from gene where known=TRUE")
		rows = self.curs.fetchall()
		for row in rows:
			self.gene_no2gene_id[row[0]] = row[2]
			self.global_gene_to_go_dict[row[0]] = []
			go_functions_list = row[1][1:-1].split(',')
			for go_no in go_functions_list:
				self.global_gene_to_go_dict[row[0]].append(int(go_no))
		self.no_of_genes = len(self.global_gene_to_go_dict)
		#setup gene_no2tf
		self.curs.execute("select g.gene_no, gt.tf from gene_backup g, gene2tf gt where \
			g.gene_id=gt.gene_id and gt.organism='%s'"%self.org_short2long[self.organism])
		rows = self.curs.fetchall()
		for row in rows:
			if self.gene_no2tf.has_key(row[0]):
				self.gene_no2tf[row[0]].append(row[1])
			else:
				self.gene_no2tf[row[0]] = [row[1]]
		
		sys.stderr.write("Done\n")
		
	def run(self):
		self.curs.execute("begin")
		self.curs.execute("DECLARE crs CURSOR FOR select mcl_id,vertex_set from %s"%self.table)
		self.curs.execute("fetch 5000 from crs")
		rows = self.curs.fetchall()
		while rows:
			for row in rows:
				mcl_id = row[0]
				vertex_set = row[1]
				self._cluster_stat(mcl_id, vertex_set)
				
			if self.report:
				sys.stderr.write("%s%s"%("\x08"*20,self.no_of_records))
			self.curs.execute("fetch 5000 from crs")
			rows = self.curs.fetchall()
		sys.stderr.write('\n\tTotal %d records.\n'%self.no_of_records)
		
	def _cluster_stat(self, mcl_id, vertex_set):
		vertex_list = vertex_set[1:-1].split(',')
		vertex_list = map(int, vertex_list)
		vertex_list_gene_symbol = []
		for vertex in vertex_list:
			vertex_list_gene_symbol.append(self.gene_no2gene_id[vertex])
		cluster_size = len(vertex_list)
		local_go_no_dict = self.local_go_no_dict_construct(vertex_list)
		
		if local_go_no_dict == {}:
			self.logfile.write('%d %s: local_go_no_dict empty\n'%(mcl_id, repr(vertex_set)))
			return
		for go_no in local_go_no_dict:
			if self.wu:
			# code after 'or' deals with the situation that Jasmine's strategy is applied to whole gene-set(unknown included)
				x = len(local_go_no_dict[go_no])
				m = self.global_go_no_to_size_dict[go_no]
				n = self.no_of_genes - m
				k = cluster_size
			else:
				pass
			if self.bonferroni:
				p_value = r.phyper(x-1,m,n,k,lower_tail = r.FALSE)*len(local_go_no_dict)
			else:
				p_value = r.phyper(x-1,m,n,k,lower_tail = r.FALSE)
			if x >1:
				#this function must have more than one gene associated.
				transfac_dict = self.transfac_dict_construct(local_go_no_dict[go_no])
				gene_id_list = []
				for gene_no in local_go_no_dict[go_no]:
					gene_id_list.append(self.gene_no2gene_id[gene_no])
				self.logfile.write('%d\t%s\t%d\t%d\t%d\t%d\t%f\t%s\t%s\t%s\n'%(mcl_id,\
					self.global_go_no2go_name[go_no], x, k, m, self.no_of_genes,\
					p_value, '|'.join(gene_id_list), repr(vertex_list_gene_symbol), repr(transfac_dict)))

		self.no_of_records += 1
		
	def local_go_no_dict_construct(self, vertex_list):
		'''
		construct a local go_no:size dictionary for a specific cluster.
		'''
		local_go_no_dict = {}
		for gene_no in vertex_list:
			go_no_list = self.global_gene_to_go_dict[gene_no]
			for go_no in go_no_list:
				if self.global_go_no_to_size_dict[go_no] > 1:
				#population singleton is discarded
					if go_no in local_go_no_dict:
						local_go_no_dict[go_no].append(gene_no)
					else:
						local_go_no_dict[go_no] = [gene_no]
		return local_go_no_dict

	def transfac_dict_construct(self, vertex_list):
		transfac_dict = {}
		subcluster_size = len(vertex_list)
		for gene_no in vertex_list:
			if self.gene_no2tf.has_key(gene_no):
				transfac_list = self.gene_no2tf[gene_no]
				for transfac in transfac_list:
					if transfac_dict.has_key(transfac):
						transfac_dict[transfac] += 1
					else:
						transfac_dict[transfac] = 1
		for transfac in transfac_dict.keys():
			#we will modify the dictionary during the iteration, so not iterate over transfac_dict
			#the transfactors that have < subcluster_size-1 genes associated are deleted.
			#if subcluster_size is 2, threshhold becomes subcluster_size.
			if subcluster_size == 2:
				if transfac_dict[transfac] < subcluster_size:
					del transfac_dict[transfac]
			elif subcluster_size > 2:
				if transfac_dict[transfac] < subcluster_size-1:
					del transfac_dict[transfac]
				
		return transfac_dict

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hrz:d:k:t:g:bcw", ["help", "report", "hostname=", "dbname=",\
			"table=", "schema=", "organism=", "bonferroni", "commit", "wu"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = 'mcl_result'
	organism = 'hs'
	bonferroni = 0
	commit = 0
	report = 0
	wu = 1
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
		elif opt in ("-g", "--organism"):
			organism = arg
		elif opt in ("-b", "--bonferroni"):
			bonferroni = 1
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-w", "--wu"):
			wu = 1

	if schema:
		instance = mcl_result_stat2file(hostname, dbname, schema, table, organism, bonferroni, report, wu, commit)
		instance.dstruc_loadin()
		instance.run()

	else:
		print __doc__
		sys.exit(2)
