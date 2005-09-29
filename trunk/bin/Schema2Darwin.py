#!/usr/bin/env python
"""
Usage: Schema2Darwin.py -k SCHEMA -f OFNAME -a ACC_CUTOFF -o OUTPUT_DIR [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-f ...,	the name of the initial cluster file(algorithm output)
	-a ...,	accuracy cutoff used to prediction, 0.6(default)
	-o ...,	output_dir
	-g ...,	organism, (mm, default)
	-b,	debug version.
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	Schema2Darwin.py -k mm_fim_97 -f mm_fim_97m4x200 -a 0.6 -o /tmp/yuhuang/
	
Description:
	Program to convert results of one schema into darwin format.
	Including TF, cluster, prediction from p_gene_table, gene_p_table, good_cluster_table
	and cluster_bs_table.

"""

import sys, os, csv, getopt
sys.path += [os.path.expanduser('~/script/annot/bin')]
from codense.common import db_connect, get_gene_no2gene_id, get_mt_no2tf_name, \
	get_mcl_id2tf_set, dict_map, get_go_no2name, org_short2long, org2tax_id, get_gene_id2gene_symbol, dict_transfer
from threading import *

class tf_darwin_format(Thread):
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, ofname=None, acc_cut_off=None, \
		output_fname=None, gene_no2id=None, go_no2id=None, mt_no2tf_name=None, debug=0, report=0):
		Thread.__init__(self)
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.ofname = ofname
		self.acc_cut_off = float(acc_cut_off)
		self.output_fname = output_fname
		self.gene_no2id = gene_no2id
		self.go_no2id = go_no2id
		self.mt_no2tf_name = mt_no2tf_name
		self.debug = int(debug)
		self.report = int(report)
	
	def _tf_darwin_format(self, curs, good_cluster_table, output_fname, gene_no2id, mcl_id2tf_set):
		"""
		format:
			r:=[
			[mcl_id, [gene1, gene2, ...], [ [TF1], [hyper_p_value] ], [ [TF2, TF3], [hyper_p_value] ], ... ],
			[...],
			[]]:
		"""
		sys.stderr.write("TF...\n")
		of = open(output_fname, 'w')
		of.write('r:=[\n')
		curs.execute("DECLARE crs CURSOR FOR select mcl_id, vertex_set from %s"%good_cluster_table)
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		while rows:
			for row in rows:
				mcl_id, vertex_set = row
				if mcl_id in mcl_id2tf_set:
					vertex_set = vertex_set[1:-1].split(',')
					vertex_set = map(int, vertex_set)
					vertex_set = dict_map(gene_no2id, vertex_set)
					tf_list = list(mcl_id2tf_set[mcl_id])
					tf_list = map(list, tf_list)	#first transform to list, so will have []
					for i in range(len(tf_list)):
						tf_list[i] = map(list, tf_list[i])	#one tf_list[i] is (tf_name_tuple, ratio_tuple)
					tf_list = map(repr, tf_list)	#second transform inner list to string
					row = [repr(mcl_id), repr(vertex_set)] + tf_list
					of.write('[%s],\n'%(','.join(row)))
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		of.write('[]]:\n')	#add the last blank list
		del of
		sys.stderr.write("TF darwin format done.\n")
		
	def run(self):
		if self.ofname and self.acc_cut_off:
			splat_table = 'splat_%s'%self.ofname
			mcl_table = 'mcl_%s'%self.ofname
			
			p_gene_table = 'p_gene_%s_e5'%self.ofname
			acc_int=int(self.acc_cut_off*100)
			gene_p_table='gene_p_%s_e5_a%s'%(self.ofname, acc_int)
			good_cluster_table = 'good_cl_%s_e5_a%s'%(self.ofname, acc_int)
			cluster_bs_table = 'cluster_bs_%s_e5_a%s'%(self.ofname, acc_int)
			
		else:
			sys.stderr.write("ofname: %s and acc_cut_off: %s, NOT VALID\n"%(self.ofname, self.acc_cut_off))
			sys.exit(2)
		conn, curs = db_connect(self.hostname, self.dbname, self.schema)
		mcl_id2tf_set = get_mcl_id2tf_set(curs, cluster_bs_table, self.mt_no2tf_name)
		self._tf_darwin_format(curs, good_cluster_table, self.output_fname, self.gene_no2id, mcl_id2tf_set)
		del conn, curs

class cluster_darwin_format(Thread):
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, ofname=None, acc_cut_off=None, \
		output_fname=None, gene_no2id=None, go_no2id=None, mt_no2tf_name=None, debug=0, report=0):
		Thread.__init__(self)
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.ofname = ofname
		self.acc_cut_off = float(acc_cut_off)
		self.output_fname = output_fname
		self.gene_no2id = gene_no2id
		self.go_no2id = go_no2id
		self.mt_no2tf_name = mt_no2tf_name
		self.debug = int(debug)
		self.report = int(report)

	
	def _cluster_darwin_format(self, curs, good_cluster_table, gene_no2id, go_no2id, output_fname):
		"""
		format:
			r:=[
			[mcl_id, vertex_set, recurrence_array, recurrence, connectivity, unknown_ratio, size, go_id_list, p_value_list]
			[...],
			...
			[]]:
			
		"""
		sys.stderr.write("cluster...\n")
		of = open(output_fname, 'w')
		of.write('r:=[\n')
		curs.execute("DECLARE crs CURSOR FOR select mcl_id, vertex_set, recurrence_array, recurrence, \
			connectivity, unknown_ratio, size, go_no_list, p_value_list from %s"%good_cluster_table)
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		while rows:
			for row in rows:
				mcl_id, vertex_set, recurrence_array, recurrence, connectivity, unknown_ratio,\
					size, go_no_list, p_value_list = row
				vertex_set = vertex_set[1:-1].split(',')
				vertex_set = map(int, vertex_set)
				vertex_set = dict_map(gene_no2id, vertex_set)
				recurrence_array = '[' + recurrence_array[1:-1] + ']'
				go_no_list = go_no_list[1:-1].split(',')
				go_no_list = map(int, go_no_list)
				go_id_list = dict_map(go_no2id, go_no_list)
				p_value_list = '[' + p_value_list[1:-1] + ']'
				of.write('[%s, %s, %s, %s, %s, %s, %s, %s, %s],\n'%(mcl_id, repr(vertex_set), recurrence_array,\
					recurrence, connectivity, unknown_ratio, size, repr(go_id_list), p_value_list))
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		of.write('[]]:\n')	#add the last blank list
		del of
		sys.stderr.write("cluster darwin format done.\n")
		
	def run(self):
		if self.ofname and self.acc_cut_off:
			splat_table = 'splat_%s'%self.ofname
			mcl_table = 'mcl_%s'%self.ofname
			
			p_gene_table = 'p_gene_%s_e5'%self.ofname
			acc_int=int(self.acc_cut_off*100)
			gene_p_table='gene_p_%s_e5_a%s'%(self.ofname, acc_int)
			good_cluster_table = 'good_cl_%s_e5_a%s'%(self.ofname, acc_int)
			cluster_bs_table = 'cluster_bs_%s_e5_a%s'%(self.ofname, acc_int)
		else:
			sys.stderr.write("ofname: %s and acc_cut_off: %s, NOT VALID\n"%(self.ofname, self.acc_cut_off))
			sys.exit(2)
		conn, curs = db_connect(self.hostname, self.dbname, self.schema)
		self._cluster_darwin_format(curs, good_cluster_table, self.gene_no2id, self.go_no2id, self.output_fname)
		del conn, curs


class prediction_darwin_format(Thread):
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, ofname=None, acc_cut_off=None, \
		output_fname=None, gene_no2id=None, go_no2id=None, mt_no2tf_name=None, debug=0, report=0):
		Thread.__init__(self)
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.ofname = ofname
		self.acc_cut_off = float(acc_cut_off)
		self.output_fname = output_fname
		self.gene_no2id = gene_no2id
		self.go_no2id = go_no2id
		self.mt_no2tf_name = mt_no2tf_name
		self.debug = int(debug)
		self.report = int(report)
	
	def _prediction_darwin_format(self, curs, p_gene_table, gene_p_table, new_p_gene_table, gene_no2id, go_no2id, output_fname):
		"""
		format:
			r:=[
			[gene_id, go_id, is_correct_lca, p_value, mcl_id, lca_list],
			[...],
			[]]:
		"""
		sys.stderr.write("prediction...\n")
		of = open(output_fname, 'w')
		of.write('r:=[\n')
		curs.execute("select p.* into %s from %s p, %s  g where g.p_gene_id = p.p_gene_id"%\
			(new_p_gene_table, p_gene_table, gene_p_table))
		curs.execute("DECLARE crs CURSOR FOR select gene_no, go_no, is_correct_lca, avg_p_value, mcl_id, lca_list from %s"%new_p_gene_table)
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		while rows:
			for row in rows:
				gene_no, go_no, is_correct_lca, p_value, mcl_id, lca_list = row
				lca_list = lca_list[1:-1].split(',')
				lca_list = map(int, lca_list)
				lca_list = dict_map(go_no2id, lca_list)
				of.write('[%s, %s, %s, %s, %s, %s],\n'%(gene_no2id[gene_no], go_no2id[go_no], is_correct_lca,\
					p_value, mcl_id, repr(lca_list)))
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		of.write('[]]:\n')	#add the last blank list
		del of
		sys.stderr.write("prediction darwin format done.\n")
	
	def run(self):
		if self.ofname and self.acc_cut_off:
			splat_table = 'splat_%s'%self.ofname
			mcl_table = 'mcl_%s'%self.ofname
			
			p_gene_table = 'p_gene_%s_e5'%self.ofname
			acc_int=int(self.acc_cut_off*100)
			new_p_gene_table = 'p_gene_%s_e5_a%s'%(self.ofname, acc_int)	#different from the others
			gene_p_table='gene_p_%s_e5_a%s'%(self.ofname, acc_int)
			good_cluster_table = 'good_cl_%s_e5_a%s'%(self.ofname, acc_int)
			cluster_bs_table = 'cluster_bs_%s_e5_a%s'%(self.ofname, acc_int)
			
		else:
			sys.stderr.write("ofname: %s and acc_cut_off: %s, NOT VALID\n"%(self.ofname, self.acc_cut_off))
			sys.exit(2)
		conn, curs = db_connect(self.hostname, self.dbname, self.schema)
		self._prediction_darwin_format(curs, p_gene_table, gene_p_table, new_p_gene_table, self.gene_no2id, self.go_no2id, self.output_fname)
		conn.commit()	#for the new_p_gene_table
		del conn, curs
	
class Schema2Darwin:
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, ofname=None, acc_cut_off=None, \
		output_dir=None, organism='mm', debug=0, report=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.ofname = ofname
		self.acc_cut_off = float(acc_cut_off)
		self.output_dir = output_dir
		self.organism = org_short2long(organism)
		self.debug = int(debug)
		self.report = int(report)
	
	def run(self):
		"""
		09-28-05
		"""
		if self.ofname and self.acc_cut_off:
			splat_table = 'splat_%s'%self.ofname
			mcl_table = 'mcl_%s'%self.ofname
			
			p_gene_table = 'p_gene_%s_e5'%self.ofname
			acc_int=int(self.acc_cut_off*100)
			gene_p_table='gene_p_%s_e5_a%s'%(self.ofname, acc_int)
			good_cluster_table = 'good_cl_%s_e5_a%s'%(self.ofname, acc_int)
			cluster_bs_table = 'cluster_bs_%s_e5_a%s'%(self.ofname, acc_int)
			
			tf_darwin_ofname = os.path.join(self.output_dir, '%s_e5_a%s.tf.darwin'%(self.ofname, acc_int))
			cluster_darwin_ofname = os.path.join(self.output_dir, '%s_e5_a%s.cluster.darwin'%(self.ofname, acc_int))
			prediction_darwin_ofname = os.path.join(self.output_dir, '%s_e5_a%s.prediction.darwin'%(self.ofname, acc_int))
		else:
			sys.stderr.write("ofname: %s and acc_cut_off: %s, NOT VALID\n"%(self.ofname, self.acc_cut_off))
			sys.exit(2)
		
		if not os.path.isdir(self.output_dir):
			os.makedirs(self.output_dir)
		conn, curs = db_connect(self.hostname, self.dbname, self.schema)
		
		tax_id = org2tax_id(self.organism)
		#gene_no2id = get_gene_no2gene_id(curs)	#Watch, if unigene, should use this.
		gene_id2symbol = get_gene_id2gene_symbol(curs, tax_id)
		#gene_no2symbol = dict_transfer(gene_no2id, gene_id2symbol)
			#Jasmine wants the gene symbol 09-28-05
			#gene_id is integer in gene.gene table and same as gene_no, so just use it.
		go_no2name = get_go_no2name(curs)	#09-28-05 Jasmine wants the go_name, not go_id
		mt_no2tf_name = get_mt_no2tf_name()
		
		instance1 =tf_darwin_format(self.hostname, self.dbname, self.schema, self.ofname, self.acc_cut_off, \
			tf_darwin_ofname, gene_id2symbol, go_no2name, mt_no2tf_name, debug, report)
		instance2 = cluster_darwin_format(self.hostname, self.dbname, self.schema, self.ofname, self.acc_cut_off, \
			cluster_darwin_ofname, gene_id2symbol, go_no2name, mt_no2tf_name, debug, report)
		#instance3 = prediction_darwin_format(self.hostname, self.dbname, self.schema, self.ofname, self.acc_cut_off, \
			#prediction_darwin_ofname, gene_id2symbol, go_no2name, mt_no2tf_name, debug, report)
		instance1.start()
		instance2.start()
		#instance3.start()
		instance1.join()
		instance2.join()
		#instance3.join()



if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:f:a:o:g:br", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	ofname = None
	acc_cut_off = 0.6
	output_dir = None
	organism = 'mm'
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
		elif opt in ("-f"):
			ofname = arg
		elif opt in ("-a"):
			acc_cut_off = float(arg)
		elif opt in ("-o"):
			output_dir = arg
		elif opt in ("-g"):
			organism = arg
		elif opt in ("-b"):
			debug = 1
		elif opt in ("-r"):
			report = 1
	if schema and ofname and acc_cut_off and output_dir:
		instance = Schema2Darwin(hostname, dbname, schema, ofname, acc_cut_off, output_dir,\
			organism, debug, report)
		instance.run()
		
	else:
		print __doc__
		sys.exit(2)
