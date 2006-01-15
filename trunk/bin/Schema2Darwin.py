#!/usr/bin/env python
"""
Usage: Schema2Darwin.py -k SCHEMA -f OFNAME -a ACC_CUTOFF -o OUTPUT_DIR [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-f ...,	the name of the initial cluster file(algorithm output)
	-l ...,	lm_bit, 111(default)
	-a ...,	accuracy cutoff used to prediction, 0.6(default)
	-o ...,	output_dir
	-g ...,	organism
	-n,	running bit, (tf_darwin_format, cluster_darwin_format, prediction_darwin_format, 
		pattern_darwin_format), 1111(default)
	-b,	debug version.
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	ssh $HOSTNAME Schema2Darwin.py -k mm_fim_97 -f mm_fim_97m4x200 -a 0.6 -l 111 
		-o /tmp/yuhuang/ -g mm
	
Description:
	Program to convert results of one schema into darwin format.
	Including TF, cluster, prediction from p_gene_table, gene_p_table, good_cluster_table
	and cluster_bs_table.

	ssh $HOSTNAME is used for qsub system because it doesn't allow thread.

"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import sys, os, csv, getopt
from codense.common import db_connect, get_gene_no2gene_id, get_mt_no2tf_name, \
	get_mcl_id2tf_set, dict_map, get_go_no2name, org_short2long, org2tax_id, \
	get_gene_id2gene_symbol, dict_transfer, form_schema_tables
from threading import *
from sets import Set

class tf_darwin_format(Thread):
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, ofname=None, lm_bit='111', acc_cut_off=None, \
		output_fname=None, gene_no2id=None, go_no2id=None, mt_no2tf_name=None, debug=0, report=0):
		Thread.__init__(self)
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.ofname = ofname
		self.lm_bit = lm_bit
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
					vertex_set = dict_map(gene_no2id, vertex_set, type=2)
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
		if self.ofname and self.acc_cut_off and self.lm_bit:
			schema_instance = form_schema_tables(self.ofname, self.acc_cut_off, self.lm_bit)
			
		else:
			sys.stderr.write("ofname: %s and acc_cut_off: %s and lm_bit %s, NOT VALID\n"%(self.ofname, self.acc_cut_off, self.lm_bit))
			sys.exit(2)
		conn, curs = db_connect(self.hostname, self.dbname, self.schema)
		mcl_id2tf_set = get_mcl_id2tf_set(curs, schema_instance.cluster_bs_table, self.mt_no2tf_name)
		self._tf_darwin_format(curs, schema_instance.good_cluster_table, self.output_fname, self.gene_no2id, mcl_id2tf_set)
		del conn, curs

class cluster_darwin_format(Thread):
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, ofname=None, lm_bit='111', acc_cut_off=None, \
		output_fname=None, gene_no2id=None, go_no2id=None, mt_no2tf_name=None, debug=0, report=0):
		Thread.__init__(self)
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.ofname = ofname
		self.lm_bit = lm_bit
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
				vertex_set = dict_map(gene_no2id, vertex_set, type=2)
				recurrence_array = '[' + recurrence_array[1:-1] + ']'
				go_no_list = go_no_list[1:-1].split(',')
				go_no_list = map(int, go_no_list)
				go_id_list = dict_map(go_no2id, go_no_list, type=2)
				p_value_list = '[' + p_value_list[1:-1] + ']'
				of.write('[%s, %s, %s, %s, %s, %s, %s, %s, %s],\n'%(mcl_id, repr(vertex_set), recurrence_array,\
					recurrence, connectivity, unknown_ratio, size, repr(go_id_list), p_value_list))
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		of.write('[]]:\n')	#add the last blank list
		del of
		sys.stderr.write("cluster darwin format done.\n")
		
	def run(self):
		if self.ofname and self.acc_cut_off and self.lm_bit:
			schema_instance = form_schema_tables(self.ofname, self.acc_cut_off, self.lm_bit)
			
		else:
			sys.stderr.write("ofname: %s and acc_cut_off: %s and lm_bit %s, NOT VALID\n"%(self.ofname, self.acc_cut_off, self.lm_bit))
			sys.exit(2)
		conn, curs = db_connect(self.hostname, self.dbname, self.schema)
		self._cluster_darwin_format(curs, schema_instance.good_cluster_table, self.gene_no2id, self.go_no2id, self.output_fname)
		del conn, curs


class pattern_darwin_format(Thread):
	"""
	12-19-05
		output all patterns from pattern_table into darwin format
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, ofname=None, lm_bit='111', acc_cut_off=None, \
		output_fname=None, gene_no2id=None, go_no2id=None, mt_no2tf_name=None, debug=0, report=0):
		Thread.__init__(self)
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.ofname = ofname
		self.lm_bit = lm_bit
		self.acc_cut_off = float(acc_cut_off)
		self.output_fname = output_fname
		self.gene_no2id = gene_no2id
		self.go_no2id = go_no2id
		self.mt_no2tf_name = mt_no2tf_name
		self.debug = int(debug)
		self.report = int(report)

	def get_mcl_id_set_from_good_cluster_table(self, curs, good_cluster_table):
		"""
		01-14-06
			temporarily add, 
		"""
		sys.stderr.write("Getting mcl_id_set from good_cluster_table...\n")
		curs.execute("DECLARE crs CURSOR FOR select mcl_id from %s"%good_cluster_table)
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		mcl_id_set = Set()
		while rows:
			for row in rows:
				mcl_id_set.add(row[0])
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		curs.execute("close crs")
		sys.stderr.write("mcl_id_set done.\n")
		return mcl_id_set
	
	def _pattern_darwin_format(self, curs, pattern_table, gene_no2id, go_no2id, output_fname, mcl_id_set=None):
		"""
		format:
			r:=[
			[mcl_id, vertex_set, recurrence_array, recurrence, connectivity, unknown_ratio]
			[...],
			...
			[]]:
			
		"""
		sys.stderr.write("pattern...\n")
		of = open(output_fname, 'w')
		of.write('r:=[\n')
		curs.execute("DECLARE crs CURSOR FOR select id, vertex_set, recurrence_array, recurrence, \
			connectivity, unknown_gene_ratio from %s"%pattern_table)
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		while rows:
			for row in rows:
				mcl_id, vertex_set, recurrence_array, recurrence, connectivity, unknown_ratio = row
				if mcl_id_set and mcl_id not in mcl_id_set:
					continue
				vertex_set = vertex_set[1:-1].split(',')
				vertex_set = map(int, vertex_set)
				vertex_set = dict_map(gene_no2id, vertex_set, type=2)
				recurrence_array = '[' + recurrence_array[1:-1] + ']'
				of.write('[%s, %s, %s, %s, %s, %s],\n'%(mcl_id, repr(vertex_set), recurrence_array,\
					recurrence, connectivity, unknown_ratio))
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		of.write('[]]:\n')	#add the last blank list
		del of
		sys.stderr.write("pattern darwin format done.\n")
		
	def run(self):
		if self.ofname and self.acc_cut_off and self.lm_bit:
			schema_instance = form_schema_tables(self.ofname, self.acc_cut_off, self.lm_bit)
			
		else:
			sys.stderr.write("ofname: %s and acc_cut_off: %s and lm_bit %s, NOT VALID\n"%(self.ofname, self.acc_cut_off, self.lm_bit))
			sys.exit(2)
		conn, curs = db_connect(self.hostname, self.dbname, self.schema)
		#mcl_id_set = self.get_mcl_id_set_from_good_cluster_table(curs, schema_instance.good_cluster_table)
		mcl_id_set = None	#01-14-06
		self._pattern_darwin_format(curs, schema_instance.pattern_table, self.gene_no2id, self.go_no2id, self.output_fname, mcl_id_set)
		del conn, curs

class prediction_darwin_format(Thread):
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, ofname=None, lm_bit='111', acc_cut_off=None, \
		output_fname=None, gene_no2id=None, go_no2id=None, mt_no2tf_name=None, debug=0, report=0):
		Thread.__init__(self)
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.ofname = ofname
		self.lm_bit = lm_bit
		self.acc_cut_off = float(acc_cut_off)
		self.output_fname = output_fname
		self.gene_no2id = gene_no2id
		self.go_no2id = go_no2id
		self.mt_no2tf_name = mt_no2tf_name
		self.debug = int(debug)
		self.report = int(report)
	
	def _prediction_darwin_format(self, curs, p_gene_table, gene_p_table, gene_no2id, go_no2id, output_fname):
		"""
		12-01-05
			deal with lca_list={}
		
		format:
			r:=[
			[gene_id, go_id, is_correct_lca, p_value, mcl_id, lca_list],
			[...],
			[]]:
		"""
		sys.stderr.write("prediction...\n")
		of = open(output_fname, 'w')
		of.write('r:=[\n')
		curs.execute("DECLARE crs CURSOR FOR select p.gene_no, p.go_no, p.is_correct_lca, p.avg_p_value, p.mcl_id, p.lca_list\
			from %s p, %s g where g.p_gene_id = p.p_gene_id"%(p_gene_table, gene_p_table))
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		while rows:
			for row in rows:
				gene_no, go_no, is_correct_lca, p_value, mcl_id, lca_list = row
				if lca_list and len(lca_list)>2:	#12-01-05 lca_list={} just blank
					lca_list = lca_list[1:-1].split(',')
					lca_list = map(int, lca_list)
					lca_list = dict_map(go_no2id, lca_list, type=2)
				else:
					lca_list = []
				of.write("['%s', '%s', %s, %s, %s, %s],\n"%(gene_no2id.get(gene_no) or gene_no, go_no2id[go_no], is_correct_lca,\
					p_value, mcl_id, repr(lca_list)))
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		of.write('[]]:\n')	#add the last blank list
		del of
		sys.stderr.write("prediction darwin format done.\n")
	
	def run(self):
		if self.ofname and self.acc_cut_off and self.lm_bit:
			schema_instance = form_schema_tables(self.ofname, self.acc_cut_off, self.lm_bit)
			
		else:
			sys.stderr.write("ofname: %s and acc_cut_off: %s and lm_bit %s, NOT VALID\n"%(self.ofname, self.acc_cut_off, self.lm_bit))
			sys.exit(2)
		conn, curs = db_connect(self.hostname, self.dbname, self.schema)
		self._prediction_darwin_format(curs, schema_instance.p_gene_table, schema_instance.gene_p_table, \
			self.gene_no2id, self.go_no2id, self.output_fname)
		del conn, curs
	
class Schema2Darwin:
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, ofname=None, lm_bit='111',\
		acc_cut_off=None, output_dir=None, organism=None, running_bit='111', debug=0, report=0):
		"""
		12-01-05
			add running_bit
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.ofname = ofname
		self.lm_bit = lm_bit
		self.acc_cut_off = float(acc_cut_off)
		self.output_dir = output_dir
		self.organism = org_short2long(organism)
		self.running_bit = running_bit
		self.debug = int(debug)
		self.report = int(report)
	
	def run(self):
		"""
		09-28-05
		12-19-05
			use class_list and output_fname_list to ease program writing
		12-30-05
			fix a bug in indexing darwin_instance_list
		"""
		if self.ofname and self.acc_cut_off and self.lm_bit:
			schema_instance = form_schema_tables(self.ofname, self.acc_cut_off, self.lm_bit)
			
			tf_darwin_ofname = os.path.join(self.output_dir, '%s.tf.darwin'%schema_instance.lm_suffix)
			cluster_darwin_ofname = os.path.join(self.output_dir, '%s.cluster.darwin'%schema_instance.lm_suffix)
			prediction_darwin_ofname = os.path.join(self.output_dir, '%s.prediction.darwin'%schema_instance.lm_suffix)
			pattern_darwin_ofname = os.path.join(self.output_dir, '%s.pattern.darwin'%schema_instance.lm_suffix)
			
		else:
			sys.stderr.write("ofname: %s and acc_cut_off: %s and lm_bit %s, NOT VALID\n"%(self.ofname, self.acc_cut_off, self.lm_bit))
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
		
		class_list = [tf_darwin_format, cluster_darwin_format, prediction_darwin_format, pattern_darwin_format]
		output_fname_list = [tf_darwin_ofname, cluster_darwin_ofname, prediction_darwin_ofname, pattern_darwin_ofname]
			#self.ofname is schema table's ofname.
		darwin_instance_list = []
		for i in range(len(self.running_bit)):
			if self.running_bit[i] == '1':
				darwin_instance_list.append(class_list[i](self.hostname, self.dbname, self.schema, self.ofname, self.lm_bit, self.acc_cut_off, \
					output_fname_list[i], gene_id2symbol, go_no2name, mt_no2tf_name, debug, report))
				current_pos = len(darwin_instance_list)-1 #12-30-05
				darwin_instance_list[current_pos].start()
			
		for i in range(len(darwin_instance_list)):
			darwin_instance_list[i].join()


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:f:l:a:o:g:n:br", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	ofname = None
	lm_bit = '111'
	acc_cut_off = 0.6
	output_dir = None
	organism = None
	running_bit = '111'
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
		elif opt in ("-f",):
			ofname = arg
		elif opt in ("-l",):
			lm_bit = arg
		elif opt in ("-a",):
			acc_cut_off = float(arg)
		elif opt in ("-o",):
			output_dir = arg
		elif opt in ("-g",):
			organism = arg
		elif opt in ("-n",):
			running_bit = arg
		elif opt in ("-b",):
			debug = 1
		elif opt in ("-r",):
			report = 1
	if schema and ofname and lm_bit and acc_cut_off and output_dir and organism:
		instance = Schema2Darwin(hostname, dbname, schema, ofname, lm_bit, acc_cut_off, output_dir,\
			organism, running_bit, debug, report)
		instance.run()
		
	else:
		print __doc__
		sys.exit(2)
