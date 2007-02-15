#!/usr/bin/env python
"""
Usage: Schema2Darwin.py -k SCHEMA -p -s -i  -o OUTPUT_DIR [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-p ...,	pattern_table
	-s ...,	cluster_bs_table
	-i ...,	prediction input file, p_gene_table(2007-01-04)
	-l ...,	lm_bit, 111(default, IGNORE), gene_p_table(2007-01-04)
	-a ...,	accuracy cutoff used to prediction, 0.6(default, IGNORE)
	-o ...,	output_dir
	-g ...,	organism, for gene_id2symbol
	-n,	running bit, (tf_darwin_format, cluster_darwin_format, prediction_darwin_format, 
		pattern_darwin_format), 1111(default)
	-b,	debug version.
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	ssh $HOSTNAME Schema2Darwin.py -k mm_fim_97 -p -s -i mm_fim_97m4x200
		-o /tmp/yuhuang/ -g mm
	
	Schema2Darwin.py -k hs_fim_65 -p pattern_hs_fim_65_n2s175_m5x65s4l5 -s bs_hs_fim_65_n2s175_m5x65s4l5_ft2_e5_000001a60p01y1fe0_1 -i p_gene_hs_fim_65_n2s175_m5x65s4l5_ft2_e5 -l gene_p_hs_fim_65_n2s175_m5x65s4l5_ft2_e5_000001a60 -o /tmp/yuhuang -g hs -n 0011

	Schema2Darwin.py -k hs_fim_65 -p pattern_hs_fim_65_n2s175_m5x65s4l5bug -s good_cl_hs_fim_65_n2s175_m5x65s4l5bug_ft2_e5_000001a60 -i p_gene_hs_fim_65_n2s175_m5x65s4l5bug_ft2_e5 -l gene_p_hs_fim_65_n2s175_m5x65s4l5bug_ft2_e5_000001a60  -o /tmp/yuhuang -g hs -n 00001
	
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

class darwin_format:
	def __init__(self, *arguments):
		self.hostname, self.dbname, self.schema, self.pattern_table, self.cluster_bs_table, self.input_fname, \
		self.lm_bit, self.acc_cut_off, self.output_fname,	\
		self.gene_no2id, self.go_no2id, self.mt_no2tf_name, self.debug, self.report = arguments
		self.acc_cut_off = float(self.acc_cut_off)
		self.debug = int(self.debug)
		self.report = int(self.report)
	
class tf_darwin_format(Thread, darwin_format):
	"""
	2006-09-25
		use class inheritance
	"""
	def __init__(self, *arguments):
		Thread.__init__(self)
		darwin_format.__init__(self, *arguments)
	
	def _tf_darwin_format(self, curs, good_cluster_table, output_fname, gene_no2id, mcl_id2tf_set):
		"""
		2006-09-25
			change good_cluster_table to be pattern_table
			
		format:
			r:=[
			[mcl_id, [gene1, gene2, ...], [ [TF1], [hyper_p_value] ], [ [TF2, TF3], [hyper_p_value] ], ... ],
			[...],
			[]]:
		"""
		sys.stderr.write("TF...\n")
		of = open(output_fname, 'w')
		of.write('r:=[\n')
		curs.execute("DECLARE crs CURSOR FOR select id, vertex_set from %s"%good_cluster_table)	#2006-09-25
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
		"""
		2006-09-25
			use self.cluster_bs_table and self.pattern_table
		"""
		conn, curs = db_connect(self.hostname, self.dbname, self.schema)
		mcl_id2tf_set = get_mcl_id2tf_set(curs, self.cluster_bs_table, self.mt_no2tf_name)
		self._tf_darwin_format(curs, self.pattern_table, self.output_fname, self.gene_no2id, mcl_id2tf_set)
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


class pattern_darwin_format(Thread, darwin_format):
	"""
	12-19-05
		output all patterns from pattern_table into darwin format
	2006-09-25
		use class inheritance
	"""
	def __init__(self, *arguments):
		Thread.__init__(self)
		darwin_format.__init__(self, *arguments)

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
		2007-01-07
			add edge_set
		
		format:
			r:=[
			[mcl_id, vertex_set, edge_set, recurrence_array, recurrence, connectivity, unknown_ratio]
			[...],
			...
			[]]:
			
		"""
		sys.stderr.write("pattern...\n")
		of = open(output_fname, 'w')
		of.write('r:=[\n')
		curs.execute("DECLARE crs CURSOR FOR select id, vertex_set, edge_set, recurrence_array, recurrence, \
			connectivity, unknown_gene_ratio from %s"%pattern_table)
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		while rows:
			for row in rows:
				mcl_id, vertex_set, edge_set, recurrence_array, recurrence, connectivity, unknown_ratio = row
				if mcl_id_set and mcl_id not in mcl_id_set:
					continue
				vertex_set = vertex_set[1:-1].split(',')
				vertex_set = map(int, vertex_set)
				vertex_set = dict_map(gene_no2id, vertex_set, type=2)
				edge_set = edge_set[2:-2].split('},{')
				for i in range(len(edge_set)):
					edge = edge_set[i].split(',')
					edge = map(int, edge)
					edge = dict_map(gene_no2id, edge, type=2)
					edge_set[i] = edge
				recurrence_array = '[' + recurrence_array[1:-1] + ']'
				of.write('[%s, %s, %s, %s, %s, %s, %s],\n'%(mcl_id, repr(vertex_set), repr(edge_set), recurrence_array,\
					recurrence, connectivity, unknown_ratio))
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		of.write('[]]:\n')	#add the last blank list
		del of
		sys.stderr.write("pattern darwin format done.\n")
		
	def run(self):
		"""
		2006-09-25
			use self.pattern_table
		"""
		conn, curs = db_connect(self.hostname, self.dbname, self.schema)
		#mcl_id_set = self.get_mcl_id_set_from_good_cluster_table(curs, schema_instance.good_cluster_table)
		mcl_id_set = None	#01-14-06
		self._pattern_darwin_format(curs, self.pattern_table, self.gene_no2id, self.go_no2id, self.output_fname, mcl_id_set)
		del conn, curs

class prediction_darwin_format(Thread, darwin_format):
	def __init__(self, *arguments):
		Thread.__init__(self)
		darwin_format.__init__(self, *arguments)
	
	def _prediction_darwin_format(self, curs, p_gene_table, gene_p_table, gene_no2id, go_no2id, output_fname):
		"""
		12-01-05
			deal with lca_list={}
		03-01-06
			add no_of_distinct_funcitons_from_gene_p_table in the output
		2006-09-25 now defunct
		
		format:
			r:=[
			[gene_id, go_id, is_correct_lca, p_value, mcl_id, lca_list, no_of_distinct_funcitons_from_gene_p_table],
			[...],
			[]]:
		"""
		sys.stderr.write("prediction...\n")
		#03-01-06 firstly get the gene_no2p_gene_id_src_set
		curs.execute("DECLARE crs_1 CURSOR FOR SELECT p.gene_no, g.p_gene_id_src from %s p, %s g\
			where p.p_gene_id=g.p_gene_id_src"%(p_gene_table, gene_p_table))
		curs.execute("fetch 5000 from crs_1")
		rows = curs.fetchall()
		gene_no2p_gene_id_src_set = {}
		while rows:
			for row in rows:
				gene_no, p_gene_id_src = row
				if gene_no not in gene_no2p_gene_id_src_set:
					gene_no2p_gene_id_src_set[gene_no] = Set()
				gene_no2p_gene_id_src_set[gene_no].add(p_gene_id_src)
			curs.execute("fetch 5000 from crs_1")
			rows = curs.fetchall()
		curs.execute("close crs_1")
		
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
				of.write("['%s', '%s', %s, %s, %s, %s, %s],\n"%(gene_no2id.get(gene_no) or gene_no, go_no2id[go_no], is_correct_lca,\
					p_value, mcl_id, repr(lca_list), len(gene_no2p_gene_id_src_set[gene_no])))	#03-01-06
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		of.write('[]]:\n')	#add the last blank list
		del of
		curs.execute("close crs")
		sys.stderr.write("prediction darwin format done.\n")
	
	def _prediction_darwin_format_from_file(self, input_fname, gene_no2id, go_no2id, output_fname):
		"""
		2006-09-25
		format:
			r:=[
			[gene_id, pattern_id, go_id, is_known, is_correct, markov_score, p_value],
			[...],
			[]]:
		"""
		sys.stderr.write("prediction...\n")
		of = open(output_fname, 'w')
		of.write('r:=[\n')
		reader = csv.reader(open(input_fname, 'r'), delimiter='\t')
		for row in reader:
			gene_id, pattern_id, go_id, is_known, is_correct, markov_score, p_value = row
			gene_no = int(gene_id)
			of.write("['%s', %s, '%s', %s, %s, %s, %s],\n"%(gene_no2id.get(gene_no) or gene_no, pattern_id, \
				go_id, is_known, is_correct, markov_score, p_value))
		del reader
		of.write('[]]:\n')	#add the last blank list
		del of
		sys.stderr.write("prediction darwin format done.\n")
		
	def run(self):
		"""
		2006-09-25
			use self.input_fname
		"""
		conn, curs = db_connect(self.hostname, self.dbname, self.schema)
		self._prediction_darwin_format(curs, self.input_fname, self.lm_bit, \
			self.gene_no2id, self.go_no2id, self.output_fname)
		#self._prediction_darwin_format_from_file(self.input_fname, self.gene_no2id, self.go_no2id, self.output_fname)
		del conn, curs

class context_prediction_csv_format(Thread, darwin_format):
	"""
	2007-02-08
	"""
	def __init__(self, *arguments):
		Thread.__init__(self)
		darwin_format.__init__(self, *arguments)
	
	def get_mcl_id_set_from_good_cluster_table(self, curs, good_cluster_table):
		"""
		2007-02-08
			copied from pattern_darwin_format
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
	
	def get_mcl_id2vertex_edge_recurrence(self, curs, pattern_table, gene_no2id, go_no2id, mcl_id_set):
		"""
		2007-02-08
		
		"""
		sys.stderr.write("Getting mcl_id2vertex_edge_recurrence ...\n")
		mcl_id2vertex_edge_recurrence = {}
		curs.execute("DECLARE crs CURSOR FOR select id, vertex_set, edge_set, recurrence_array from %s"%pattern_table)
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		recurrence_func = lambda x: int(float(x)>=0.8)
		while rows:
			for row in rows:
				mcl_id, vertex_set, edge_set, recurrence_array = row
				if mcl_id in mcl_id_set:
					vertex_set = vertex_set[1:-1].split(',')
					vertex_set = map(int, vertex_set)
					vertex_set = dict_map(gene_no2id, vertex_set, type=2)
					edge_set = edge_set[2:-2].split('},{')
					for i in range(len(edge_set)):
						edge = edge_set[i].split(',')
						edge = map(int, edge)
						edge = dict_map(gene_no2id, edge, type=2)
						edge_set[i] = edge
					recurrence_array = recurrence_array[1:-1].split(',')
					recurrence_array = map(recurrence_func, recurrence_array)
					new_recurrence_array = []
					for i in range(len(recurrence_array)):
						if recurrence_array[i]==1:
							new_recurrence_array.append(i+1)
					mcl_id2vertex_edge_recurrence[mcl_id] = [vertex_set, edge_set, new_recurrence_array]
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		curs.execute("close crs")
		sys.stderr.write("done.\n")
		return mcl_id2vertex_edge_recurrence
	
	def _prediction_csv_format(self, curs, p_gene_table, gene_p_table, gene_no2id, go_no2id, output_fname, mcl_id2vertex_edge_recurrence):
		"""
		2007-02-08
			copied from prediction_darwin_format
		"""
		sys.stderr.write("prediction...\n")
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		header = ['gene', 'go_function', 'is_correct',\
			'hyper-geometric p_value', 'pattern id', 'recurrence', 'dataset id list', 'connectivity',\
			'cluster_size', 'unknown_ratio', 'avg_degree', 'network_topo_score', 'network neighbors',\
			'network topology (edge_set)' ]
		writer.writerow(header)
		curs.execute("DECLARE crs CURSOR FOR select p.gene_no, p.go_no, p.is_correct_lca, p.avg_p_value, p.mcl_id, p.recurrence_cut_off,\
			p.connectivity_cut_off, p.cluster_size_cut_off, p.unknown_cut_off, p.vertex_gradient, p.edge_gradient\
			from %s p, %s g where g.p_gene_id = p.p_gene_id order by gene_no, go_no"%(p_gene_table, gene_p_table))
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		while rows:
			for row in rows:
				gene_no, go_no, is_correct_lca, p_value, mcl_id, recurrence, connectivity, cluster_size, unknown_ratio, avg_degree, network_topo_score = row
				
				writer.writerow([gene_no2id.get(gene_no) or gene_no, go_no2id[go_no], is_correct_lca,\
					p_value, mcl_id, recurrence, mcl_id2vertex_edge_recurrence[mcl_id][2], connectivity,\
					cluster_size, unknown_ratio, avg_degree, network_topo_score, mcl_id2vertex_edge_recurrence[mcl_id][0],\
					mcl_id2vertex_edge_recurrence[mcl_id][1] ])
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		del writer
		curs.execute("close crs")
		sys.stderr.write("prediction csv format done.\n")
		
	def run(self):
		"""
		2007-02-08
		"""
		conn, curs = db_connect(self.hostname, self.dbname, self.schema)
		mcl_id_set = self.get_mcl_id_set_from_good_cluster_table(curs, self.cluster_bs_table)
		mcl_id2vertex_edge_recurrence = self.get_mcl_id2vertex_edge_recurrence(curs, self.pattern_table, self.gene_no2id, self.go_no2id, mcl_id_set)
		self._prediction_csv_format(curs, self.input_fname, self.lm_bit, \
			self.gene_no2id, self.go_no2id, self.output_fname, mcl_id2vertex_edge_recurrence)
		del conn, curs

class Schema2Darwin:
	def __init__(self, *arguments):
		"""
		12-01-05
			add running_bit
		2006-09-25 use pattern_table, cluster_bs_table, input_fname
			use *arguments
		"""
		self.hostname, self.dbname, self.schema, self.pattern_table, self.cluster_bs_table, self.input_fname, \
		self.lm_bit, self.acc_cut_off, self.output_dir, self.organism, self.running_bit, self.debug, self.report = arguments
		
		self.organism = org_short2long(self.organism)
		self.acc_cut_off = float(acc_cut_off)
		self.debug = int(debug)
		self.report = int(report)
		
	
	def replace_prime_in_gene_id2symbol(self, gene_id2symbol):
		"""
		01-26-06
			some symbols contain ', which darwin can't handle
		"""
		for gene_id, symbol in gene_id2symbol.iteritems():
			symbol = symbol.replace("'", '\prime')
			gene_id2symbol[gene_id] = symbol
		return gene_id2symbol
	
	def run(self):
		"""
		09-28-05
		12-19-05
			use class_list and output_fname_list to ease program writing
		12-30-05
			fix a bug in indexing darwin_instance_list
		2006-09-25
		2007-02-08
			add context_prediction_csv_format
		"""
		tf_darwin_ofname = os.path.join(self.output_dir, '%s.tf.darwin'%self.cluster_bs_table)
		cluster_darwin_ofname = os.path.join(self.output_dir, '%s.cluster.darwin'%os.path.basename(self.input_fname))
		prediction_darwin_ofname = os.path.join(self.output_dir, '%s.prediction.darwin'%os.path.basename(self.input_fname))
		pattern_darwin_ofname = os.path.join(self.output_dir, '%s.pattern.darwin'%self.pattern_table)
		
		if not os.path.isdir(self.output_dir):
			os.makedirs(self.output_dir)
		conn, curs = db_connect(self.hostname, self.dbname, self.schema)
		
		tax_id = org2tax_id(self.organism)
		#gene_no2id = get_gene_no2gene_id(curs)	#Watch, if unigene, should use this.
		gene_id2symbol = get_gene_id2gene_symbol(curs, tax_id)
		
		gene_id2symbol = self.replace_prime_in_gene_id2symbol(gene_id2symbol)	#01-26-06
		
		#gene_no2symbol = dict_transfer(gene_no2id, gene_id2symbol)
		#Jasmine wants the gene symbol 09-28-05
		#gene_id is integer in gene.gene table and same as gene_no, so just use it.
		go_no2name = get_go_no2name(curs)	#09-28-05 Jasmine wants the go_name, not go_id
		
		#2006-09-25 use gene_id2symbol to replace mt_no2tf_name
		#mt_no2tf_name = get_mt_no2tf_name()
		mt_no2tf_name = gene_id2symbol
		
		class_list = [tf_darwin_format, cluster_darwin_format, prediction_darwin_format, pattern_darwin_format, context_prediction_csv_format]
		context_prediction_csv_fname = os.path.join(self.output_dir, '%s.context.csv'%self.input_fname)
		output_fname_list = [tf_darwin_ofname, cluster_darwin_ofname, prediction_darwin_ofname, pattern_darwin_ofname, context_prediction_csv_fname]
		darwin_instance_list = []
		for i in range(len(self.running_bit)):
			if self.running_bit[i] == '1':
				darwin_instance_list.append(class_list[i](self.hostname, self.dbname, self.schema, self.pattern_table,\
					self.cluster_bs_table, self.input_fname, self.lm_bit, self.acc_cut_off, \
					output_fname_list[i], gene_id2symbol, go_no2name, mt_no2tf_name, debug, report))	#2006-09-25
				current_pos = len(darwin_instance_list)-1 #12-30-05
				darwin_instance_list[current_pos].start()
			
		for i in range(len(darwin_instance_list)):
			darwin_instance_list[i].join()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:p:s:i:l:a:o:g:n:br", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	pattern_table = ''
	cluster_bs_table = ''
	input_fname = ''
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
		elif opt in ("-p",):
			pattern_table = arg
		elif opt in ("-s",):
			cluster_bs_table = arg
		elif opt in ("-i",):
			input_fname = arg
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
	if schema and cluster_bs_table and pattern_table and input_fname and lm_bit and acc_cut_off and output_dir and organism:
		instance = Schema2Darwin(hostname, dbname, schema, pattern_table, cluster_bs_table, input_fname, lm_bit, \
			acc_cut_off, output_dir, organism, running_bit, debug, report)
		instance.run()
		
	else:
		print __doc__
		sys.exit(2)
