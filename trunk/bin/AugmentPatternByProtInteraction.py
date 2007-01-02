#!/usr/bin/env python
"""
Usage: AugmentPatternByProtInteraction.py -k SCHEMA -i xxx -g xx -c xxx -e xxx -p xxx
	-o xxx [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-i ..., input table, either good_cluster or pattern_table
	-o ...,	output table, eg, aug_pi_xxx(running_type=1) or output file(running_type=2)
	-n ...,	prot_interaction_table, 'mrinal_pi.intact_interaction'(default)
	-y ...,	type of input table, 1(default, good_cluster) or 2(pattern_table)
	-p ...,	running type, 1(default), 2, 3, 4 (details see below)
	-x ...,	tax_id, 9606(default)
	-c,	commit
	-b,	debug version.
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	python2.4 ./script/annot/bin/AugmentPatternByProtInteraction.py -k hs_fim_65
	-i good_cl_hs_fim_65_n2s175_m5x65s4l5_ft2_e5_000001a60
	-o aug_pi_hs_fim_65_n2s175_m5x65s4l5_ft2_e5_000001a60 -c -r
	
	AugmentPatternByProtInteraction.py -k hs_fim_65 
		-i good_cl_hs_fim_65_n2s175_m5x65s4l5_ft2_e5_000001a60 
		-o /tmp/good_cl_hs_fim_65_n2s175_m5x65s4l5_ft2_e5_000001a60.prot_int.hg.p_value -r -p 2
	AugmentPatternByProtInteraction.py -k hs_fim_65 -i good_cl_hs_fim_65_n2s175_m5x65s4l5_ft2_e5_000001a60 -o ~/tmp/good_cl_hs_fim_65_n2s175_m5x65s4l5_ft2_e5_000001a60.hiv_int.hg.p_value -r -p 3 -n /usr/local/research_data/ncbi/gene_2006_12_19/GeneRIF/hiv_interactions
	AugmentPatternByProtInteraction.py -k hs_fim_65 -i good_cl_hs_fim_65_n2s175_m5x65s4l5_ft2_e5_000001a60 -o ~/tmp/good_cl_hs_fim_65_n2s175_m5x65s4l5_ft2_e5_000001a60.int.hg.p_value -r -p 4 -n /usr/local/research_data/ncbi/gene_2006_12_19/GeneRIF/interactions
	
Description:
	Running type:
	 1: To find all protein interaction among nodes in the pattern via shortest_path.
	 submit to a new database table
	 2: calculate protein interaction enrichment p-value (hyper-geometric) for all
	 patterns and output to file
	 3: similar to 2 but replacing protein interaction set with hiv interaction set
	 4: similar to 2 but replacing protein interaction set with interaction set from GeneRIF
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
from codense.common import db_connect, pg_1d_array2python_ls, return_vertex_set_string, return_edge_set_string
from sets import Set
import networkx as nx
from DrawPredGOExptTFCompTF_Patterns import DrawPredGOExptTFCompTF_Patterns
import rpy

class AugmentPatternByProtInteraction:
	def __init__(self,hostname='zhoudb', dbname='graphdb', schema=None, \
		input_table=None, output_table=None, prot_interaction_table=None, \
		input_type=1, running_type=1, tax_id=9606, need_commit=0, debug=0, report=0):
		"""
		2006-12-18
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_table = input_table
		self.output_table = output_table
		self.prot_interaction_table = prot_interaction_table
		self.input_type = int(input_type)
		self.running_type = int(running_type)
		self.tax_id = int(tax_id)
		self.need_commit = int(need_commit)
		self.debug = int(debug)
		self.report = int(report)
	
	def augment_pattern_by_prot_interaction(self, vertex_set, prot_interaction_graph):
		g = nx.Graph()
		no_of_nodes = len(vertex_set)
		for m in range(no_of_nodes):
			for n in range(m+1, no_of_nodes):
				u = vertex_set[m]
				v = vertex_set[n]
				if prot_interaction_graph.has_node(u) and prot_interaction_graph.has_node(v):
					shortest_path_list = nx.shortest_path(prot_interaction_graph, u, v)
					if shortest_path_list:	#check the whole shortest path
						for i in range(len(shortest_path_list)-1):
							if not g.has_edge(shortest_path_list[i], shortest_path_list[i+1]):
								g.add_edge(shortest_path_list[i], shortest_path_list[i+1])
		new_vertex_set = []
		for v in g:
			new_vertex_set.append(v)
		new_vertex_set.sort()
		new_edge_set = []
		for (u, v) in g.edges():
			if u<=v:
				edge_tuple = [u, v]
			else:
				edge_tuple = [v, u]
			new_edge_set.append(edge_tuple)
		new_edge_set.sort()
		return new_vertex_set, new_edge_set
	
	def create_aug_pi_table(self, curs, table):
		sys.stderr.write("Creating table %s..."%table)
		curs.execute("create table %s(\
			id	serial,\
			mcl_id	integer primary key,\
			vertex_set	integer[],\
			edge_set	integer[][])"%table)
		sys.stderr.write("done.\n")
	
	def submit_aug_pi_graph(self, curs, table, mcl_id, vertex_set, edge_set):
		if self.debug:
			sys.stderr.write("\nSubmitting augmented graph to %s..."%table)
		curs_sentence = "insert into %s(mcl_id, vertex_set, edge_set) values(%s, '%s', '%s')"%\
			(table, mcl_id, return_vertex_set_string(vertex_set), return_edge_set_string(edge_set))
		curs.execute(curs_sentence)
		if self.debug:
			sys.stderr.write("submission done.\n")
	
	def get_total_vertex_set_from_patterns(self, curs, input_table):
		"""
		2006-12-28
			get a distinct set of vertices from all patterns for hyper-geometric test
		"""
		sys.stderr.write("Getting total vertex_set from %s...\n"%input_table)
		curs.execute("select distinct g.gene_no from gene g, %s c where g.gene_no=any(c.vertex_set)"%input_table)
		rows = curs.fetchall()
		total_vertex_set = Set()
		for row in rows:
			gene_no = row[0]
			total_vertex_set.add(gene_no)
		sys.stderr.write("done.\n")
		return total_vertex_set
	
	def get_hiv_interaction_set(self, prot_interaction_table, tax_id):
		"""
		2006-12-28
			the data is from NCBI Gene / GeneRIF / hiv_interactions.gz
			
		"""
		sys.stderr.write("Getting HIV interaction set ...")
		reader = csv.reader(open(prot_interaction_table), delimiter='\t')
		hiv_interaction_set = Set()
		for row in reader:
			tax_id1, gene_id1, acc_ver1, product_name1, interaction_phrase, tax_id2, \
				gene_id2, acc_ver2, product_name2, pmid_list, update_timestamp, generif_text = row
			tax_id1 = int(tax_id1)
			tax_id2 = int(tax_id2)
			if tax_id1==tax_id:
				hiv_interaction_set.add(int(gene_id1))
			if tax_id2==tax_id:
				hiv_interaction_set.add(int(gene_id2))
		del reader
		sys.stderr.write("done.\n")
		return hiv_interaction_set
	
	def get_interaction_set(self, prot_interaction_table, tax_id):
		"""
		2006-12-28
			the data is from NCBI Gene / GeneRIF / interactions.gz
		"""
		sys.stderr.write("Getting interaction set ...")
		reader = csv.reader(open(prot_interaction_table), delimiter='\t')
		interaction_set = Set()
		reader.next()
		for row in reader:
			tax_id1, gene_id1, acc_ver1, name1, keyphrase1, tax_id2, interactant_id, interactant_id_type, acc_ver2, \
			name2, complex_id, complex_id_type, complex_name, pubmed_id_list, last_mod, generif_text, \
			interaction_id, interaction_id_type = row
			tax_id1 = int(tax_id1)
			if tax_id1==tax_id:
				interaction_set.add(int(gene_id1))
			if tax_id2!='-':
				tax_id2 = int(tax_id2)
				if tax_id2==tax_id and interactant_id_type=='GeneID':
					interaction_set.add(int(interactant_id))
		del reader
		sys.stderr.write("done.\n")
		return interaction_set
	
	def cal_pi_hg_p_value(self, vertex_set, prot_int_vertex_set, total_vertex_set):
		"""
		2006-12-28
			hyper-geometric p-value, taking log
			prot_int_vertex_set doesn't contain vertices outside total_vertex_set
		"""
		if self.debug:
			sys.stderr.write("\nCalculating protein interaction enrichment p-value...")
		vertex_set = Set(vertex_set)
		x = len(vertex_set&prot_int_vertex_set)
		n = len(prot_int_vertex_set)
		m = len(total_vertex_set) - n
		k = len(vertex_set)
		if self.debug:
			sys.stderr.write("Done.\n")
		return rpy.r.phyper(x, n, m, k, lower_tail=rpy.r.FALSE, log_p=rpy.r.TRUE)
	
	def batch_augment_all_patterns(self, curs, output_table, prot_interaction_graph, \
			running_type, total_vertex_set, prot_int_vertex_set):
		"""
		2006-12-28
			add running_type==2
		"""
		sys.stderr.write("Augment patterns ...\n")
		curs.execute("fetch 1000 from crs")
		rows = curs.fetchall()
		counter = 0
		log_p_value_mcl_id_list = []
		while rows:
			for row in rows:
				mcl_id, vertex_set = row
				vertex_set = vertex_set[1:-1].split(',')
				vertex_set = map(int, vertex_set)
				if running_type==1:
					new_vertex_set, new_edge_set = self.augment_pattern_by_prot_interaction(vertex_set, prot_interaction_graph)
					self.submit_aug_pi_graph(curs, output_table, mcl_id, new_vertex_set, new_edge_set)
				elif running_type in [2,3,4]:
					log_p_value = self.cal_pi_hg_p_value(vertex_set, prot_int_vertex_set, total_vertex_set)
					log_p_value_mcl_id_list.append([log_p_value, mcl_id])
				counter += 1
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, counter))
			curs.execute("fetch 1000 from crs")
			rows = curs.fetchall()
		
		sys.stderr.write("done.\n")
		return log_p_value_mcl_id_list
	
	def output_to_file(self, log_p_value_mcl_id_list, output_table):
		sys.stderr.write("Outputting log_p_value_mcl_id_list...")
		writer = csv.writer(open(output_table, 'w'), delimiter='\t')
		log_p_value_mcl_id_list.sort()
		for row in log_p_value_mcl_id_list:
			writer.writerow(row)
		del writer
		sys.stderr.write("done.\n")
	
	def run(self):
		"""
		2006-12-25
			
		--db_connect()
		--create_aug_pi_table()
		--get_prot_interaction_graph()
		--batch_augment_all_patterns()
			--augment_pattern_by_prot_interaction()
			--submit_aug_pi_graph()
		"""
		(conn, curs) = db_connect(self.hostname, self.dbname, self.schema)
		
		if self.running_type==1:
			self.create_aug_pi_table(curs, self.output_table)
			total_vertex_set = None	#2006-12-28 set this to None
		elif self.running_type in [2,3,4]:
			total_vertex_set = self.get_total_vertex_set_from_patterns(curs, self.input_table)
		else:
			sys.stderr.write("Unsupported running_type=%s.\n"%self.running_type)
		
		if self.input_type==1:
			curs.execute("DECLARE crs CURSOR FOR select mcl_id, vertex_set from %s"%self.input_table)
		elif self.input_type==2:
			curs.execute("DECLARE crs CURSOR FOR select id, vertex_set from %s"%self.input_table)
		
		if self.running_type in [1,2]:
			DrawPredGOExptTFCompTF_Patterns_instance = DrawPredGOExptTFCompTF_Patterns()
			prot_interaction_graph = DrawPredGOExptTFCompTF_Patterns_instance.get_prot_interaction_graph(\
				curs, self.prot_interaction_table, self.tax_id)
		else:
			prot_interaction_graph = None
		
		if self.running_type == 2:
			prot_int_vertex_set = Set(prot_interaction_graph.nodes())&total_vertex_set
		elif self.running_type==3:
			hiv_interaction_set = self.get_hiv_interaction_set(self.prot_interaction_table, self.tax_id)
			prot_int_vertex_set = hiv_interaction_set&total_vertex_set
		elif self.running_type==4:
			interaction_set = self.get_interaction_set(self.prot_interaction_table, self.tax_id)
			prot_int_vertex_set = interaction_set&total_vertex_set
		else:
			prot_int_vertex_set = None	#2006-12-28 set this to None
		
		log_p_value_mcl_id_list = self.batch_augment_all_patterns(curs, self.output_table, prot_interaction_graph, \
			self.running_type, total_vertex_set, prot_int_vertex_set)
		
		curs.execute("close crs")
		if self.need_commit:
			curs.execute("end")
		
		if self.running_type in [2,3,4]:
			self.output_to_file(log_p_value_mcl_id_list, self.output_table)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:o:n:y:p:x:cbr", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	input_table = None
	output_table = None
	prot_interaction_table = 'mrinal_pi.intact_interaction'
	input_type = 1
	running_type = 1
	tax_id = 9606
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
		elif opt in ("-i",):
			input_table = arg
		elif opt in ("-o",):
			output_table = arg
		elif opt in ("-n",):
			prot_interaction_table = arg
		elif opt in ("-y",):
			input_type = int(arg)
		elif opt in ("-p",):
			running_type = int(arg)
		elif opt in ("-x",):
			tax_id = int(arg)
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-b",):
			debug = 1
		elif opt in ("-r",):
			report = 1
	if schema and input_table and output_table:
			instance = AugmentPatternByProtInteraction(hostname, dbname, schema, \
				input_table, output_table, prot_interaction_table, input_type, \
				running_type, tax_id, commit, debug, report)
			instance.run()
	else:
		print __doc__
		sys.exit(2)