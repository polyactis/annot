#!/usr/bin/env python
"""
Usage: gene_p_translate.py -k SCHEMA -t P_GENE_TABLE -n GENE_P_TABLE -m MCL_TABLE [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --p_gene_table=...	the p_gene table
	-n ..., --gene_p_table=...	update the gene_p table, needed if needcommit
	-m ..., --mcl_table=...	the mcl_table
	-j ..., --judger_type=...	how to judge predicted functions, 0(default), 1, 2
	-c, --commit	commit this database transaction (IGNORE)
	-r, --report	report flag
	-u, --debug debug flag
	-h, --help              show this help

Examples:
	gene_p_translate.py -k sc_54 -t p_gene_repos_2_e5 -m mcl_result_repos_2
		-n gene_p_repos_2_e5 >/tmp/gene_p_translate.out
	
Description:
	03-08-05
	it simply translates the gene-p table into a human-readable table.

	
"""

import sys, os, getopt, csv
from codense.common import db_connect, get_gene_no2direct_go, get_go_no2name
from codense.common import get_gene_no2gene_id, get_go_no2go_id, dict_map
from p_gene_analysis import accuracy_struc
from numarray import array, greater
from sets import Set

class function_struc:
	#data structure for p_functions_struc_dict in gene_prediction
	def __init__(self):
		self.is_correct = 0
		self.is_correct_L1 = 0
		self.is_correct_L1 = 0
		self.p_value_list = []
		self.cluster_array = []
		self.context_dict = {}
		
class gene_p_translate:
	"""
	03-03-05
		run()
			--db_connect()
			--data_fetch()
				--two_maps_setup()
			--output()
				(loop begins)
				--output_one_gene()
				--output_function_group()
	"""
	def __init__(self, hostname=None, dbname=None, schema=None, p_gene_table=None, \
		gene_p_table=None, mcl_table=None, judger_type=0, needcommit=0, report=0, debug=0):
		
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.p_gene_table = p_gene_table
		self.gene_p_table = gene_p_table
		self.mcl_table = mcl_table
		self.judger_type = int(judger_type)
		self.needcommit = int(needcommit)
		self.report = int(report)		
		self.debug = int(debug)
		
		self.known_gene_no2p_gene_id_src = {}		#key is gene_no, value is p_gene_id_src list
		self.unknown_gene_no2p_gene_id_src = {}
		self.p_gene_id_src_map = {}	#key is p_gene_id_src, value is a dictionary of function_struc(a function group regarded as one predicted function)
		
		self.is_correct_dict = {0: 'is_correct',
			1: 'is_correct_L1',
			2: 'is_correct_lca'}
	
	def data_fetch(self, curs, p_gene_table, gene_p_table, mcl_table):
		"""
		03-02-05
			borrowed from p_gene_lm.py
			
			--two_maps_setup()
		"""
		sys.stderr.write("Setting up two maps...\n")
		curs.execute("DECLARE crs CURSOR FOR select g.p_gene_id_src, p.p_gene_id, p.gene_no, \
			p.go_no, p.is_correct, p.is_correct_L1, p.is_correct_lca, p.avg_p_value, p.mcl_id, m.vertex_set\
			from %s g, %s p, %s m where g.p_gene_id=p.p_gene_id and p.mcl_id=m.mcl_id"%\
			(gene_p_table, p_gene_table, mcl_table))
		no_of_records = 0
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		while rows:
			for row in rows:
				self.two_maps_setup(row)
				no_of_records += 1
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, no_of_records))
			
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		sys.stderr.write("Done\n")
	
	def two_maps_setup(self, row):
		"""
		03-03-05
		"""
		p_gene_id_src = row[0]
		p_gene_id = row[1]
		gene_no = row[2]
		go_no = row[3]
		is_correct = row[4]
		is_correct_L1= row[5]
		is_correct_lca = row[6]
		p_value = row[7]
		mcl_id = row[8]
		vertex_set = row[9]
		vertex_set = vertex_set[1:-1].split(',')
		vertex_set = map(int, vertex_set)
		if is_correct == -1:
			if gene_no not in self.unknown_gene_no2p_gene_id_src:
				self.unknown_gene_no2p_gene_id_src[gene_no] = Set([p_gene_id_src])
			else:
				self.unknown_gene_no2p_gene_id_src[gene_no].add(p_gene_id_src)		
		else:
			if gene_no not in self.known_gene_no2p_gene_id_src:
				self.known_gene_no2p_gene_id_src[gene_no] = Set([p_gene_id_src])
			else:
				self.known_gene_no2p_gene_id_src[gene_no].add(p_gene_id_src)
		
		if p_gene_id_src not in self.p_gene_id_src_map:
			self.p_gene_id_src_map[p_gene_id_src] = {}

		if go_no not in self.p_gene_id_src_map[p_gene_id_src]:
			self.p_gene_id_src_map[p_gene_id_src][go_no] = function_struc()
		#pass it to ease programming
		unit = self.p_gene_id_src_map[p_gene_id_src][go_no]
		unit.is_correct = is_correct
		unit.is_correct_L1 = is_correct_L1
		unit.is_correct_lca = is_correct_lca
		unit.p_value_list.append(p_value)
		unit.cluster_array.append(mcl_id)
		for vertex in vertex_set:
			if vertex in unit.context_dict:
				unit.context_dict[vertex] += 1
			else:
				unit.context_dict[vertex] = 1
	
	def get_go_no2accuracy(self, curs, p_gene_table, gene_p_table):
		"""
		03-02-05
			return (go_no2accuracy, go_no2accuracy_pair)
			
		"""
		sys.stderr.write("Getting go_no2accuracy, go_no2accuracy_pair...")
		go_no2accuracy = {}
		go_no2accuracy_pair = {}
		gene_no_go_no_pair_dict = {}
		curs.execute("select g.p_gene_id, p.gene_no, p.go_no, p.%s \
			from %s g, %s p where g.p_gene_id=p.p_gene_id"%\
			(self.is_correct_dict[self.judger_type], gene_p_table, p_gene_table))
		rows = curs.fetchall()
		for row in rows:
			gene_no = row[1]
			go_no = row[2]
			is_correct = row[3]
			if go_no not in go_no2accuracy:
				go_no2accuracy[go_no] = accuracy_struc()
			go_no2accuracy[go_no].correct += (is_correct+1)/2	#increment when is_correct=1
			go_no2accuracy[go_no].known += int(is_correct>=0)	#increment when is_correct=0 or 1
			go_no2accuracy[go_no].unknown += -(is_correct-1)/2	#increment when is_correct = -1
			if (gene_no, go_no) not in gene_no_go_no_pair_dict:
				gene_no_go_no_pair_dict[(gene_no, go_no)] = 1
				if go_no not in go_no2accuracy_pair:
					go_no2accuracy_pair[go_no] = accuracy_struc()
				go_no2accuracy_pair[go_no].correct += (is_correct+1)/2	#increment when is_correct=1
				go_no2accuracy_pair[go_no].known += int(is_correct>=0)	#increment when is_correct=0 or 1
				go_no2accuracy_pair[go_no].unknown += -(is_correct-1)/2	#increment when is_correct = -1
		
		for go_no in go_no2accuracy:
			unit = go_no2accuracy[go_no]
			unit1 = go_no2accuracy_pair[go_no]
			if unit.known !=0:
				unit.ratio = unit.correct/unit.known
			else:
				unit.ratio = -1
			if unit1.known !=0:
				unit1.ratio = unit1.correct/unit1.known
			else:
				unit1.ratio = -1
		sys.stderr.write("Done\n")
		return (go_no2accuracy, go_no2accuracy_pair)

	def output(self, curs, outf, known_gene_no2p_gene_id_src, unknown_gene_no2p_gene_id_src, p_gene_id_src_map):
		"""
		03-03-05
			loop over gene_no2p_gene_id_src and p_gene_id_src_map
		03-13-05
			add a column, #clusters in the output file
		"""
		#three dictionaries
		gene_no2gene_id = get_gene_no2gene_id(curs)
		gene_no2direct_go = get_gene_no2direct_go(curs)
		go_no2go_id = get_go_no2go_id(curs)
		go_no2name = get_go_no2name(curs)
		go_no2accuracy, go_no2accuracy_pair = self.get_go_no2accuracy(curs, self.p_gene_table, self.gene_p_table)
		
		sys.stderr.write("Outputing prediction table...")
		writer = csv.writer(outf, delimiter='\t')
		#first output the known genes
		for (gene_no, p_gene_id_src_list) in known_gene_no2p_gene_id_src.iteritems():
			self.output_one_gene(curs, writer, gene_no, gene_no2gene_id, gene_no2direct_go)
			row = ['go_no', 'go_id', 'go_name', 'is_correct', 'is_correct_L1', 'is_correct_lca', 'p_value_list', '#clusters', 'mcl_id_list', \
				'e_acc', 'e_acc_pair', 'cluster_context']
			writer.writerow(row)
			for p_gene_id_src in p_gene_id_src_list:
				self.output_function_group(curs, writer, p_gene_id_src_map[p_gene_id_src], gene_no2gene_id,\
					go_no2go_id, go_no2name, go_no2accuracy, go_no2accuracy_pair)
			writer.writerow([])
		#second output the unknown genes
		for (gene_no, p_gene_id_src_list) in unknown_gene_no2p_gene_id_src.iteritems():
			self.output_one_gene(curs, writer, gene_no, gene_no2gene_id, gene_no2direct_go)
			row = ['go_no', 'go_id', 'go_name', 'is_correct', 'is_correct_L1', 'is_correct_lca', 'p_value_list', '#clusters', 'mcl_id_list', \
				'e_acc', 'e_acc_pair', 'cluster_context']
			writer.writerow(row)
			for p_gene_id_src in p_gene_id_src_list:
				self.output_function_group(curs, writer, p_gene_id_src_map[p_gene_id_src], gene_no2gene_id,\
					go_no2go_id, go_no2name, go_no2accuracy, go_no2accuracy_pair)
			writer.writerow([])
		del writer
		sys.stderr.write("Done\n")

	def output_one_gene(self, curs, writer, gene_no, gene_no2gene_id, gene_no2direct_go):
		"""
		03-02-05
			output one row information of a gene.
		"""
		row = ['gene_no', 'gene_id', 'direct functions', 'all functions']
		writer.writerow(row)
		row = [gene_no, gene_no2gene_id[gene_no]]
		direct_go_id_list = gene_no2direct_go.get(gene_no)
		if direct_go_id_list:
			row.append(';'.join(direct_go_id_list))
		else:
			row.append('')
		
		curs.execute("select g.gene_no,go.go_id from gene g, go where go.go_no=any(g.go_functions) and g.gene_no=%d"%(gene_no))
		rows = curs.fetchall()
		go_functions_list = []
		for sub_row in rows:
			go_functions_list.append(sub_row[1])
		row.append(';'.join(go_functions_list))
		writer.writerow(row)
	
	def output_function_group(self, curs, writer, function_struc_dict, gene_no2gene_id, go_no2go_id, go_no2name, go_no2accuracy, go_no2accuracy_pair):
		"""
		03-02-05
			output one group of functions which are parent-child, regarded as one distinct prediction
		"""
		for (go_no, function_struc) in function_struc_dict.iteritems():
			#transform to character type
			p_value_list = map(repr, function_struc.p_value_list)
			mcl_id_list = map(repr, function_struc.cluster_array)
			context_list = []
			for (gene_no,frequency) in function_struc.context_dict.iteritems():
				context_list.append('%s/%d'%(gene_no2gene_id[gene_no], frequency))
			row = [go_no, go_no2go_id[go_no], go_no2name[go_no], function_struc.is_correct,\
				function_struc.is_correct_L1, function_struc.is_correct_lca, ';'.join(p_value_list),\
				len(mcl_id_list), ';'.join(mcl_id_list), go_no2accuracy[go_no].ratio, \
				go_no2accuracy_pair[go_no].ratio, ';'.join(context_list)]
				
			writer.writerow(row)
		#a blank line
		writer.writerow([])
	
	def stat_output(self, outf, known_gene_no2p_gene_id_src, unknown_gene_no2p_gene_id_src):
		"""
		03-09-05
			give an overview stats for distinct function group of each gene
		"""
		sys.stderr.write("Outputting stats ... ")
		#make some blank lines
		outf.write("\n\n")
		list_of_function_groups_of_known_genes = map(len, known_gene_no2p_gene_id_src.values())
		list_of_function_groups_of_unknown_genes = map(len, unknown_gene_no2p_gene_id_src.values())
		no_of_known_genes = len(known_gene_no2p_gene_id_src)
		known_genes_with_multiple_function_groups = sum(greater(list_of_function_groups_of_known_genes,1))
		no_of_unknown_genes = len(unknown_gene_no2p_gene_id_src)
		unknown_genes_with_multiple_function_groups = sum(greater(list_of_function_groups_of_unknown_genes,1))
		no_of_genes = no_of_known_genes +no_of_unknown_genes
		avg_functions_groups_per_known_gene = sum(list_of_function_groups_of_known_genes)/float(no_of_known_genes)
		avg_functions_groups_per_unknown_gene = sum(list_of_function_groups_of_unknown_genes)/float(no_of_unknown_genes)
		avg_functions_groups_per_gene = \
			sum(list_of_function_groups_of_known_genes + list_of_function_groups_of_unknown_genes)/float(no_of_genes)
		outf.write("Total predicted genes: %s.\n"%no_of_genes)
		outf.write("\tAverage number of function groups per gene: %s.\n"%avg_functions_groups_per_gene)
		outf.write("Total known genes: %s. %s of them with multiple function groups\n"%\
			(no_of_known_genes, known_genes_with_multiple_function_groups))
		outf.write("\tAverage number of function groups per known gene: %s.\n"%avg_functions_groups_per_known_gene)
		outf.write("Total unknown genes: %s. %s of them with multiple function groups\n"%\
			(no_of_unknown_genes, unknown_genes_with_multiple_function_groups))
		outf.write("\tAverage number of function groups per unknown gene: %s.\n"%avg_functions_groups_per_unknown_gene)

		sys.stderr.write("Done.\n")

	
	def run(self):
		"""
		03-02-05
			initial
		"""
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		curs.execute("begin")	#because of cursor usage
		
		self.data_fetch(curs, self.p_gene_table, self.gene_p_table, self.mcl_table)
		self.output(curs, sys.stdout, self.known_gene_no2p_gene_id_src, self.unknown_gene_no2p_gene_id_src, self.p_gene_id_src_map)
		self.stat_output(sys.stdout, self.known_gene_no2p_gene_id_src, self.unknown_gene_no2p_gene_id_src)
	

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)	
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:n:m:j:cru", ["help", "hostname=", \
			"dbname=", "schema=", "p_gene_table=", "gene_p_table=", \
			"mcl_table=", "judger_type=", "commit", "report", "debug"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	p_gene_table = None
	gene_p_table = None
	mcl_table = None
	judger_type = 0
	commit = 0
	report = 0
	debug = 0
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
		elif opt in ("-t", "--p_gene_table"):
			p_gene_table = arg
		elif opt in ("-n", "--gene_p_table"):
			gene_p_table = arg
		elif opt in ("-m", "--mcl_table"):
			mcl_table = arg
		elif opt in ("-j", "--judger_type"):
			judger_type = int(arg)
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-u", "--debug"):
			debug = 1
	if schema and p_gene_table and gene_p_table and mcl_table:
		instance = gene_p_translate(hostname, dbname, schema, p_gene_table, gene_p_table,\
			mcl_table, judger_type, commit, report, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
