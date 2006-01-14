#!/usr/bin/env python
"""
Usage: CoexprFromCooccu.py -k SCHEMA -m MCL_TABLE [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	the splat_result table (edge_set)
	-m ..., --mcl_table=...	the mcl_result table (vertex_set)
	-e ..., --e_mcl_table=...	the e_mcl_table stores the 2nd-order clusters of mcl_table
		if not given, it's $mcl_table+'e'
	-r, --report	report the progress(a number)
	-h, --help	show this help
	
Examples:
	CoexprFromCooccu.py -k sc_54_6661 -m mcl_result
	
Description:
	This program finds the co-expression modules which are from the same cooccurrence module.

"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path += [os.path.expanduser('~/lib64/python')]
else:   #32bit
	sys.path += [os.path.expanduser('~/lib/python')]
import psycopg, getopt, csv, numarray, re
from codense.common import db_connect
from sets import Set

class CoexprFromCooccu:
	"""
	04-11-05
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, \
		table=None, mcl_table='', e_mcl_table=None, p_value_cut_off=0.01, report=0, \
		judger_type=0, needcommit=0, gene_table='p_gene', lm_table=None, \
		stat_table_fname=None, debug=0):
		
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema		
		self.table = table
		self.mcl_table = mcl_table
		self.e_mcl_table = e_mcl_table
		if self.e_mcl_table==None:
			self.e_mcl_table = self.mcl_table+'e'
		self.p_value_cut_off = float(p_value_cut_off)
		self.report = int(report)
		self.judger_type = int(judger_type)
		self.needcommit = int(needcommit)
		self.gene_table = gene_table
		self.lm_table = lm_table
		self.stat_table_fname = stat_table_fname
		#debugging flag
		self.debug = int(debug)
		
		self.no_of_records = 0
		self.pre_2nd_cc_hierarchy = {}
		
	def data_fetch(self, curs, mcl_table, e_mcl_table):
		"""
		04-11-05
		
		04-19-05
		"""
		sys.stderr.write("Constructing the pre_2nd_cc_hierarchy from database...\n")
		curs.execute("DECLARE crs CURSOR FOR select m1.mcl_id, m2.mcl_id, m2.cooccurrent_cluster_id \
			from %s m1, %s m2 where m1.cooccurrent_cluster_id=m2.cooccurrent_cluster_id"%(mcl_table, e_mcl_table))
		
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		while rows:
			for row in rows:
				self._coexpr_from_cocccu(row)	
				self.no_of_records+=1

			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, self.no_of_records))
			
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		
		sys.stderr.write("Done\n")		
		return self.pre_2nd_cc_hierarchy	#return for the purpose of cluster_info.py
	
	def _coexpr_from_cocccu(self, row):
		"""
		04-11-05
		
		04-19-05
			
		"""
		mcl_id = row[0]
		mcl_id_2nd_order = row[1]
		cooccurrent_cluster_id = row[2]
		pregraph_id,hhu_2nd_order_id = map(int, cooccurrent_cluster_id.split('.'))
		if pregraph_id not in self.pre_2nd_cc_hierarchy:
			self.pre_2nd_cc_hierarchy[pregraph_id] = {}
		if mcl_id_2nd_order not in self.pre_2nd_cc_hierarchy[pregraph_id]:
			self.pre_2nd_cc_hierarchy[pregraph_id][mcl_id_2nd_order] = Set()
		self.pre_2nd_cc_hierarchy[pregraph_id][mcl_id_2nd_order].add(mcl_id)
	
	def output(self, pre_2nd_cc_hierarchy):
		"""
		04-11-05
		
		04-19-05
			
		"""
		writer = csv.writer(sys.stdout, delimiter='\t')
		for pregraph_id,mcl_id_2nd_order_dict in pre_2nd_cc_hierarchy.iteritems():
			writer.writerow([pregraph_id])
			for mcl_id_2nd_order,mcl_id_set in mcl_id_2nd_order_dict.iteritems():
				writer.writerow(['', mcl_id_2nd_order])	#one indentation
				for mcl_id in mcl_id_set:
					writer.writerow(['', '', mcl_id])	#two indentations
		del writer
		
	def run(self):
		"""
		04-11-05
		
		04-19-05
			generate pre_2nd_cc_hierarchy.
		"""
		(conn, curs) = db_connect(self.hostname, self.dbname, self.schema)
		self.data_fetch(curs, self.mcl_table, self.e_mcl_table)
		self.output(self.pre_2nd_cc_hierarchy)


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema=", "table=", "mcl_table=", "e_mcl_table=", "p_value_cut_off=",\
		"judger_type=", "report", "commit", "gene_table=", "lm_table=", "debug", "accuracy_cut_off=",\
		"gene_p_table=",  "recurrence_gap_size=", "connectivity_gap_size="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:m:e:p:j:rcg:l:ba:n:x:y:", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = None
	mcl_table = ''
	e_mcl_table = None
	p_value_cut_off = 0.001
	judger_type = 0
	report = 0
	commit = 0
	gene_table = 'p_gene'
	lm_table = None
	debug = 0
	accuracy_cut_off = 0
	stat_table_fname = None
	gene_p_table = None
	recurrence_gap_size = 2
	connectivity_gap_size = 2
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
		elif opt in ("-e", "--e_mcl_table"):
			e_mcl_table = arg
		elif opt in ("-p", "--p_value_cut_off"):
			p_value_cut_off = float(arg)
		elif opt in ("-j", "--judger_type"):
			judger_type = int(arg)
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-g", "--gene_table"):
			gene_table = arg
		elif opt in ("-l", "--lm_table"):
			lm_table = arg
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-a", "--accuracy_cut_off="):
			accuracy_cut_off = float(arg)
		elif opt in ("-n", "--gene_p_table="):
			gene_p_table = arg
		elif opt in ("-x", "--recurrence_gap_size"):
			recurrence_gap_size = int(arg)
		elif opt in ("-y", "--connectivity_gap_size"):
			connectivity_gap_size = int(arg)

			
	if schema and mcl_table:
		instance = CoexprFromCooccu(hostname, dbname, schema, table, mcl_table, e_mcl_table, p_value_cut_off,\
			report, judger_type, commit, gene_table, lm_table, stat_table_fname, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
