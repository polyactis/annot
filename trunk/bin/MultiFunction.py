#!/usr/bin/env python
"""
Usage: MultiFunction.py -k SCHEMA -f FNAME -s SUPPORT [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	schema list separated by coma
	-f ..., --netmine_fname=...	netmine_fname list separated by coma
	
	-s ..., --support=...	the number of schemas in which one gene or one function appear
	-p ..., --prefix=...	the output filename prefix, MultiFunction(default)
	
	-r, --report	report flag
	-b, --debug debug flag
	-h, --help              show this help

Examples:
	MultiFunction.py -k mm_oxi_stress_7t1 -f fmos_7t1g1e3d40q20s200c50z0001c8 -s 4
	
Description:
	Program to examine the multi functionality betwee different schemas.
"""

import sys, os, psycopg, getopt, csv, math
from sets import Set
from codense.common import db_connect, get_gene_no2id, get_go_id2name, get_go_term_id2depth, get_go_id2term_id

def gene_no2go_id_set_from_gene_p_table(input_fname, hostname='zhoudb', dbname='graphdb', schema='sc_new_38'):
	"""
	07-06-05
		copied from misc.py
	"""
	from sets import Set
	gene_no2go_id_set = {}
	p_gene_table = "p_gene_%s_e5"%input_fname
	gene_p_table = "gene_p_%s_e5_p01"%input_fname
	from codense.common import db_connect
	import psycopg
	conn, curs = db_connect(hostname, dbname, schema)
	curs.execute("select p.gene_no,go.go_id from %s p, go, %s g where p.p_gene_id=g.p_gene_id and go.go_no=p.go_no"%\
		(p_gene_table, gene_p_table))
	rows = curs.fetchall()
	for row in rows:
		gene_no = row[0]
		go_id = row[1]
		if gene_no not in gene_no2go_id_set:
			gene_no2go_id_set[gene_no] = Set()
		gene_no2go_id_set[gene_no].add(go_id)
	return gene_no2go_id_set


class MultiFunction:
	"""
	07-06-05
	"""
	def __init__(self, hostname=None, dbname=None, schema=None, netmine_fname=None,\
		support=1, prefix='MultiFunction', report=0, debug=0):
		self.hostname = hostname
		self.dbname = dbname
		if schema:
			self.schema_list = schema.split(',')
		else:
			self.schema_list = []
		if netmine_fname:
			self.netmine_fname_list = netmine_fname.split(',')
		else:
			self.netmine_fname_list = []
		self.support = int(support)
		self.prefix = prefix
		self.report = int(report)
		self.debug = int(debug)
	
	def fetch_predicted_function(self, schema_list, netmine_fname_list):
		"""
		07-06-05
		"""
		sys.stderr.write("Getting predictions...")
		gene_no2go_id_set_list = []
		go_id_set_list = []
		for i in range(len(schema_list)):
			gene_no2go_id_set = gene_no2go_id_set_from_gene_p_table(netmine_fname_list[i], schema=schema_list[i])
			go_id_set = Set()
			for value in gene_no2go_id_set.values():
				go_id_set |= value
			gene_no2go_id_set_list.append(gene_no2go_id_set)
			go_id_set_list.append(go_id_set)
		sys.stderr.write("Done.\n")
		return (gene_no2go_id_set_list, go_id_set_list)
	
	def output(self, gene_no2go_id_set_list, go_id_set_list, support, prefix, gene_no2id, go_id2name, schema_list):
		"""
		07-06-05
		"""
		sys.stderr.write("Outputing...")
		
		#get the total set
		total_gene_no_set = Set()
		total_go_id_set = Set()
		for i in range(len(gene_no2go_id_set_list)):
			total_gene_no_set |= Set(gene_no2go_id_set_list[i].keys())
			total_go_id_set |= go_id_set_list[i]
		
		gene_ofname = '%s.gene'%prefix
		function_ofname = '%s.function'%prefix
		gene_writer = csv.writer(open(gene_ofname,'w'), delimiter='\t')
		function_writer = csv.writer(open(function_ofname, 'w'), delimiter='\t')
		gene_writer.writerow(['']+schema_list)
		function_writer.writerow([''] + schema_list)
		
		from gene_p_map_redundancy import gene_p_map_redundancy
		node_distance_class = gene_p_map_redundancy()
		
		#output the gene-oriented information
		for gene_no in total_gene_no_set:
			freq = 0
			for i in range(len(gene_no2go_id_set_list)):
				if gene_no in gene_no2go_id_set_list[i]:
					freq += 1
			if freq == support:
				row = [gene_no2id[gene_no]]
				for i in range(len(gene_no2go_id_set_list)):
					if gene_no in gene_no2go_id_set_list[i]:
						row.append(';'.join(gene_no2go_id_set_list[i][gene_no]))
					else:
						row.append('')
				gene_writer.writerow(row)
		
		#output the function_oriented information
		for go_id in total_go_id_set:
			freq = 0
			for i in range(len(go_id_set_list)):
				if go_id in go_id_set_list[i]:
					freq += 1
			if freq == support:
				row = ['%s(%s)'%(go_id2name[go_id],go_id)]
				for i in range(len(go_id_set_list)):
					if go_id in go_id_set_list[i]:
						row.append('1')
					else:
						row.append('0')
				function_writer.writerow(row)
		
		
		sys.stderr.write("Done.\n")
	
	def run(self):
		"""
		07-07-05
		"""
		conn, curs  = db_connect(hostname, dbname)
		gene_no2id = get_gene_no2id(curs)
		go_id2name =get_go_id2name(curs)
		gene_no2go_id_set_list, go_id_set_list = self.fetch_predicted_function(self.schema_list,self.netmine_fname_list)
		self.output(gene_no2go_id_set_list,go_id_set_list, support, prefix, gene_no2id, go_id2name, self.schema_list)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:f:t:s:m:l:a:p:j:b:ocru", ["help", "hostname=", \
			"dbname=", "schema=", "netmine_fname=","support=", "prefix=", \
			"report", "debug"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	netmine_fname = None
	support = 1
	prefix = "MultiFunction"
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
		elif opt in ("-f", "--netmine_fname"):
			netmine_fname = arg
		elif opt in ("-s", "--support"):
			support = int(arg)
		elif opt in ("-p", "--prefix"):
			prefix = arg
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-u", "--debug"):
			debug = 1

	if schema and netmine_fname:
		instance = MultiFunction(hostname, dbname, schema, netmine_fname, support, prefix, report, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
