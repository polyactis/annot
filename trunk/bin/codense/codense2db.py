#!/usr/bin/env python
"""
Usage: codense2db.py -k SCHEMA -p MAPPING_FILE [OPTION] DATAFILE

Option:
	DATAFILE usually is a components file, output by codense.
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	the splat_result table (edge_set)
	-m ..., -mcl_table=...	the mcl_result table (vertex_set)
	-p ..., --mapping_file=...	the file to get the mapping between haiyan's index and my gene_no
	-o ..., --cor_cut_off=...	the cor_cut_off for an edge to be valid, 0.6(default)
	-c, --commit	commit this database transaction
	-r, --report	report the progress(a number)
	-h, --help	show this help
	
Examples:
	codense2db.py -k sc_54 -t splat_result_1 -m mcl_result_1 -p sc_54_gene_id2no  -c -r components.txt
	
Description:
	Parse the codense results and import into schema.splat_result and mcl_result tables.
	This program must be run after the edge_cor_vector is setup.
	
	The recurrence_array's and recurrence_cut_off's in all the tables of schema
	are changed to be type float.

"""


import sys,os,psycopg,getopt,csv, numarray
from common import *

class cluster_dstructure:
	def __init__(self):
		self.cluster_id = None
		self.vertex_set = None
		#string form is for database submission
		self.vertex_set_string = None
		self.edge_set = None
		#string form is for database submission
		self.edge_set_string = None
		self.no_of_edges = None
		self.recurrence_array = None
		self.connectivity = None

		
class codense2db:
	'''
	'''
	def __init__(self, infname, hostname, dbname, schema, table, mcl_table, mapping_file, cor_cut_off,\
			report, needcommit=0):		
		self.inf = csv.reader(open(infname, 'r'), delimiter='\t')
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.table = table
		self.mcl_table = mcl_table
		self.mapping_file = mapping_file
		self.cor_cut_off = float(cor_cut_off)
		self.report = int(report)
		self.needcommit = int(needcommit)		
		self.haiyan_no2gene_no = {}
		
	def parse(self, row):
		cluster = cluster_dstructure()
		cluster.cluster_id = int(row[0])
		cluster.connectivity = float(row[1])
		cluster.vertex_set = row[2][1:-2].split(';')
		cluster.vertex_set = map(int, cluster.vertex_set)
		cluster.vertex_set = dict_map(self.haiyan_no2gene_no, cluster.vertex_set)
		#in ascending order
		cluster.vertex_set.sort()
		#cluster.vertex_set_string = row[2][1:-2].replace(';',',')
		#cluster.vertex_set_string = '{' + cluster.vertex_set + '}'
		#cluster.edge_set_string = row[3][2:-3].replace(' );(', '},{')
		#cluster.edge_set_string = '{{' + cluster.edge_set_string + '}}'
		edge_list = row[3][2:-4].split(' );(')
		cluster.edge_set = []
		for edge in edge_list:
			edge = edge.split(',')
			edge = map(int, edge)
			edge = dict_map(self.haiyan_no2gene_no, edge)
			#in ascending order
			edge.sort()
			cluster.edge_set.append(edge)
		cluster.no_of_edges = len(cluster.edge_set)
		combined_cor_vector = self.get_combined_cor_vector(cluster.edge_set)
		cluster.recurrence_array = self.parse_recurrence(combined_cor_vector, cluster.no_of_edges)
		
		return cluster
	
	def get_combined_cor_vector(self, edge_set):
		combined_cor_vector = []
		for edge in edge_set:
			edge_string = '{' + repr(edge)[1:-1] + '}'
			self.curs.execute("select cor_vector from edge_cor_vector where edge_name='%s'"%edge_string)
			rows = self.curs.fetchall()
			if len(rows) == 0:
				sys.stderr.write('%s not found in edge_cor_vector\n'%edge_string)
				sys.exit(1)
			cor_vector = rows[0][0][1:-1].split(',')
			cor_vector = map(float, cor_vector)
			combined_cor_vector += cor_vector
		return combined_cor_vector

	def parse_recurrence(self, combined_cor_vector, no_of_edges):
		cor_array = numarray.array(combined_cor_vector)
		y_dimension = len(cor_array)/no_of_edges
		cor_array = numarray.reshape(cor_array, (no_of_edges, y_dimension))
		recurrence_array = []
		for i in range(y_dimension):
			#regard the correlations >= self.cor_cut_off to be 1, others 0
			edge_cor_in_one_dataset = numarray.greater_equal(cor_array[:,i], self.cor_cut_off)
			recurrence_array.append(sum(edge_cor_in_one_dataset)/float(no_of_edges))
		#handle this in gene_stat_plot.py
		#recurrence_array = numarray.greater_equal(recurrence_array, self.subgraph_cut_off)
		return recurrence_array
	
	def create_tables(self):
		#create tables if necessary
		if self.table != 'splat_result':
			try:
				self.curs.execute("create table %s(\
					splat_id		serial primary key,\
					no_of_edges	integer,\
					recurrence_pattern	bit varying(200),\
					recurrence_array	float[],\
					edge_set	integer[][],\
					connectivity	float)"%self.table)
			except:
				sys.stderr.write("Error occurred when creating table %s\n"%self.table)
		if self.mcl_table != 'mcl_result':
			try:
				self.curs.execute("create table %s(\
					mcl_id	serial primary key,\
					splat_id	integer,\
					vertex_set	integer[],\
					parameter	varchar,\
					connectivity	float,\
					p_value_min	float,\
					go_no_vector	integer[],\
					unknown_gene_ratio	float,\
					recurrence_array	float[])"%self.mcl_table)
			except:
				sys.stderr.write("Error occurred when creating table %s\n"%self.mcl_table)	

	def db_submit(self, cluster):
		try:
			#inserting into the splat_table
			self.curs.execute("insert into %s(splat_id, no_of_edges, \
						recurrence_array, edge_set, connectivity) values (%d, %d, ARRAY%s, ARRAY%s, %f)"%\
						(self.table, cluster.cluster_id, cluster.no_of_edges, repr(cluster.recurrence_array),\
						repr(cluster.edge_set), cluster.connectivity))
			#inserting into the mcl_table
			self.curs.execute("insert into %s(mcl_id, splat_id, vertex_set, connectivity, recurrence_array)\
							values (%d, %d, ARRAY%s, %f, ARRAY%s)"%\
							(self.mcl_table, cluster.cluster_id, cluster.cluster_id, repr(cluster.vertex_set),\
							cluster.connectivity, repr(cluster.recurrence_array)) )
		except:
			sys.stderr.write('Error occurred when inserting pattern. Aborted.\n')
			self.conn.rollback()
			sys.exit(1)

	def run(self):
		#setup the self.haiyan_no2gene_no	
		self.haiyan_no2gene_no = get_haiyan_no2gene_no(self.mapping_file)
		
		self.create_tables()
		no = 0
		for row in self.inf:
			cluster = self.parse(row)
			self.db_submit(cluster)			
			no+=1
			if self.report and no%1000==0:
				sys.stderr.write('%s%d'%('\x08'*20, no))
		if self.report:
			sys.stderr.write('%s%d'%('\x08'*20, no))
		if self.needcommit:
			self.conn.commit()
		sys.stderr.write('\n\tTotal patterns: %d\n'%no)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:m:p:o:cr", ["help", "hostname=", \
			"dbname=", "schema=", "table=", "mcl_table=", "mapping_file=", "cor_cut_off=",\
			"commit", "report"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = 'splat_result'
	mcl_table = 'mcl_result'
	mapping_file = None
	cor_cut_off = 0.6
	commit = 0
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
		elif opt in ("-t", "--table"):
			table = arg
		elif opt in ("-m", "--mcl_table"):
			mcl_table = arg
		elif opt in ("-p", "--mapping_file"):
			mapping_file = arg
		elif opt in ("-o", "--cor_cut_off"):
			cor_cut_off = float(arg)
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-r", "--report"):
			report = 1
	if schema and len(args)==1:
		instance = codense2db(args[0], hostname, dbname, schema, table, mcl_table, mapping_file, cor_cut_off,\
			report, commit)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
