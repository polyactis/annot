#!/usr/bin/env python
"""
Usage: codense2db.py -k SCHEMA -p MAPPING_FILE [OPTION] DATAFILE

Option:
	DATAFILE usually is a components file, output by codense.
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	the splat_result table (edge_set)
	-m ..., --mcl_table=...	the mcl_result table (vertex_set)
	-p ..., --mapping_file=...	the file to get the mapping between haiyan's index and my gene_no
		(02-19-05)if not given, regard the input graph uses gene_no instead of haiyan's index
	-o ..., --cor_cut_off=...	the cor_cut_off for an edge to be valid, 0(default)
		NOTICE: 0 means the binary conversion won't be used, just summing the floats.
	-y ..., --parser_type=...	the type of parser to use, 1(copath, default), 2(codense)
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
		self.splat_connectivity = None	#03-03-05 this is the connectivity of the summary subgraph
		self.connectivity = None	#03-03-05 this is by averaging connectivities of the summary subgraph in each individual dataset

		
class codense2db:
	'''
	run()
		--init
		--get_haiyan_no2gene_no
		--create_tables
		--parser_dict[parser_type]
		--db_submit
	'''
	def __init__(self, infname=None, hostname='zhoudb', dbname='graphdb', schema=None, table=None,\
			mcl_table=None, mapping_file=None, cor_cut_off=0,\
			parser_type=1, needcommit=0, report=0):
		"""
		02-25-05
			modify the interface of the class and 2 member functions(get_combined_cor_vector, parse_recurrence)
			to make it module-independent
		"""
		self.infname = infname
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.table = table
		self.mcl_table = mcl_table
		self.mapping_file = mapping_file
		self.cor_cut_off = float(cor_cut_off)
		self.parser_type = int(parser_type)
		self.needcommit = int(needcommit)
		self.report = int(report)
	
		self.debug = 0
		self.parser_dict = {1:self.copath_parser,
			2: self.codense_parser}

	
	def copath_parser(self, row, argument=None, argument2=None):
		"""
		03-03-05
			get recurrence_array by passing the combined_sig_vetor
			get a fair recurrence_array(see 'redefine recurrence and connectivity' in log_05)
			
			the connectivity of the summary subgraph is passed to splat_connectivity
			
			--get_combined_cor_vector
			--parse_recurrence
			--parse_connectivity
		"""
		haiyan_no2gene_no = argument
		curs = argument2
		cluster = cluster_dstructure()
		cluster.cluster_id = int(row[0])
		cluster.splat_connectivity = float(row[1])
		cluster.vertex_set = row[2][1:-2].split(';')
		cluster.vertex_set = map(int, cluster.vertex_set)
		if self.mapping_file != None:
			cluster.vertex_set = dict_map(haiyan_no2gene_no, cluster.vertex_set)
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
			if self.mapping_file!=None:
				edge = dict_map(haiyan_no2gene_no, edge)
			#in ascending order
			edge.sort()
			cluster.edge_set.append(edge)
		cluster.no_of_edges = len(cluster.edge_set)
		(combined_cor_vector, combined_sig_vetor) = self.get_combined_cor_vector(curs, cluster.edge_set)
		cluster.connectivity = self.parse_connectivity(combined_sig_vetor, cluster.no_of_edges, len(cluster.vertex_set))
		
		cluster.recurrence_array = self.parse_recurrence(combined_sig_vetor, cluster.no_of_edges, self.cor_cut_off)
		
		return cluster
	
	def get_combined_cor_vector(self, curs, edge_set):
		"""
		03-03-05
			return (combined_cor_vector, combined_sig_vetor)
		"""
		combined_cor_vector = []
		combined_sig_vetor = []
		for edge in edge_set:
			edge_string = '{' + repr(edge)[1:-1] + '}'
			curs.execute("select cor_vector, sig_vector from edge_cor_vector where edge_name='%s'"%edge_string)
			rows = curs.fetchall()
			if len(rows) == 0:
				sys.stderr.write('%s not found in edge_cor_vector\n'%edge_string)
				sys.exit(1)
			cor_vector = rows[0][0][1:-1].split(',')
			cor_vector = map(float, cor_vector)
			combined_cor_vector += cor_vector
			
			sig_vector = rows[0][1][1:-1].split(',')
			sig_vector = map(int, sig_vector)
			combined_sig_vetor += sig_vector			
		return (combined_cor_vector, combined_sig_vetor)

	def parse_recurrence(self, combined_vector, no_of_edges, cor_cut_off=0):
		"""
		03-03-05
			replace the combined_vector with combined_sig_vetor
			get a fair recurrence_array(see 'redefine recurrence and connectivity' in log_05)
		"""
		cor_array = numarray.array(combined_vector)
		y_dimension = len(cor_array)/no_of_edges
		cor_array = numarray.reshape(cor_array, (no_of_edges, y_dimension))
		recurrence_array = []
		for i in range(y_dimension):
			#regard the correlations >= self.cor_cut_off to be 1, others 0
			if cor_cut_off>0:
				edge_cor_in_one_dataset = numarray.greater_equal(cor_array[:,i], cor_cut_off)
			else:
				edge_cor_in_one_dataset = cor_array[:,i]
			recurrence_array.append(sum(edge_cor_in_one_dataset)/float(no_of_edges))
			if self.debug:
				print edge_cor_in_one_dataset
				print recurrence_array
		#handle this in gene_stat_plot.py
		#recurrence_array = numarray.greater_equal(recurrence_array, self.subgraph_cut_off)
		return recurrence_array
	
	def parse_connectivity(self, combined_vector, no_of_edges, no_of_nodes):
		"""
		03-03-05
			parsing the combined_sig_vetor and return connectivity
		"""
		sig_array = numarray.array(combined_vector)
		y_dimension = len(sig_array)/no_of_edges
		sig_array = numarray.reshape(sig_array, (no_of_edges, y_dimension))
		connectivity_list = []
		for i in range(y_dimension):
			no_of_edges_in_one_dataset = sum(sig_array[:,i])
			if self.debug:
				print sig_array[:,i]
				
			connectivity_list.append(2.0*no_of_edges_in_one_dataset/(no_of_nodes*(no_of_nodes-1)) )
		if self.debug:
			print connectivity_list
		return sum(connectivity_list)/y_dimension
		
	
	def codense_parser(self, row, argument=None, argument2=None):
		"""
		03-04-05
			parser for codense results
		03-18-05
			previous codense results were processed by haiyan's
			program. Modify to parse real codense results.
		"""
		gene_id2gene_no = argument
		cluster = cluster_dstructure()
		cluster.cluster_id = int(row[0])
		no_of_nodes = int(row[1])
		no_of_edges = int(row[2])
		cluster.splat_connectivity = 2*float(no_of_edges)/(no_of_nodes*(no_of_nodes-1))	#row[3] should be equal to the below one got from no_of_edges and no_of_nodes
		cluster.connectivity = 2*float(no_of_edges)/(no_of_nodes*(no_of_nodes-1))
		cluster.no_of_edges = no_of_edges
		cluster.vertex_set = []
		for i in range(3,len(row)):
			if row[i]:
				#to avoid the last blank field
				cluster.vertex_set.append(int(row[i]))
		cluster.vertex_set = dict_map(gene_id2gene_no, cluster.vertex_set)
		cluster.edge_set = [[1,1]]	#dummy edge
		cluster.recurrence_array = [1]	#dummy recurrence
		
		return cluster
	
	def create_tables(self, curs, table, mcl_table):
		#create tables if necessary
		if table != 'splat_result':
			try:
				curs.execute("create table %s(\
					splat_id		serial primary key,\
					no_of_edges	integer,\
					recurrence_pattern	bit varying(200),\
					recurrence_array	float[],\
					edge_set	integer[][],\
					connectivity	float)"%table)
			except:
				sys.stderr.write("Error occurred when creating table %s\n"%table)
		if mcl_table != 'mcl_result':
			try:
				curs.execute("create table %s(\
					mcl_id	serial primary key,\
					splat_id	integer,\
					vertex_set	integer[],\
					parameter	varchar,\
					connectivity	float,\
					p_value_min	float,\
					go_no_vector	integer[],\
					unknown_gene_ratio	float,\
					recurrence_array	float[])"%mcl_table)
			except:
				sys.stderr.write("Error occurred when creating table %s\n"%mcl_table)	

	def db_submit(self, curs, cluster, table, mcl_table):
		"""
		03-03-05
			splat table's connectivity is the splat_connectivity, see doc above
		"""
		#try:
		#inserting into the splat_table
		curs.execute("insert into %s(splat_id, no_of_edges, \
					recurrence_array, edge_set, connectivity) values (%d, %d, ARRAY%s, ARRAY%s, %f)"%\
					(table, cluster.cluster_id, cluster.no_of_edges, repr(cluster.recurrence_array),\
					repr(cluster.edge_set), cluster.splat_connectivity))
		#inserting into the mcl_table
		curs.execute("insert into %s(mcl_id, splat_id, vertex_set, connectivity, recurrence_array)\
						values (%d, %d, ARRAY%s, %f, ARRAY%s)"%\
						(mcl_table, cluster.cluster_id, cluster.cluster_id, repr(cluster.vertex_set),\
						cluster.connectivity, repr(cluster.recurrence_array)) )
		"""
		except:
			sys.stderr.write('Error occurred when inserting pattern. Aborted.\n')
			self.conn.rollback()
			sys.exit(1)
		"""

	def run(self):
		"""
		03-18-05
			mapping_dict all changed to haiyan_no2gene_no
			
		"""

		inf = csv.reader(open(self.infname, 'r'), delimiter='\t')
		(conn, curs) = db_connect(self.hostname, self.dbname, self.schema)
		
		#setup the haiyan_no2gene_no
		if self.mapping_file != None:
			haiyan_no2gene_no = get_haiyan_no2gene_no(self.mapping_file)
		else:
			haiyan_no2gene_no = {}	#a blank dictionary, 
		gene_id2gene_no = get_gene_id2gene_no(curs)
		mapping_dict = {1:haiyan_no2gene_no,
			2:haiyan_no2gene_no}
		self.create_tables(curs, self.table, self.mcl_table)
		no = 0
		for row in inf:
			cluster = self.parser_dict[self.parser_type](row, mapping_dict[self.parser_type], curs)
			self.db_submit(curs, cluster, self.table, self.mcl_table)
			no+=1
			if self.report and no%1000==0:
				sys.stderr.write('%s%d'%('\x08'*20, no))
		if self.report:
			sys.stderr.write('%s%d'%('\x08'*20, no))
		if self.needcommit:
			conn.commit()
		sys.stderr.write('\n\tTotal patterns: %d\n'%no)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:m:p:o:y:cr", ["help", "hostname=", \
			"dbname=", "schema=", "table=", "mcl_table=", "mapping_file=", "cor_cut_off=",\
			"parser_type=", "commit", "report"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = 'splat_result'
	mcl_table = 'mcl_result'
	mapping_file = None
	cor_cut_off = 0
	parser_type = 1
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
		elif opt in ("-y", "--parser_type"):
			parser_type = int(arg)
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-r", "--report"):
			report = 1
	if schema and len(args)==1:
		instance = codense2db(args[0], hostname, dbname, schema, table, mcl_table, mapping_file, cor_cut_off,\
			parser_type, commit, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
