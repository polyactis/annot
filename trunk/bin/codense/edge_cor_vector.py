#!/usr/bin/env python
"""
Usage: edge_cor_vector.py  -k SCHEMA -e COMPONENT_FILE -p MAPPING_FILE [OPTION] DATADIR

Option:
	DATADIR is the directory of the datasets to be choosen
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	the edge_cor_vector table
	-e ..., --edge_file=...
	-p ..., --mapping_file=...	the file to get the mapping between haiyan's index and my gene_no
	-c, --commit	commit this database transaction
	-r, --report	report the progress(a number)
	-h, --help              show this help
	
Examples:
	edge_cor_vector.py -k sc_54 -e components.txt -p sc_54_gene_id2no -t edge_cor_vector_1 gph_result/sc_54

Description:

"""

import sys, pickle, os, psycopg, getopt, csv, re
from sets import Set
from common import *

class edge_cor_vector:
	'''

	'''
	def __init__(self, hostname, dbname, schema, table, edge_file, mapping_file, needcommit, report, dir):
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.table = table
		self.edge_file = csv.reader(open(edge_file, 'r'), delimiter='\t')
		self.mapping_file = mapping_file
		self.report = int(report)
		self.needcommit = int(needcommit)
		self.dir = dir
		self.files = os.listdir(self.dir)
		
		self.haiyan_no2gene_no = {}
		#to get the number of the dataset
		self.p_no = re.compile(r'\d+$')
		self.gene_id2no = {}
		self.edge2cor_vector = {}
		
		
	def dstruc_loadin(self):
		'''
		This method loads in the data structures from database.
		'''
		sys.stderr.write("Loading Data STructure...")
		
		#setup self.gene_no2gene_id
		self.curs.execute("select gene_id, gene_no from gene")
		rows = self.curs.fetchall()
		for row in rows:
			self.gene_id2no[row[0]] = row[1]
		
		#setup the self.haiyan_no2gene_no
		self.haiyan_no2gene_no = get_haiyan_no2gene_no(self.mapping_file)

		#setup the keys in self.edge2cor_vector
		for row in self.edge_file:
			edge_list = row[3][2:-4].split(' );(')
			for edge in edge_list:
				edge = edge.split(',')
				edge = map(int, edge)
				edge = dict_map(self.haiyan_no2gene_no, edge)
				#in ascending order
				if edge[0] < edge[1]:
					self.edge2cor_vector[(edge[0], edge[1])] = [0]*len(self.files)
				else:
					self.edge2cor_vector[(edge[1], edge[0])] = [0]*len(self.files)

		sys.stderr.write("Done\n")

	def db_submit(self):
		sys.stderr.write("Database transacting...")
		
		#create tables if necessary
		if self.table != 'edge_cor_vector':
			try:
				self.curs.execute("create table %s(\
					edge_id serial	primary key,\
					edge_name	integer[],\
					cor_vector	float[])"%self.table)
			except:
				sys.stderr.write("Error occurred when creating table %s\n"%self.table)
		
		for edge in self.edge2cor_vector:
			sys.stdout.write('%s\t%s\n'%(edge, self.edge2cor_vector[edge]))
		
		if self.needcommit:				
			self.curs.execute("end")
		sys.stderr.write("done.\n")

	
	def run(self):
		self.dstruc_loadin()
		
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(self.files))
		for f in self.files:
			sys.stderr.write("%d/%d:\t%s"%(self.files.index(f)+1,len(self.files),f))
			
			no = int(self.p_no.search(self.infname).group())	#integer conversion
			f_path = os.path.join(self.dir, f)
			reader = csv.reader(file(f_path), delimiter='\t')
			for row in reader:
				if row[0] == 'e':
					gene_one = self.gene_id2no[row[1]]
					gene_two = self.gene_id2no[row[2]]
					correlation = float(row[3])
					if gene_one < gene_two:
						self.edge2cor_vector[(gene_one, gene_two)][no-1] = correlation
					else:
						self.edge2cor_vector[(gene_two, gene_one)][no-1] = correlation
			#output in a readable format
			del reader
		
		self.db_submit()
	
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:e:p:rc", ["help", "hostname=", "dbname=", "schema=",\
			"table=", "edge_file=", "mapping_file=", "report", "commit"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = 'edge_cor_vector'
	edge_file = None
	mapping_file = None
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
		elif opt in ("-e", "--edge_file"):
			edge_file = arg
		elif opt in ("-p", "--mapping_file"):
			mapping_file = arg
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-r", "--report"):
			report = 1	
	if schema and edge_file and len(args) == 1:
		instance = edge_cor_vector(hostname, dbname, schema, table, edge_file, mapping_file, commit, report, args[0])
		instance.run()
	else:
		print __doc__
		sys.exit(2)
