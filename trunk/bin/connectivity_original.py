#!/usr/bin/env python
"""
Usage:	connectivity_original.py -k SCHEMA -t TABLE [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	a 'mcl_result or fim_result' -like table, compute which table's connectivity
	-e ..., --edge_table=...	the table containing all edge information. 'edge_cor_vector'(default)
	-w ..., --min_weight=...	the minimum weight of an edge. 6(default)
	-s ..., --dict_threshold=...	the threshold of the internal dictionary, 10,000,000(default) 26.5% of 6G
	-b, --debug	enable debugging, no debug by default
	-c, --commit	commit the database transaction
	-r, --report	report the progress(a number)
	-h, --help	show this help
	
Examples:
	connectivity_original.py -k sc_54_6661 -c -t mcl_result
		:compute the connectivity of mcl_result

Description:
	This program computes the original connectivity of a cluster, which is the connectivity
	of the vertex set of the cluster in the summary graph.

"""

import sys, os, psycopg, re, getopt
from codense.common import db_connect

class connectivity_original:
	"""
	04-01-05
	
	--run
		--db_connect
		--alter_table
		--data_fetch
			--_connectivity_original
				--get_edge
		--update_connectivity_original
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, \
		table=None, edge_table='edge_cor_vector', min_weight=6, \
		dict_threshold=10000000, debug=0, report=0, needcommit=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.table = table
		self.edge_table = edge_table
		self.min_weight = int(min_weight)
		self.dict_threshold = int(dict_threshold)
		self.debug = int(debug)
		self.report = int(report)
		self.needcommit = int(needcommit)
	
		#a dictionary limited by dict_threshold, key is (vertex1,vertex2), ascending order, value is its weight.
		self.edge_dict = {}
		#mcl_id v.s. connectivity_original mapping.
		self.mcl_id2connectivity = {}
		#counter
		self.no_of_records = 0
		
	def data_fetch(self, curs, table, edge_table='edge_cor_vector', min_weight=6):
		"""
		04-01-05
			_connectivity_original
		"""
		sys.stderr.write("Computing cluster's density in its source graph...\n")
		curs.execute("DECLARE crs CURSOR FOR select mcl_id, vertex_set \
			from %s"%(table))
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		while rows:
			for row in rows:
				self._connectivity_original(curs, row, edge_table, min_weight)
				self.no_of_records+=1

			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, self.no_of_records))
			
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		
		sys.stderr.write("Done\n")		
	
	def _connectivity_original(self, curs, row, edge_table, min_weight):
		"""
		04-01-05
			input: mcl_id and vertex_set(in row)
			output: update its connectivity_original
			
			--get_edge()
		"""
		mcl_id = row[0]
		if self.debug:
			print "mcl_id: %s"%mcl_id
		vertex_set = row[1][1:-1].split(',')
		vertex_set = map(int, vertex_set)
		#sort it first.
		vertex_set.sort()
		no_of_nodes = len(vertex_set)
		no_of_edges = 0
		for i in range(no_of_nodes):
			for j in range(i+1, no_of_nodes):
				edge_weight = self.get_edge(curs, (vertex_set[i], vertex_set[j]), edge_table, min_weight)
				if edge_weight:	#Not None, means such an edge with weight>=min_weight exists
					if self.debug:
						print "edge (%s, %s) weight: %s"%(vertex_set[i], vertex_set[j], edge_weight)
					no_of_edges+=1
		connectivity = float(2*no_of_edges)/(no_of_nodes*(no_of_nodes-1))
		if self.debug:
			print "connectivity: %s"%connectivity
			raw_input("WAIT:\t")
		self.mcl_id2connectivity[mcl_id] = connectivity
	
	def get_edge(self, curs, edge, edge_table='edge_cor_vector', min_weight=6):
		"""
		04-01-05
			input: two vertices of an edge
			output: return None if no such edge, its weight if it exists.
				
		"""
		if edge in self.edge_dict:
			#the internal dictionary has it
			weight = self.edge_dict[edge]
			if weight>=min_weight:
				if self.debug:
					print "edge %s in internal dictionary with weight %s"%(repr(edge), weight)
				return weight
			else:
				if self.debug:
					print "edge %s in internal dictionary but with weight %s<%s"%(repr(edge), weight,min_weight)
				return None
		else:		
			#check database
			weight = self._get_edge(curs, edge, edge_table)
			if weight:
				if len(self.edge_dict)>=self.dict_threshold:
					#over the threshold, throw away the first item.
					self.edge_dict.popitem()
				self.edge_dict[edge] = weight
				
				if weight>=min_weight:
					#bigger than the min_weight
					if self.debug:
						print "edge %s in database with weight %s"%(repr(edge), weight)
					return weight
				else:
					if self.debug:
						print "edge %s in database but with weight %s<%s"%(repr(edge), weight, min_weight)
					return None
			else:
				if self.debug:
					print "edge %s not found"%(repr(edge))
				return None
	
	def _get_edge(self, curs, edge, edge_table='edge_cor_vector'):
		"""
		04-03-05
			split from get_edge(), to be class independent
		"""
		curs.execute("select sig_vector from %s where edge_name='{%s,%s}'"%(edge_table, edge[0], edge[1]))
		rows = curs.fetchall()
		if len(rows)>0:
			#edge is present
			edge_data = rows[0][0][1:-1]
			edge_data = edge_data.split(',')
			edge_data = map(int, edge_data)
			weight = sum(edge_data)
		else:
			weight = None
		return weight
	
	def alter_table(self, curs, table):
		"""
		04-01-05
			add a column, connectivity_original to the table
		"""
		sys.stderr.write("Adding connectivity_original to %s..."%table)
		curs.execute("alter table %s add connectivity_original float"%table)
		sys.stderr.write("Done.\n")
		
	def update_connectivity_original(self, curs, table, mcl_id2connectivity):
		"""
		04-01-05
			update connectivity_original in table
		"""
		sys.stderr.write("Database updating %s..."%table)
		for mcl_id,connectivity in mcl_id2connectivity.iteritems():
			curs.execute("update %s set connectivity_original=%s where mcl_id=%d"%\
			(table, connectivity, mcl_id))
		sys.stderr.write("Done.\n")
	
	def run(self):
		"""
		04-01-05
		"""
		(conn, curs) = db_connect(self.hostname, self.dbname, self.schema)
		if self.needcommit:
			self.alter_table(curs, self.table)
		self.data_fetch(curs, self.table, self.edge_table, self.min_weight)
		if self.needcommit:
			self.update_connectivity_original(curs, self.table, self.mcl_id2connectivity)
			curs.execute("end")
		
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:e:w:s:brc", \
			["help", "hostname=", "dbname=", "schema=", "table=", "edge_table=", \
			"min_weight=", "dict_threshold=", "debug", "report", "commit"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = ''
	edge_table = 'edge_cor_vector'
	min_weight = 6
	dict_threshold = 10000000
	debug = 0
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
		elif opt in ("-e", "--edge_table"):
			edge_table = arg
		elif opt in ("-w", "--min_weight"):
			min_weight = int(arg)
		elif opt in ("-s", "--dict_threshold"):
			dict_threshold = int(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-r", "--report"):
			report = 1


	if schema and table:
		instance = connectivity_original(hostname, dbname, schema, table, \
			edge_table, min_weight, dict_threshold, debug, report, commit)
		instance.run()

	else:
		print __doc__
		sys.exit(2)
