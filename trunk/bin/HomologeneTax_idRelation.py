#!/usr/bin/env python
"""
Usage: TFBindingSiteParse.py [OPTIONS] -k SCHEMA -i INPUT_FILE -y PARSER_TYPE

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database, homologene(default)
	-s ...,	source tax_id, 9606(default, human)
	-b,	debug mode
	-r,	report the progress
	-c,	commit this database transaction
	-h,	show this help

Examples:
		HomologeneTax_idRelation.py -c -b >/tmp/stdout

Description:
	Program to fill up homologene.tax_id_index and homologene.tax_id_relationship
	
"""

import sys, getopt, os, csv
sys.path += [os.path.join(os.path.expanduser('~/script/annot/bin'))]
sys.path += [os.path.join(os.path.expanduser('~/script/microarray/bin'))]
from codense.common import db_connect, org_short2long, org2tax_id
from sets import Set
from graphlib import Graph, GraphAlgo

class HomologeneTax_idRelation:
	def __init__(self, hostname='zhoudb', dbname='mdb', schema='', src_tax_id=9606, \
		debug=0, report=0, commit=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.src_tax_id = int(src_tax_id)
		self.debug = int(debug)
		self.report = int(report)
		self.commit = int(commit)
	
	def construct_tax_tree(self, curs, nodes_table='taxonomy.nodes'):
		sys.stderr.write("Constructing tax_tree...")
		tax_tree = Graph.Graph()
		curs.execute("DECLARE crs CURSOR FOR select tax_id, parent_tax_id from %s"%\
			(nodes_table))
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		while rows:
			for row in rows:
				if row[0] ==1 and row[1] == 1:	#remove this loop node
					continue
				tax_tree.add_edge(row[1], row[0])	#from parent to child
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		curs.execute("close crs")
		sys.stderr.write("Done.\n")
		return tax_tree
	
	def get_tax_id_set(self, curs, homologene_table='homologene.homologene'):
		sys.stderr.write("Getting tax_id_set...")
		curs.execute("select distinct tax_id from %s"%homologene_table)
		rows = curs.fetchall()
		tax_id_set = Set()
		for row in rows:
			tax_id_set.add(row[0])
		sys.stderr.write("Done.\n")
		return tax_id_set
	
	def get_tax_id2index_given_tax_id_set(self, tax_tree, tax_id_set):
		sys.stderr.write("Getting tax_id2index...")
		tax_id2index = {}
		for tax_id in tax_id_set:
			inc_nbrs_list = tax_tree.inc_nbrs(tax_id)
			index = [tax_id]	#itself is the end of the index
			while inc_nbrs_list:
				if len(inc_nbrs_list)>=2:
					sys.stderr.write("tax_id: %s, index: %s, inc_nbrs_list: %s, >1 parents.\n"%\
						(tax_id, index, inc_nbrs_list))
					sys.exit(3)
				if self.debug:
					print inc_nbrs_list
					print index
				index.insert(0, inc_nbrs_list[0])	#ancestor in time order
				inc_nbrs_list = tax_tree.inc_nbrs(inc_nbrs_list[0])
			tax_id2index[tax_id] = index
		sys.stderr.write("Done.\n")
		return tax_id2index
	
	def submit2tax_id_index(self, curs, tax_id2index, tax_id_index_table='homologene.tax_id_index'):
		sys.stderr.write('Submitting tax_id2index to tax_id_index_table...')
		for tax_id, index in tax_id2index.iteritems():
			curs.execute("insert into %s(tax_id, index, depth) values(%s, ARRAY%s, %s)"%\
				(tax_id_index_table, tax_id, repr(index), len(index)))
		sys.stderr.write("Done.\n")
	
	def common_ancestor_of_2indices(self, index1, index2):
		"""
		12-12-05 correct a bug when for loop is not breaked
		"""
		length = min(len(index1), len(index2))
		for i in range(length):
			if index1[i] != index2[i]:
				break
		if index1[i] == index2[i]:	#the for loop is not breaked
			i = i+1
		return index1[:i], i	#common ancestor and its depth
	
	def submit_common_ancestor(self, src_tax_id, tax_id2index, curs, \
		tax_id_relationship_table = 'homologene.tax_id_relationship'):
		
		sys.stderr.write("Submitting common ancestor...")
		tax_id_pair2common_ancestor = {}
		for tax_id, index in tax_id2index.iteritems():
			common_ancestor, depth = self.common_ancestor_of_2indices(tax_id2index[src_tax_id], index)
			curs. execute("insert into %s(src_tax_id, tg_tax_id, common_ancestor_index, common_ancestor_depth)\
				values(%s, %s, ARRAY%s, %s)"%(tax_id_relationship_table, src_tax_id, tax_id, repr(common_ancestor), depth))
		sys.stderr.write("Done.\n")
	
	def run(self):
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		tax_tree = self.construct_tax_tree(curs)
		tax_id_set = self.get_tax_id_set(curs)
		tax_id2index = self.get_tax_id2index_given_tax_id_set(tax_tree, tax_id_set)
		self.submit2tax_id_index(curs, tax_id2index)
		self.submit_common_ancestor(self.src_tax_id, tax_id2index, curs)
		if self.commit:
			curs.execute("end")

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "schema=", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:s:brch", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = 'homologene'
	src_tax_id = 9606
	debug = 0
	report = 0
	commit = 0
	
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
		elif opt in ("-s"):
			src_tax_id = int(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-c", "--commit"):
			commit = 1

	if hostname and dbname and schema and src_tax_id:
		instance = HomologeneTax_idRelation(hostname, dbname, schema, src_tax_id,\
			debug, report, commit)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
