#!/usr/bin/env python
"""
Usage:	go_informative_node.py -a GO_ASSOCIATION -i GO_INDEX [OPTIONS}
	go_informative_node.py -k SCHEMA -b [OPTIONS]

Option:
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-a ..., -go_association=...	the file mapping go_id to gene_id
	-i ..., --go_index=...	the file mapping go_id to go_index, OR the go_graph file
	-s ..., --size=...	the size of the informative node, 60(default)
	-t, --stat	output in stats format
	-b, --bfs	construct by BFS instead of index
	-h, --help      show this help
	
Examples:
	go_informative_node.py -i process.pathindex.annotation.txt -a go2ug.txt
	go_informative_node.py -i term2term.txt -a go2ug.txt -b
	go_informative_node.py -i term2term.txt -a go2ug.txt -b -t
	go_informative_node.py -k hs_yh60_42 -b

Description:
	First usage:
	this program constructs informative nodes based on go_association file
	and go_index file or go_graph file.
	go_association file format:	go_id	gene_id
	go_index file format:	go_id	go_index
	go_graph file format:	term1_id	term2_id	go_id1	go_id2
	Default Output is stdout.
	go_association file should not have biological_process unknown class.
	
	Second usage:
	this program constructs go_association and go_graph directly from database.
	It depends on schema.gene, graph.association, go.term and go.term2term.
"""

import sys, os, psycopg, getopt, csv
from kjbuckets import *

class go_informative_node:
	def __init__(self, go_association, go_index, size, stat):
		self.goa_f = csv.reader(open(go_association, 'r'), delimiter='\t')
		self.goi_f = csv.reader(open(go_index, 'r'), delimiter='\t')
		self.size = int(size)
		self.stat = int(stat)
		self.go_index2gene_id_dict = {}
		#value utilizes kjSet data structure
		self.go_id2go_index_dict = {}
		self.go_index2go_id_dict = {}
		self.go_id2gene_id_dict = {}
		#value utilizes kjSet data structure
		self.informative_node_dict = {}
		self.log_file = open('/tmp/go_informative_node.log', 'w')
		
	def dstruc_loadin(self):
		for row in self.goa_f:
		#setup the go_id2gene_id_dict structure
			go_id = row[0]
			gene_id = row[1]
			if go_id not in self.go_id2gene_id_dict:
				self.go_id2gene_id_dict[go_id] = kjSet([gene_id])
			else:
				self.go_id2gene_id_dict[go_id].add(gene_id)
		
		for row in self.goi_f:
		#setup the go_id2go_index_dict structure
			go_id = 'GO:' + row[0]
			go_index = row[1]
			if go_index not in self.go_index2go_id_dict:
				self.go_index2go_id_dict[go_index] = go_id
			else:
				self.log_file.write("%s points to >1 go_ids\n"%go_index)
			if go_id not in self.go_id2go_index_dict:
				self.go_id2go_index_dict[go_id] = [go_index]
			else:
				self.go_id2go_index_dict[go_id].append(go_index)
				
		for go_id in self.go_id2go_index_dict:
		#setup the go_index2gene_id_dict structure
			go_index_list = self.go_id2go_index_dict[go_id]
			for go_index in go_index_list:
				if go_index not in self.go_index2gene_id_dict:
					if go_id in self.go_id2gene_id_dict:
						self.go_index2gene_id_dict[go_index] = self.go_id2gene_id_dict[go_id]
					else:
						self.go_index2gene_id_dict[go_index] = kjSet()
				else:
					self.log_file.write("Error: %s has >1 corresponding go_id\n"%go_index)
	
	def index_parent(self, go_index):
		index_list = go_index.split(',')
		index_list.pop()
		if index_list:
			return ','.join(index_list)
		else:
			return None
			
	def run(self):
		self.dstruc_loadin()
		for go_index in self.go_index2gene_id_dict:
		#collect the associated genes from descendent terms
			go_index_p = self.index_parent(go_index)
			while go_index_p:
				if go_index_p in self.go_index2gene_id_dict:
					self.go_index2gene_id_dict[go_index_p] += self.go_index2gene_id_dict[go_index]
				else:
					self.log_file.write("""%s doesn't exist\n"""%go_index_p)
				go_index_p = self.index_parent(go_index_p)
		
		go_index_list = self.go_index2gene_id_dict.keys()
		go_index_list.sort()
		for go_index in go_index_list:
		#find the informative nodes
			if len(self.go_index2gene_id_dict[go_index]) >= self.size:
				go_index_p = self.index_parent(go_index)
				while go_index_p:
					if go_index_p in self.informative_node_dict:
						del self.informative_node_dict[go_index_p]
					go_index_p = self.index_parent(go_index_p)
				self.informative_node_dict[go_index] = 1
				
		go_index_list = self.informative_node_dict.keys()
		self.informative_node_dict = {}
		for go_index in go_index_list:
		#transform go_index to go_id, cause go_index is redundant.
			self.log_file.write('informative node: %s\n'%go_index)
			go_id = self.go_index2go_id_dict[go_index]
			self.informative_node_dict[go_id] = 1
			
		self.output()
	
	def output(self):
		for go_id in self.informative_node_dict:
			go_index = self.go_id2go_index_dict[go_id][0]
			#take the first go_index associated with the go_id
			gene_id_list = self.go_index2gene_id_dict[go_index].items()
			if self.stat:
				sys.stdout.write("%s\t%d\n"%(go_id, len(gene_id_list)))
			else:
				for gene_id in gene_id_list:
					sys.stdout.write('%s\t%s\t%s\n'%(go_id, go_index, gene_id))


class go_informative_node_bfs:
	def __init__(self, dbname, schema, go_association, go_graph, size, stat):
		self.schema = schema
		if self.schema:
		#input from database
			self.conn = psycopg.connect('dbname=%s'%dbname)
			self.curs = self.conn.cursor()
			self.curs.execute("set search_path to %s"%self.schema)
		else:
		#input from files
			self.goa_f = csv.reader(open(go_association, 'r'), delimiter='\t')
			self.gog_f = csv.reader(open(go_graph, 'r'), delimiter='\t')
		self.size = int(size)
		self.stat = int(stat)
		self.root = 'GO:0008150'
		self.go_graph = kjGraph()
		self.go_id2gene_id_dict = {}
		#value utilizes kjSet data structure
		self.go_id_descendent2gene_id_dict = {}
		self.informative_node_dict = {}
		self.log_file = open('/tmp/go_informative_node.log', 'w')

	def dstruc_loadin(self):
		sys.stderr.write("Loading Data STructure...")
		for row in self.goa_f:
		#setup the go_id2gene_id_dict structure
			go_id = row[0]
			gene_id = row[1]
			if go_id not in self.go_id2gene_id_dict:
				self.go_id2gene_id_dict[go_id] = kjSet([gene_id])
			else:
				self.go_id2gene_id_dict[go_id].add(gene_id)
		
		for row in self.gog_f:
		#setup the go_graph structure
			self.go_graph.add((row[2], row[3]))
		self.go_id_set = kjSet(self.go_graph.keys() + self.go_graph.values())
		sys.stderr.write("Done\n")

	def dstruc_loadin_from_db(self):
		sys.stderr.write("Loading Data STructure...")
		#biological_process, not obsolete, not biological_process unknown.
		self.curs.execute("select a.go_id, g.gene_id from graph.association a, go.term t,\
			gene g where g.gene_id=a.gene_id and t.acc=a.go_id and\
			t.term_type='biological_process' and t.acc!='GO:0000004' and t.is_obsolete=0")
		rows = self.curs.fetchall()
		for row in rows:
		#setup the go_id2gene_id_dict structure
			go_id = row[0]
			gene_id = row[1]
			if go_id not in self.go_id2gene_id_dict:
				self.go_id2gene_id_dict[go_id] = kjSet([gene_id])
			else:
				self.go_id2gene_id_dict[go_id].add(gene_id)
		#get the non-obsolete biological_process GO DAG
		self.curs.execute("select t2t.term1_id, t2t.term2_id, t1.acc, t2.acc from \
			go.term2term t2t, go.term t1, go.term t2 where t2t.term1_id=t1.id and \
			t2t.term2_id=t2.id and t1.is_obsolete=0 and t2.is_obsolete=0 and \
			t1.term_type='biological_process' and t2.term_type='biological_process' ")
		rows = self.curs.fetchall()
		for row in rows:
		#setup the go_graph structure
			self.go_graph.add((row[2], row[3]))
		self.go_id_set = kjSet(self.go_graph.keys() + self.go_graph.values())
		sys.stderr.write("Done\n")

	def run(self):
		if self.schema:
		#input from database	
			self.dstruc_loadin_from_db()
		else:
		#input from files
			self.dstruc_loadin()
		for go_id in self.go_id_set.items():
			#collect the associated genes from descendent terms
			self.informative_node_dict[go_id] = 0
			descendent_set = self.go_graph.reachable(go_id)
			if go_id in self.go_id2gene_id_dict:
				self.go_id_descendent2gene_id_dict[go_id] = self.go_id2gene_id_dict[go_id]
			else:
				self.go_id_descendent2gene_id_dict[go_id] = kjSet()
			for descendent in descendent_set.items():
				if descendent in self.go_id2gene_id_dict:
					self.go_id_descendent2gene_id_dict[go_id] += self.go_id2gene_id_dict[descendent]
		#below is analogous to BFS,
		#informative_node_dict's role is similar to vertex coloring in real BFS
		self.informative_node_dict[self.root] = 1
		list_queue = [self.root]
		while list_queue:
			go_id = list_queue.pop(0)
			neighbor_list = self.go_graph.neighbors(go_id)
			for neighbor in neighbor_list:
				if len(self.go_id_descendent2gene_id_dict[neighbor]) >= self.size:
					self.log_file.write('%s ousted by %s\n'%(go_id,neighbor))
					self.informative_node_dict[go_id] = 0
					if self.informative_node_dict[neighbor] == 0:
					#trick is here. informative_node candidates and first touch
						list_queue.append(neighbor)
					self.informative_node_dict[neighbor] = 1
		self.output()
	
	def output(self):
		for go_id,value in self.informative_node_dict.iteritems():
			if value == 1:	
				gene_id_list = self.go_id_descendent2gene_id_dict[go_id].items()
				if self.stat:
					sys.stdout.write("%s\t%d\n"%(go_id, len(gene_id_list)))
				else:
					for gene_id in gene_id_list:
						sys.stdout.write('%s\t%s\t%s\n'%(go_id, go_id, gene_id))
				
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hd:k:a:i:s:tb", ["help", "dbname=", "schema=", "go_association=", "go_index=", "size=", "stat", "bfs"])
	except:
		print __doc__
		sys.exit(2)
	
	dbname = 'graphdb'
	schema = ''	
	go_association = ''
	go_index = ''
	size = 60
	stat = 0
	bfs = 0
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-k", "--schema"):
			schema = arg
		elif opt in ("-a", "--go_association"):
			go_association = arg
		elif opt in ("-i", "--go_index"):
			go_index = arg
		elif opt in ("-s", "--size"):
			size = int(arg)
		elif opt in ("-t", "--stat"):
			stat = 1
		elif opt in ("-b", "--bfs"):
			bfs = 1
	if bfs==1 and schema:
		instance = go_informative_node_bfs(dbname, schema, go_association, go_index, size, stat)
		instance.run()
	elif bfs==1 and go_association and go_index:
		instance = go_informative_node_bfs(dbname, schema, go_association, go_index, size, stat)
		instance.run()		
	elif go_association and go_index:
		instance = go_informative_node(go_association, go_index, size, stat)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
