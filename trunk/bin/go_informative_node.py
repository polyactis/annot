#!/usr/bin/env python
"""
Usage: go_informative_node.py -a GO_ASSOCIATION -i GO_INDEX

Option:
	-a ..., -go_association=...	the file mapping go_id to gene_id
	-i ..., --go_index=...	the file mapping go_id to go_index, OR the go_graph file
	-s ..., --size=...	the size of the informative node, 60(default)
	-b, --bfs	construct by BFS instead of index
	-h, --help      show this help
	
Examples:
	go_index.py -i process.pathindex.annotation.txt -a go2ug.txt
	go_index.py -i term2term.txt -a go2ug.txt -b
Description:
	this program constructs informative nodes based on go_association file
	and go_index file.
	go_association file format:	go_id	gene_id
	go_index file format:	go_id	go_index
	go_graph file format:	term_id1	term_id2	go_id1	go_id2
	Default Output is stdout.
"""

import sys, os, getopt, csv
from kjbuckets import *

class go_informative_node:
	def __init__(self, go_association, go_index, size):
		self.goa_f = csv.reader(open(go_association, 'r'), delimiter='\t')
		self.goi_f = csv.reader(open(go_index, 'r'), delimiter='\t')
		self.size = int(size)
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
			for gene_id in gene_id_list:
				sys.stdout.write('%s\t%s\t%s\n'%(go_id, go_index, gene_id))


class go_informative_node_bfs:
	def __init__(self, go_association, go_graph, size):
		self.goa_f = csv.reader(open(go_association, 'r'), delimiter='\t')
		self.gog_f = csv.reader(open(go_graph, 'r'), delimiter='\t')
		self.size = int(size)
		self.root = 'GO:0008150'
		self.go_graph = kjGraph()
		self.go_id2gene_id_dict = {}
		#value utilizes kjSet data structure
		self.go_id_descendent2gene_id_dict = {}
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
		
		for row in self.gog_f:
		#setup the go_graph structure
			self.go_graph.add((row[2], row[3]))
		self.go_id_set = kjSet(self.go_graph.keys() + self.go_graph.values())
		
	def run(self):
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
				for gene_id in gene_id_list:
					sys.stdout.write('%s\t%s\t%s\n'%(go_id, go_id, gene_id))
				
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "ha:i:s:b", ["help", "go_association=", "go_index=", "size=", "bfs"])
	except:
		print __doc__
		sys.exit(2)
	
	go_association = ''
	go_index = ''
	size = 60
	bfs = 0
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-a", "--go_association"):
			go_association = arg
		elif opt in ("-i", "--go_index"):
			go_index = arg
		elif opt in ("-s", "--size"):
			size = int(arg)
		elif opt in ("-b", "--bfs"):
			bfs = 1
	if bfs==1 and go_association and go_index:
		instance = go_informative_node_bfs(go_association, go_index, size)
		instance.run()		
	elif go_association and go_index:
		instance = go_informative_node(go_association, go_index, size)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
