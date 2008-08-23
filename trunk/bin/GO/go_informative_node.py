#!/usr/bin/env python
"""

Examples:
	go_informative_node.py -i process.pathindex.annotation.txt -a go2ug.txt
	go_informative_node.py -i term2term.txt -a go2ug.txt -b
	go_informative_node.py -i term2term.txt -a go2ug.txt -b -t 1
	go_informative_node.py -k hs_yh60_42 -b
	
	#2006 output informative nodes
	go_informative_node.py -k sc_38_all_no_info -n 2 -b
	
	#2008-08-22 generate informative nodes and output into file
	go_informative_node.py -d go -o /tmp/go_informative_node.tsv
	
	go_informative_node.py -d go -o /tmp/go_informative_node.tsv   -n 2 -s 150
	
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
	
	12-30-05 A snippet is commented out which could be used to do
		species-specific association tasks.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import sys, os, getopt, csv
from kjbuckets import *
from heapq import heappush, heappop, heapreplace, heapify
from sets import Set
from graphlib import Graph

class go_informative_node:
	def __init__(self, go_association, go_index, size, type):
		self.goa_f = csv.reader(open(go_association, 'r'), delimiter='\t')
		self.goi_f = csv.reader(open(go_index, 'r'), delimiter='\t')
		self.size = int(size)
		self.type = int(type)
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
			if self.type:
				sys.stdout.write("%s\t%d\n"%(go_id, len(gene_id_list)))
			else:
				for gene_id in gene_id_list:
					sys.stdout.write('%s\t%s\t%s\n'%(go_id, go_index, gene_id))


class go_informative_node_bfs(object):
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							("go_association", 0, ): ['genome.gene2go', 'a', 1, 'the table or file mapping go_id to gene_id'],\
							("go_graph", 0, ): [None, 'i', 1, 'the file mapping go_id to go_index, OR the go_graph file. if not given, get graph from database.'],\
							('size', 1, int): [60, 's', 1, 'the minimum size for each informative node'],\
							('max_size', 1, int): [100, 'e', 1, 'the maximum size for each node, only for node_type=6'],\
							('max_no_of_nodes', 0, int): [160, 'm', 1, 'only nodes bigger than size above, only for node_type=4'],\
							('type', 0, int): [0, 't', 1, 'output format, 0(full), 1(stat), 2(between)'],\
							('node_type', 0, int): [2, 'n', 1, '1(all), 2(informative), 3(level 4), 4(biggest-first-break \
		control no of nodes), 5(biggest-first-break-level-by-level, similar to 4), 6(all nodes whose size is from size to max_size but not their children)'],\
							('level', 1, int):[4, 'l', 1, 'for node_type=3 or 4, this determines the level of the initial\
		function candidates 4(default)'],\
							('tax_id', 0, int): [3702, 'x', 1, 'Taxonomy ID to get gene position and coordinates.'],\
							('branch', 0, int):[2, 'c', 1, 'branch, 1(molecular_function), 2(biological_process), 3(cellular_component)'],\
							("output_fname", 0, ): [None, 'o', 1, 'if you wanna results into a file'],\
							('exclude_IEA', 0): [0, 'I', 0, 'Exclude gene2go entries whose evidence is IEA'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2008-08-22 use ProcessOptions
		03-08-05
			add a parameter, level to indicate the level of qualified nodes for node_type==3.
		11-06-05 add node_type=4
		12-30-05
			add branch
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)
		
		output_format_dict = {0:self.output_full,
			1:self.output_stat,
			2:self.output_between}
		if self.type in output_format_dict:
			self._output = output_format_dict[self.type]
		else:
			sys.stderr.write('Type %d invalid\n'%self.type)
			sys.exit(2)
		#12-30-05
		self.root_dict = {1:'GO:0003674',
			2:'GO:0008150',
			3:'GO:0005575'}
		self.root = self.root_dict[self.branch]	#12-30-05
		self.branch_name_dict = {1:'molecular_function',
			2:'biological_process',
			3:'cellular_component'}
		self.branch_unknown_acc_dict = {1:'GO:0005554',
			2:'GO:0000004',
			3:'GO:0008372'}
		
		
		self.go_graph = kjGraph()
		#11-06-05
		self.complement_go_graph = kjGraph()	#used to check parents of a node
		
		self.go_id2gene_id_dict = {}	#value utilizes kjSet data structure
		self.go_id2go_name = {}	#mapping between go_id an it's name
		self.go_id_descendent2gene_id_dict = {}
		self.informative_node_dict = {}
		self.go_id2depth = {}
		if self.debug:
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

	def dstruc_loadin_from_db(self, go_association, tax_id, go_branch, go_unknown_id, exclude_IEA=0, go_term_table='go.term', go_term2term_table='go.term2term'):
		"""
		2008-08-23
			fit into the new db infrastructure
			association info comes from 'go_association' table
			use tax_id to restrict genes
		11-07-05 add the complement_go_graph
		12-30-05
			flexible, adjut to self.branch
		03-11-06
			count #total genes
		"""
		sys.stderr.write("Loading Data Structure ...")
		#biological_process, not obsolete, not biological_process unknown.
		"""
		#12-30-05 deal with non-schema specific association tasks, here is human
		self.curs.execute("select a.go_id, a.gene_id from graph.association a, go.term t\
			where a.organism='Homo sapiens' and t.acc=a.go_id and\
			t.term_type='%s' and t.acc!='%s' and t.is_obsolete=0"%(self.branch_name_dict[self.branch],
			self.branch_unknown_acc_dict[self.branch]))
		"""
		#not obsolete, not unknown
		query = "select a.go_id, a.gene_id from %s a, %s t\
				where a.tax_id=%s and t.acc=a.go_id and\
				t.term_type='%s' and t.acc!='%s' and t.is_obsolete=0"%(go_association, go_term_table, tax_id, go_branch,
																						go_unknown_id)
		if exclude_IEA:
			query += " and evidence!='IEA'"
		self.curs.execute(query)
		
		#03-11-06 used to count #genes in total
		self.total_gene_id_set = Set()
		
		rows = self.curs.fetchall()
		for row in rows:
		#setup the go_id2gene_id_dict structure
			go_id = row[0]
			gene_id = row[1]
			self.total_gene_id_set.add(gene_id)	#03-11-06
			if go_id not in self.go_id2gene_id_dict:
				self.go_id2gene_id_dict[go_id] = kjSet([gene_id])
			else:
				self.go_id2gene_id_dict[go_id].add(gene_id)
		
		#get the non-obsolete biological_process GO DAG
		#12-30-05
		self.curs.execute("select t2t.term1_id, t2t.term2_id, t1.acc, t2.acc,t1.depth, t2.depth from \
			%s t2t, %s t1, %s t2 where t2t.term1_id=t1.id and \
			t2t.term2_id=t2.id and t1.is_obsolete=0 and t2.is_obsolete=0 and \
			t1.term_type='%s' and t2.term_type='%s' "%(go_term2term_table, go_term_table,\
													go_term_table, go_branch, go_branch))
		
		rows = self.curs.fetchall()
		for row in rows:
		#setup the go_graph structure
			self.go_graph.add((row[2], row[3]))
			#11-06-05
			self.complement_go_graph.add((row[3], row[2]))
			#12-30-05
			if row[4]:
				depth1 = int(row[4])
			else:
				depth1 = 0
			if row[5]:
				depth2 = int(row[5])
			else:
				depth2 = 0
			self.go_id2depth[row[2]] = depth1
			self.go_id2depth[row[3]] = depth2
		self.go_id_set = kjSet(self.go_graph.keys() + self.go_graph.values())
		
		#setup self.go_id2go_name
		self.curs.execute("select acc, name from %s"%go_term_table)
		rows = self.curs.fetchall()
		for row in rows:
			self.go_id2go_name[row[0]] = row[1]
		
		sys.stderr.write("%s total genes. Done\n"%(len(self.total_gene_id_set)))
	
	def keep_source_and_singletons_from_candidates(self, go_graph, go_id_descendent2gene_id_dict, candidates_set):
		"""
		03-23-06
			to avoid parent-child relationships in candidates
		"""
		graph_of_candidates = Graph.Graph()
		for go_id in candidates_set:
			graph_of_candidates.add_node(go_id)
		candidates_list = list(candidates_set)
		for i in range(len(candidates_list)):
			for j in range(i+1, len(candidates_list)):
				if go_graph.reachable(candidates_list[j]).member(candidates_list[i]):
					graph_of_candidates.add_edge(candidates_list[j], candidates_list[i])
					if self.debug:
						sys.stderr.write("In candidates: %s => %s\n"%(candidates_list[j], candidates_list[i]))
				if go_graph.reachable(candidates_list[i]).member(candidates_list[j]):
					graph_of_candidates.add_edge(candidates_list[i], candidates_list[j])
					if self.debug:
						sys.stderr.write("In candidates: %s => %s\n"%(candidates_list[i], candidates_list[j]))
		new_candidates_set = Set()
		new_candidates_queue = []
		for go_id in graph_of_candidates.node_list():
			if graph_of_candidates.inc_degree(go_id)==0:
				if self.debug:
					sys.stderr.write("new candidate: %s\n"%go_id)
				new_candidates_set.add(go_id)
				heappush(new_candidates_queue, [-len(go_id_descendent2gene_id_dict[go_id]), go_id])
		return new_candidates_queue, new_candidates_set
	
	def node_type_4(self, init_candidates_queue, init_candidates_set, go_graph, complement_go_graph, \
		go_id_descendent2gene_id_dict, max_no_of_nodes, size):
		"""
		11-06-05 total 5 data structures to check
		1. candidates_queue: big nodes to be investigated #a min queue, each entry is (-size, go_id), sorted by -size
		2. candidates_set: all nodes(big) wanted
		
		03-23-06 overhaul it to be much simpler
			get almost identical result as node_type 2(informative node), but longer time
		"""
		sys.stderr.write("Getting nodes for node type 4...\n")
		candidates_set = Set()
		while init_candidates_queue:
			node_size, go_id = heappop(init_candidates_queue)
			if node_size>=size:
				candidates_set.add(go_id)
		candidates_queue, candidates_set = self.keep_source_and_singletons_from_candidates(go_graph, go_id_descendent2gene_id_dict, candidates_set)
		
		if self.debug:
			sys.stderr.write("no of big_nodes %s\n"%len(candidates_queue))
			sys.stderr.write("candidates_queue %s\n"%repr(candidates_queue))
			sys.stderr.write("candidates_set %s\n"%repr(candidates_set))
		#following two mappings are used to decide whether some childs should be picked back
		trash_child_node2parent = {}
		parent2trash_child_node = {}
		terminal_good_node_list = []
		while len(candidates_set)<max_no_of_nodes and candidates_queue:
			node_size, go_id = heappop(candidates_queue)
			node_size = -node_size	#WATCH minus(-) ahead of node_size
			if self.debug:
				sys.stderr.write("Checking %s(size=%s)  no_of_big_nodes:%s...\n"%(go_id, node_size, len(candidates_set) ))
			neighbor_list = go_graph.neighbors(go_id)
			valid_neighbor = 0
			for neighbor in neighbor_list:
				neighbor_size = len(go_id_descendent2gene_id_dict[neighbor])
				if neighbor_size>=size:
					if self.debug:
						sys.stderr.write("One child %s exceeds the size %s.\n"%(neighbor, neighbor_size))
					valid_neighbor = 1
					candidates_set.add(neighbor)
			if valid_neighbor:
				#remove this node
				if self.debug:
					sys.stderr.write("%s removed from candidates_set.\n"%go_id)
				candidates_set.remove(go_id)
			else:	#real good nodes, no big child
				if self.debug:
					sys.stderr.write("%s moved from candidates_set to final set.\n"%go_id)
				candidates_set.remove(go_id)
				terminal_good_node_list.append(go_id)
			candidates_queue, candidates_set = self.keep_source_and_singletons_from_candidates(go_graph, go_id_descendent2gene_id_dict, candidates_set)
		
		#add those terminal nodes back into candidates_set
		for go_id in terminal_good_node_list:
			candidates_set.add(go_id)
		
		#finally, come to fruit
		informative_node_dict = {}
		for node in candidates_set:
				informative_node_dict[node] = 1
		
		sys.stderr.write("Done getting nodes for node type 4 with %s nodes.\n"%(len(candidates_set)) )
		return informative_node_dict
	
	def node_type_5(self, init_candidates_queue, init_candidates_set, go_graph, complement_go_graph, \
		go_id_descendent2gene_id_dict, max_no_of_nodes, size):
		"""
		11-06-05 similar to node_type_4 but, level by level and stop at level 5
		"""
		sys.stderr.write("Getting nodes for node type 5...\n")
		#remove childs among init_candidates_queue, start from the smallest one
		candidates_queue = []	#a min queue, each entry is (-size, go_id), sorted by -size
		candidates_set = Set()
		no_of_big_nodes = 0
		while init_candidates_queue:
			node_size, go_id = heappop(init_candidates_queue)
			parent_set = Set(complement_go_graph.neighbors(go_id))
			if parent_set&init_candidates_set:	#there's a parent
				init_candidates_set.remove(go_id)
			else:
				if node_size >= size:
					no_of_big_nodes += 1
					candidates_queue.append([-node_size, go_id])	#WATCH minus(-) ahead of node_size
				candidates_set.add(go_id)
		candidates_queue.sort()	#diff 4
		if self.debug:
			print "no_of_big_nodes", no_of_big_nodes
			print "candidates_queue", candidates_queue
			print "candidates_set", candidates_set
		#following two mappings are used to decide whether some childs should be picked back
		trash_child_node2parent = {}
		parent2trash_child_node = {}
		while no_of_big_nodes<max_no_of_nodes and candidates_queue:
			node_size, go_id = candidates_queue.pop(0)	#diff 4
			node_size = -node_size	#WATCH minus(-) ahead of node_size
			if self.debug:
				sys.stderr.write("Checking %s(size=%s)  no_of_big_nodes:%s...\n"%(go_id, node_size, no_of_big_nodes))
			neighbor_list = go_graph.neighbors(go_id)
			valid_neighbor = 0	#a flag to check whether there's a valid child node
			for neighbor in neighbor_list:
				neighbor_size = len(go_id_descendent2gene_id_dict[neighbor])
				if neighbor_size>=size:
					if self.debug:
						sys.stderr.write("One child %s exceeds the size %s.\n"%(neighbor, neighbor_size))
					valid_neighbor = 1
					break
			if valid_neighbor:
				#1st remove this node
				if self.debug:
					sys.stderr.write("%s removed from candidates_set.\n"%go_id)
				no_of_big_nodes -= 1	#previous condition check ensures that its size >self.size
				candidates_set.remove(go_id)
				#2nd release trash_child_node if possible
				if go_id in parent2trash_child_node:
					if self.debug:
						sys.stderr.write("%s appears in parent2trash_child_node.\n"%go_id)
					for trash_child_node in parent2trash_child_node[go_id]:
						trash_child_node2parent[trash_child_node].remove(go_id)
						no_of_parents = len(trash_child_node2parent[trash_child_node])
						if no_of_parents==0:
							if self.debug:
								sys.stderr.write("trash_child_node %s  has no parents now and revives.\n"%trash_child_node)
							del trash_child_node2parent[trash_child_node]
							trash_child_node_size = len(go_id_descendent2gene_id_dict[trash_child_node])
							if trash_child_node_size>=size:
								if self.debug:
									sys.stderr.write("trash_child_node %s has %s associated genes.\n"%(trash_child_node, trash_child_node_size))
								no_of_big_nodes += 1
								#candidates_queue.append([-trash_child_node_size, trash_child_node])	#WATCH minus(-) ahead of trash_child_node_size
							if trash_child_node_size>0:	#skip empty nodes
								candidates_set.add(trash_child_node)
					#remove it
					del parent2trash_child_node[go_id]
				#3rd pick good child
				for neighbor in neighbor_list:
					neighbor_size = len(go_id_descendent2gene_id_dict[neighbor])
					if neighbor_size==0:	#skip empty nodes
						continue
					parent_set = Set(complement_go_graph.neighbors(neighbor))
					parent_in_candidates = parent_set&candidates_set
					if parent_in_candidates:
						if self.debug:
							sys.stderr.write("unfortunately, %s 's child %s has parent %s in candidates_set.\n"%(go_id, neighbor, parent_in_candidates))
						#fill in trash_child_node2parent and parent2trash_child_node
						if neighbor not in trash_child_node2parent:
							trash_child_node2parent[neighbor] = Set()
						trash_child_node2parent[neighbor] |= parent_in_candidates
						for parent in parent_in_candidates:
							if parent not in parent2trash_child_node:
								parent2trash_child_node[parent] = Set()
							parent2trash_child_node[parent].add(neighbor)
					else:
						#good child
						if self.debug:
							sys.stderr.write("fortunately, %s 's child %s has no parent in candidates_set.\n"%(go_id, neighbor))
						if neighbor_size >= size and neighbor not in candidates_set:	#some of them are already checked int by other parents
							if self.debug:
								sys.stderr.write("furthermore, %s 's child %s has %s associated genes, into candidates_queue.\n"%(go_id, neighbor, neighbor_size))
							no_of_big_nodes += 1
							#candidates_queue.append([-neighbor_size, neighbor])	#WATCH minus(-) ahead of neighbor_size
						candidates_set.add(neighbor)
						#deal with trash_child_node
						neighbor_s_child_set = Set(go_graph.neighbors(neighbor))
						trash_child_node_set = Set(trash_child_node2parent)
						neighbor_s_trash_child_set = neighbor_s_child_set&trash_child_node_set
						if neighbor_s_trash_child_set:
							if self.debug:
								sys.stderr.write("oh, %s 's child %s has %s grand_trash_children, update two dicts.\n"%(go_id, neighbor, neighbor_s_trash_child_set))
							for trash_child_node in neighbor_s_trash_child_set:
								trash_child_node2parent[trash_child_node].add(neighbor)
							if neighbor not in parent2trash_child_node:
								parent2trash_child_node[neighbor] = Set()
							parent2trash_child_node[neighbor] |= neighbor_s_trash_child_set
		#finally, come to fruit
		informative_node_dict = {}
		for node in candidates_set:
			informative_node_dict[node] = 1
		sys.stderr.write("Done getting nodes for node type 5 with %s/%s nodes.\n"%(len(candidates_set), no_of_big_nodes))
		return informative_node_dict
	
	def node_type_2(self, root, min_size):
		"""
		2008-08-22
			split out of run()
		"""
		sys.stderr.write("Generating type 2 nodes ... ")
		#below is analogous to BFS,
		#informative_node_dict's role is similar to vertex coloring in real BFS
		self.informative_node_dict[root] = 1
		list_queue = [root]
		while list_queue:
			go_id = list_queue.pop(0)
			neighbor_list = self.go_graph.neighbors(go_id)
			for neighbor in neighbor_list:
				if len(self.go_id_descendent2gene_id_dict[neighbor]) >= min_size:
					if self.debug:
						self.log_file.write('%s ousted by %s\n'%(go_id,neighbor))
					self.informative_node_dict[go_id] = 0
					if self.informative_node_dict[neighbor] == 0:
					#trick is here. informative_node candidates and first touch
						list_queue.append(neighbor)
					self.informative_node_dict[neighbor] = 1
		sys.stderr.write("Done.\n")
	
	def node_type_6(self, root, min_size, max_size, go_graph, informative_node_dict, \
				go_id2gene_id_dict):
		"""
		2008-08-22
			all nodes whose size is from size to max_size but not their children
		"""
		sys.stderr.write("Generating type 6 nodes ... ")
		#below is analogous to BFS,
		#informative_node_dict's role is similar to vertex coloring in real BFS
		informative_node_dict[root] = 0
		list_queue = [root]
		while list_queue:
			go_id = list_queue.pop(0)
			neighbor_list = go_graph.neighbors(go_id)
			for neighbor in neighbor_list:
				no_of_associated_genes = len(go_id2gene_id_dict[neighbor])
				if no_of_associated_genes >= min_size and no_of_associated_genes<=max_size:
					informative_node_dict[neighbor] = 1
				else:
					if no_of_associated_genes >= min_size:
						list_queue.append(neighbor)
		sys.stderr.write("Done.\n")
	
	def prepareData(self):
		"""
		2008-08-22
			split out of run()
		"""
		sys.stderr.write("Preparing data ... ")
		#11-06-05
		if self.node_type==4 or self.node_type==5:
			init_candidates_queue = []	#a min queue, each entry is (size, go_id), sorted by size
			init_candidates_set = Set()
		else:
			init_candidates_queue = None
			init_candidates_set = None
			
		for go_id in self.go_id_set.items():
			#collect the associated genes from descendent terms
			descendent_set = self.go_graph.reachable(go_id)
			if go_id in self.go_id2gene_id_dict:
				self.go_id_descendent2gene_id_dict[go_id] = self.go_id2gene_id_dict[go_id]
			else:
				self.go_id_descendent2gene_id_dict[go_id] = kjSet()
			for descendent in descendent_set.items():
				if descendent in self.go_id2gene_id_dict:
					self.go_id_descendent2gene_id_dict[go_id] += self.go_id2gene_id_dict[descendent]
			if self.node_type==1:
				if go_id!=self.root and len(self.go_id_descendent2gene_id_dict[go_id]) > 0:
					#not biological_process and have >0 genes associated, then they are informative nodes.
					self.informative_node_dict[go_id] = 1
			elif self.node_type==2:
				self.informative_node_dict[go_id] = 0
			elif self.node_type==3:
				if go_id!=self.root and len(self.go_id_descendent2gene_id_dict[go_id]) > 0 and self.go_id2depth[go_id]==self.level:
					#not biological_process and have >0 genes associated and depth=4(level=4)
					self.informative_node_dict[go_id] = 1
			#11-06-05
			elif self.node_type==4 or self.node_type==5:
				node_size = len(self.go_id_descendent2gene_id_dict[go_id])
				if self.go_id2depth[go_id] == self.level and node_size>0:
					heappush(init_candidates_queue, [node_size, go_id])	#WATCH NO minus(-) ahead of node_size
					init_candidates_set.add(go_id)
		sys.stderr.write("Done.\n")
		return init_candidates_queue, init_candidates_set
	
	def run(self):
		"""
		2008-08-22
			restructure
		03-06-05
			add the node_type=3, level=5
		11-06-05 add node_type=4
		11-07-05 add node_type=5
		"""
		if os.path.isfile(self.go_association):
			#input from files
			self.goa_f = csv.reader(open(go_association, 'r'), delimiter='\t')
			self.gog_f = csv.reader(open(go_graph, 'r'), delimiter='\t')
			self.dstruc_loadin()
		else:
			if self.drivername=='postgres':
				#input from database
				try:
					import psycopg
				except:
					import psycopg2 as psycopg
				self.conn = psycopg.connect('host=%s dbname=%s user=%s passwd=%s'%(hostname, dbname, self.db_user, self.db_passwd))
				self.curs = self.conn.cursor()
				if self.schema:
					self.curs.execute("set search_path to %s"%self.schema)
			elif self.drivername=='mysql':
				import MySQLdb
				self.conn = MySQLdb.connect(db=self.dbname, host=self.hostname, user=self.db_user, passwd = self.db_passwd)
				self.curs = self.conn.cursor()
			self.dstruc_loadin_from_db(self.go_association, self.tax_id, self.branch_name_dict[self.branch], \
									self.branch_unknown_acc_dict[self.branch], self.exclude_IEA)
		
		init_candidates_queue, init_candidates_set = self.prepareData()
		
		"""
		node_size_ls = [len(value) for key, value in self.go_id2gene_id_dict.iteritems()]
		import pylab
		pylab.title("Histogram of GO node size")
		pylab.hist(node_size_ls, 300)
		pylab.show()
		"""
		if self.node_type==2:
			self.node_type_2(self.root, self.size)
		elif self.node_type==4:
			self.informative_node_dict = self.node_type_4(init_candidates_queue, init_candidates_set,\
				self.go_graph, self.complement_go_graph, self.go_id_descendent2gene_id_dict, self.max_no_of_nodes, self.size)
		elif self.node_type==5:
			self.informative_node_dict = self.node_type_5(init_candidates_queue, init_candidates_set,\
				self.go_graph, self.complement_go_graph, self.go_id_descendent2gene_id_dict, self.max_no_of_nodes, self.size)
		elif self.node_type==6:
			self.node_type_6(self.root, self.size, self.max_size, self.go_graph, self.informative_node_dict, \
				self.go_id_descendent2gene_id_dict)
		if self.output_fname:
			self.output(self.output_fname, self.go_id_descendent2gene_id_dict)
		"""
		#03-12-06 to see how the informative_node's are distributed, designed for node_type=2
		self.output_go_graph_above_informative_nodes(self.root, self.go_graph, self.informative_node_dict, \
			self.go_id_descendent2gene_id_dict, '/tmp/informative_node.dot', self.size)
		"""
		
	def output(self, output_fname, go_id2gene_id_dict):
		"""
		2008-08-22
			add output_fname
		03-11-06 count #genes in informative nodes
		"""
		sys.stderr.write("Outputting nodes ...")
		self.outf = open(output_fname, 'w')
		gene_id_in_informative_nodes_set = Set()
		no_of_go_nodes = 0
		for go_id,value in self.informative_node_dict.iteritems():
			if value == 1:
				no_of_go_nodes += 1
				gene_id_list = go_id2gene_id_dict[go_id].items()
				self._output(go_id, gene_id_list)
				#03-11-06
				for gene_id in gene_id_list:
					gene_id_in_informative_nodes_set.add(gene_id)
		sys.stderr.write("%s GO nodes retained. %s-%s=%s, #genes lost.\n"%(no_of_go_nodes, len(self.total_gene_id_set), len(gene_id_in_informative_nodes_set),\
			len(self.total_gene_id_set) - len(gene_id_in_informative_nodes_set)))
	
	def output_full(self, go_id, gene_id_list):
		for gene_id in gene_id_list:
			self.outf.write('%s\t%s\t%s\n'%(go_id, self.go_id2go_name[go_id], gene_id))

	def output_stat(self, go_id, gene_id_list):
		self.outf.write("%s\t%d\n"%(go_id, len(gene_id_list)))

	def output_between(self, go_id, gene_id_list):
		self.outf.write("%s\t%s\t%s\n"%(go_id, self.go_id2go_name[go_id], '\t'.join(gene_id_list)))

	def output_go_graph_above_informative_nodes(self, root, go_graph, informative_node_dict, \
			go_id_descendent2gene_id_dict, output_fname, size):
		"""
		03-12-06
			use the boost.graph to output part of GO graph which is above informative_node
		"""
		sys.stderr.write("Layout the informative_node's...")
		import boost.graph as bgl
		g = bgl.Graph()
		edge_weight_map = g.edge_property_map('integer')
		vertex_color_map = g.vertex_property_map('string')
		vertex_label_map = g.vertex_property_map('string')
		
		#informative_node_set is the boundary of the BFS search
		informative_node_set = Set()
		for go_id, value in informative_node_dict.iteritems():
			if value==1:
				informative_node_set.add(go_id)
		
		#this anti_vertex_label_map is used for add_edge() of g
		#give each vertex a node_id for graphviz
		anti_vertex_label_map = {}
		vertex_label2node_id = {}
		list_queue = [root]
		vertex_descriptor = g.add_vertex()
		anti_vertex_label_map[root] = vertex_descriptor
		vertex_label2node_id[root] = len(vertex_label2node_id) + 1
		vertex_label_map[vertex_descriptor] = str(vertex_label2node_id[root])
		
		
		#repeat the informative_node search algorithm
		another_informative_node_dict = {}
		another_informative_node_dict[root] = 1
		while list_queue:
			go_id = list_queue.pop(0)
			neighbor_list = go_graph.neighbors(go_id)
			add_neighbor_list_to_queue = 1	#an indicator whether to continue BFS on this neighbor_list
			for neighbor in neighbor_list:
				if neighbor not in anti_vertex_label_map:
					vertex_descriptor = g.add_vertex()
					anti_vertex_label_map[neighbor] = vertex_descriptor
					vertex_label2node_id[neighbor] = len(vertex_label2node_id) + 1
					vertex_label_map[vertex_descriptor] = str(vertex_label2node_id[neighbor])
					
				if neighbor in informative_node_set:	#one child of go_id is in informative_node_set, so stop BFS from this go_id
					vertex_color_map[anti_vertex_label_map[neighbor]] = 'red'
					add_neighbor_list_to_queue = 0
				if len(go_id_descendent2gene_id_dict[neighbor])>0:	#if no genes in this node, ignore the edge
					e = g.add_edge(anti_vertex_label_map[go_id], anti_vertex_label_map[neighbor])
					#the weight is used to show the percentage of genes inherited by the child
					weight = int(float(len(go_id_descendent2gene_id_dict[neighbor]))/len(go_id_descendent2gene_id_dict[go_id])*100)
					edge_weight_map[e] = weight
					#decide whether the neighbor should be added to list_queue(similar as in run())
					if len(go_id_descendent2gene_id_dict[neighbor]) >= size:
						another_informative_node_dict[go_id] = 0
						if neighbor not in another_informative_node_dict:
							#trick is here. informative_node candidates and first touch
							list_queue.append(neighbor)
						another_informative_node_dict[neighbor] = 1
		
		#'label', not 'weight' is drawn graphviz
		g.edge_properties['label'] = edge_weight_map
		#write_graphviz() needs 'node_id' map, not 'label'
		g.vertex_properties['node_id'] = vertex_label_map
		g.vertex_properties['color'] = vertex_color_map
		#output it in dot format
		g.write_graphviz(output_fname)
		sys.stderr.write("Done.\n")


if __name__ == '__main__':
	from pymodule import ProcessOptions
	main_class = go_informative_node_bfs
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()
