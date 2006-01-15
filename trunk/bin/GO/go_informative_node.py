#!/usr/bin/env python
"""
Usage:	go_informative_node.py -a GO_ASSOCIATION -i GO_INDEX [OPTIONS}
	go_informative_node.py -k SCHEMA -b [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-a ..., -go_association=...	the file mapping go_id to gene_id
	-i ..., --go_index=...	the file mapping go_id to go_index, OR the go_graph file
	-s ..., --size=...	the size of the informative node, 60(default)
	-m ...,	max_no_of_nodes, only nodes bigger than size above, 150(default), node_type=4
	-t ..., --type=...	output format, 0(default, full), 1(stat), 2(between)
	-n ..., --node_type=...	1(all, default), 2(informative), 3(level 4), 4(biggest-first-break 
		control no of nodes), 5(biggest-first-break-level-by-level, similar to 4)
	-l ..., --level=...	for node_type=3 or 4, this determines the level of a qualified go, 4(default)
	-c ...,	branch, 1(molecular_function), 2(biological_process, default), 3(cellular_component)
	-u,	debug
	-h, --help      show this help
	
Examples:
	go_informative_node.py -i process.pathindex.annotation.txt -a go2ug.txt
	go_informative_node.py -i term2term.txt -a go2ug.txt -b
	go_informative_node.py -i term2term.txt -a go2ug.txt -b -t 1
	go_informative_node.py -k hs_yh60_42 -b
	go_informative_node.py -k sc_38_all_no_info -n 2 -b

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
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import sys, os, psycopg, getopt, csv
from kjbuckets import *
from heapq import heappush, heappop, heapreplace
from sets import Set

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


class go_informative_node_bfs:
	"""
	03-08-05
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, go_association=None, \
		go_graph=None, size=60, max_no_of_nodes=160,\
		type=0, node_type=1, level=4, branch=2, debug=0):
		"""
		03-08-05
			add a parameter, level to indicate the level of qualified nodes for node_type==3.
		11-06-05 add node_type=4
		12-30-05
			add branch
		"""
		self.schema = schema
		if self.schema:
		#input from database
			self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
			self.curs = self.conn.cursor()
			self.curs.execute("set search_path to %s"%self.schema)
		else:
		#input from files
			self.goa_f = csv.reader(open(go_association, 'r'), delimiter='\t')
			self.gog_f = csv.reader(open(go_graph, 'r'), delimiter='\t')
		self.size = int(size)
		self.max_no_of_nodes = int(max_no_of_nodes)
		self.type = int(type)
		self.node_type = int(node_type)
		self.level = int(level)
		self.branch = int(branch)
		self.debug = int(debug)
		
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
		"""
		11-07-05 add the complement_go_graph
		12-30-05
			flexible, adjut to self.branch
		"""
		sys.stderr.write("Loading Data STructure...")
		#biological_process, not obsolete, not biological_process unknown.
		"""
		#12-30-05 deal with non-schema specific association tasks, here is human
		self.curs.execute("select a.go_id, a.gene_id from graph.association a, go.term t\
			where a.organism='Homo sapiens' and t.acc=a.go_id and\
			t.term_type='%s' and t.acc!='%s' and t.is_obsolete=0"%(self.branch_name_dict[self.branch],
			self.branch_unknown_acc_dict[self.branch]))
		"""
		
		self.curs.execute("select a.go_id, g.gene_id from graph.association a, go.term t,\
			gene g where g.gene_id=a.gene_id and t.acc=a.go_id and\
			t.term_type='%s' and t.acc!='%s' and t.is_obsolete=0"%(self.branch_name_dict[self.branch],
			self.branch_unknown_acc_dict[self.branch]))
		
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
		#12-30-05
		self.curs.execute("select t2t.term1_id, t2t.term2_id, t1.acc, t2.acc,t1.depth, t2.depth from \
			go.term2term t2t, go.term t1, go.term t2 where t2t.term1_id=t1.id and \
			t2t.term2_id=t2.id and t1.is_obsolete=0 and t2.is_obsolete=0 and \
			t1.term_type='%s' and t2.term_type='%s' "%(self.branch_name_dict[self.branch],
			self.branch_name_dict[self.branch]))
		
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
		self.curs.execute("select acc, name from go.term")
		rows = self.curs.fetchall()
		for row in rows:
			self.go_id2go_name[row[0]] = row[1]
		
		sys.stderr.write("Done\n")
	
	def node_type_4(self, init_candidates_queue, init_candidates_set, go_graph, complement_go_graph, \
		go_id_descendent2gene_id_dict, max_no_of_nodes, size):
		"""
		11-06-05 total 5 data structures to check
		1. candidates_queue: big nodes to be investigated
		2. candidates_set: all good nodes(big + small)
		3. no_of_big_nodes: total big nodes
		4. trash_child_node2parent: child node discovered but has parent in candidates_set
		5. parent2trash_child_node: reverse map of 4
		"""
		sys.stderr.write("Getting nodes for node type 4...\n")
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
					heappush(candidates_queue, [-node_size, go_id])	#WATCH minus(-) ahead of node_size
				candidates_set.add(go_id)
		if self.debug:
			print "no_of_big_nodes", no_of_big_nodes
			print "candidates_queue", candidates_queue
			print "candidates_set", candidates_set
		#following two mappings are used to decide whether some childs should be picked back
		trash_child_node2parent = {}
		parent2trash_child_node = {}
		while no_of_big_nodes<max_no_of_nodes and candidates_queue:
			node_size, go_id = heappop(candidates_queue)
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
								heappush(candidates_queue, [-trash_child_node_size, trash_child_node])	#WATCH minus(-) ahead of trash_child_node_size
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
							heappush(candidates_queue, [-neighbor_size, neighbor])	#WATCH minus(-) ahead of neighbor_size
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
		sys.stderr.write("Done getting nodes for node type 4 with %s/%s nodes.\n"%(len(candidates_set), no_of_big_nodes))
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
		sys.stderr.write("Done getting nodes for node type 4 with %s/%s nodes.\n"%(len(candidates_set), no_of_big_nodes))
		return informative_node_dict
	
	def run(self):
		"""
		03-06-05
			add the node_type=3, level=5
		11-06-05 add node_type=4
		11-07-05 add node_type=5
		"""
		if self.schema:
		#input from database	
			self.dstruc_loadin_from_db()
		else:
		#input from files
			self.dstruc_loadin()
		#11-06-05
		if self.node_type==4 or self.node_type==5:
			init_candidates_queue = []	#a min queue, each entry is (size, go_id), sorted by size
			init_candidates_set = Set()
			
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
		
		if self.node_type==2:
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
		elif self.node_type==4:
			self.informative_node_dict = self.node_type_4(init_candidates_queue, init_candidates_set,\
				self.go_graph, self.complement_go_graph, self.go_id_descendent2gene_id_dict, self.max_no_of_nodes, self.size)
		elif self.node_type==5:
			self.informative_node_dict = self.node_type_5(init_candidates_queue, init_candidates_set,\
				self.go_graph, self.complement_go_graph, self.go_id_descendent2gene_id_dict, self.max_no_of_nodes, self.size)
		self.output()
	
	def output(self):
		for go_id,value in self.informative_node_dict.iteritems():
			if value == 1:	
				gene_id_list = self.go_id_descendent2gene_id_dict[go_id].items()
				self._output(go_id, gene_id_list)
	
	def output_full(self, go_id, gene_id_list):
		for gene_id in gene_id_list:
			sys.stdout.write('%s\t%s\t%s\n'%(go_id, self.go_id2go_name[go_id], gene_id))

	def output_stat(self, go_id, gene_id_list):
		sys.stdout.write("%s\t%d\n"%(go_id, len(gene_id_list)))

	def output_between(self, go_id, gene_id_list):
		sys.stdout.write("%s\t%s\t%s\n"%(go_id, self.go_id2go_name[go_id], '\t'.join(gene_id_list)))

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:a:i:s:m:t:n:l:c:bu", ["help", "hostname=", \
			"dbname=", "schema=", "go_association=", "go_index=", "size=", "type=", "node_type=", "level=", "bfs"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''	
	go_association = ''
	go_index = ''
	size = 60
	max_no_of_nodes = 160
	type = 0
	node_type = 1
	level = 4
	branch = 2
	bfs = 1
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
		elif opt in ("-a", "--go_association"):
			go_association = arg
		elif opt in ("-i", "--go_index"):
			go_index = arg
		elif opt in ("-s", "--size"):
			size = int(arg)
		elif opt in ("-m",):
			max_no_of_nodes = int(arg)
		elif opt in ("-t", "--type"):
			type = int(arg)
		elif opt in ("-n", "--node_type"):
			node_type = int(arg)
		elif opt in ("-l", "--level"):
			level = int(arg)
		elif opt in ("-c",):
			branch = int(arg)
		elif opt in ("-b", "--bfs"):
			bfs = 1
		elif opt in ("-u",):
			debug = 1
	if schema or (go_association and go_index):
		instance = go_informative_node_bfs(hostname, dbname, schema, go_association, \
			go_index, size, max_no_of_nodes, type, node_type, level, branch, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
