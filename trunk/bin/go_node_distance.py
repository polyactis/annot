#!/usr/bin/env python
"""
Usage: go_node_distance.py -k SCHEMA [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database to get the node_list
	-t ..., --table=...	table to store the node distances, node_distance(default)
	-b ..., --branch=...	which branch of GO, 0, 1(default) or 2
	-a, --all	compute distances between all nodes in that branch,
		instead of getting node_list from schema.go.
	-n,	--new_table	table is new. (create it first)
	-r, --report	report the progress(a number) IGNORE
	-c, --commit	commit the database transaction
	-l, --log	enable logging
	-h, --help              show this help

Examples:
	go_node_distance.py -k sc_38_all_no_info -l -t node_dist -n -r -c

Description:
	Program to compute three kinds of distances between two GO nodes and
	store them into a go table.
	branch illustration:
	0:	molecular_function
	1:	biological_process
	2:	cellular_component
"""

import sys, os, psycopg, getopt
from graphlib import Graph, GraphAlgo
from sets import Set

class distance_of_2nodes:
	def __init__(self):
		self.common_ancestor_set = Set()
		self.raw_distance = 10000
		self.lee_distance = 10000
		self.jasmine_distance = 10000

class go_node_distance:
	def __init__(self, hostname, dbname, schema, table, branch, all, new_table, report=0, \
		needcommit=0, log=0):
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.table = table
		
		self.all = int(all)
		self.new_table = int(new_table)
		self.report = int(report)
		self.needcommit = int(needcommit)
		self.log = int(log)
		
		#mapping for the branches
		self.branch_dict = {0:'molecular_function',
			1:'biological_process',
			2:'cellular_component'}
		self.branch = self.branch_dict[int(branch)]
		if self.log:
			self.log_file = open('/tmp/go_node_distance.log','w')
		#mapping between go term id and its index list
		self.go_id2index = {}
		#mapping between go_id and go_acc
		self.go_id2acc = {}
		#GO DAG (directed)
		self.go_digraph = Graph.Graph()
		#GO undirected Graph, to compute distance between two nodes
		self.go_graph = Graph.Graph()
		#the node_list contains the nodes to compute pairwise distances
		self.node_list = Set()
		#the list containing all indices
		self.index_list = Set()
		#key structure, mapping between go id and its associated distances
		self.go_id2distance = {}
		
	def dstruc_loadin(self):
		sys.stderr.write("Loading Data STructure...")
		
		#setup self.go_id2acc
		self.curs.execute("select id, acc from go.term where term_type='%s'"%(self.branch))
		rows = self.curs.fetchall()
		for row in rows:
			self.go_id2acc[row[0]] = row[1]
		
		#get the non-obsolete biological_process GO DAG
		self.curs.execute("select t2t.term1_id, t2t.term2_id, t1.acc, t2.acc from \
			go.term2term t2t, go.term t1, go.term t2 where t2t.term1_id=t1.id and \
			t2t.term2_id=t2.id and t1.is_obsolete=0 and t2.is_obsolete=0 and \
			t1.term_type='%s' and t2.term_type='%s' "%(self.branch, self.branch))
		rows = self.curs.fetchall()
		for row in rows:
			#setup the go_digraph and go_graph structure
			self.go_digraph.add_edge(row[0], row[1])
			self.go_graph.add_edge(row[0], row[1])
			self.go_graph.add_edge(row[1], row[0])

		#find the root based on different branches
		self.curs.execute("select id from go.term where name='%s'"%self.branch)
		rows = self.curs.fetchall()
		self.root = rows[0][0]

		#setup self.go_index2id and self.go_id2index
		self.go_index_setup()

		#setup the node list
		#biological_process unknown is not included.
		self.curs.execute("select t.id from go g, go.term t where g.go_id=t.acc and g.go_no!=0")
		rows = self.curs.fetchall()
		for row in rows:
			self.node_list |= self.go_id2index[row[0]]
		#convert the Set to list
		if self.all:
			self.node_list = list(self.index_list)
		else:
			self.node_list = list(self.node_list)
		
		sys.stderr.write("Done\n")
		
	def go_index_setup(self):
		'''
		each index is in a tuple form
		the indices of one node is in a Set form
		'''
		list_queue = [self.root]
		#one go term id may have multiple indices
		#',' is the trick, singleton tuple must have the trailing ','
		self.go_id2index[self.root] = Set( [(self.root,)] )

		while list_queue:
			go_id = list_queue.pop(0)
			neighbor_list = self.go_digraph.out_nbrs(go_id)
			for neighbor in neighbor_list:
				for index in self.go_id2index[go_id]:
					#convert index to list, append the neighbor, convert it back to tuple
					index_list = list(index)
					index_list.append(neighbor)
					new_index = tuple(index_list)
					#note: biological_process's index is not included
					self.index_list.add(new_index)
					#1. a node has two parents which are on the same level
					#2. a node has another parent which is on the same level of itself
					#3. a node has another parent which is on the level below its own level
					if neighbor in self.go_id2index:
						self.go_id2index[neighbor].add(new_index)
					else:
						self.go_id2index[neighbor] = Set([new_index])
				if neighbor not in list_queue:
					#it might have been put into the queue by previous nodes
					#because one node have multiple parents
					list_queue.append(neighbor)
			

	def run(self):
		if self.new_table:
			#first create the target_table
			try:
				self.curs.execute("create table go.%s(\
					go_id1	integer,\
					go_id2	integer,\
					raw_distance	integer,\
					lee_distance	integer,\
					jasmine_distance	integer,\
					common_ancestor_list	integer[])"%self.table)

			except:
				sys.stderr.write("Error occurred when creating table %s\n"%self.table)
				self.curs.execute("set search_path to %s"%schema)
		
		no_of_nodes = len(self.node_list)
		for i in range(no_of_nodes):
			for j in range(i+1, no_of_nodes):
				go_index1 = self.node_list[i]
				go_index2 = self.node_list[j]
				self.process_2indices(go_index1, go_index2)
				'''
				#deprecated
				lc_ancestors = self.lowest_common_ancestor(go_id1, go_id2)
				raw_dist = self.raw_distance(go_id1, go_id2)
				lee_dist = self.lee_distance(lc_ancestors)
				jasmine_dist = self.jasmine_distance(go_id1, go_id2, lc_ancestors)
				'''
		self.submit()
	
	def process_2indices(self, index1, index2):
		#compute three distances for these two indices
		go_id1 = index1[-1]
		go_id2 = index2[-1]
		if go_id1 == go_id2:
			#indices pointing to same go id, stop!
			return
		if self.log:
			self.log_file.write("index of %d: %s\n"%(go_id1, repr(index1)) )
			self.log_file.write("index of %d: %s\n"%(go_id2, repr(index2)) )
		if go_id1 < go_id2:
			#arrange the key in ascending order
			key = (go_id1, go_id2)
		else:
			key = (go_id2, go_id1)
		if key not in self.go_id2distance:
			self.go_id2distance[key] = distance_of_2nodes()
		depth1 = len(index1)
		depth2 = len(index2)
		min_depth = min(depth1, depth2)
		#find the length of the longest common head sequence
		for i in range(min_depth):
			if index1[i] != index2[i]:
				break
		#i points to the first different column
		lca = index1[i-1]
		#add the lca to the common_ancestor set of this pair
		self.go_id2distance[key].common_ancestor_set.add(lca)
		#lenth of the longest common sequence is lca's depth
		depth_of_lca = i
		distance1_from_lca = depth1 - depth_of_lca
		distance2_from_lca = depth2 - depth_of_lca
		raw_distance = distance1_from_lca + distance2_from_lca
		lee_distance = 15 - depth_of_lca
		jasmine_distance = min(distance1_from_lca, distance2_from_lca)
		
		if raw_distance < self.go_id2distance[key].raw_distance:
			self.go_id2distance[key].raw_distance = raw_distance
		if lee_distance < self.go_id2distance[key].lee_distance:
			self.go_id2distance[key].lee_distance = lee_distance
		if jasmine_distance < self.go_id2distance[key].jasmine_distance:
			self.go_id2distance[key].jasmine_distance = jasmine_distance

		if self.log:
			self.log_file.write(">%d\t%d\t%d\t%d\t%d\t%d\n"%(go_id1, go_id2, raw_distance, lee_distance,\
				jasmine_distance, lca))
			self.log_file.write("%s and %s:\n"%(self.go_id2acc[go_id1], self.go_id2acc[go_id2]))
			self.log_file.write("\traw_distance: %d\n"%raw_distance)
			self.log_file.write("\tlee_distance: %d\n"%lee_distance)
			self.log_file.write("\tjasmine_distance: %d\n"%jasmine_distance)
			self.log_file.write("\tOne lowest common ancestor: %s\n"%(self.go_id2acc[lca]))
					
	def submit(self):
		sys.stderr.write("Database transacting...")
		for key,value in self.go_id2distance.iteritems():
			go_id1 = key[0]
			go_id2 = key[1]
			raw_distance = value.raw_distance
			lee_distance = value.lee_distance
			jasmine_distance = value.jasmine_distance
			common_ancestor_list = list(value.common_ancestor_set)
			self.curs.execute("insert into go.%s(go_id1, go_id2, raw_distance, lee_distance, \
				jasmine_distance, common_ancestor_list) values(%d, %d, %d, %d, %d, ARRAY%s)"%\
				(self.table, go_id1, go_id2, raw_distance, lee_distance, jasmine_distance,\
				repr(common_ancestor_list)) )

		if self.needcommit:
			self.curs.execute("end")
		sys.stderr.write("done.\n")
			
	'''
	code below is deprecated. It's doing the same thing via graph algorithm. But much slow.
	Code above using the indices(actually the path tuples) is much faster.
	'''
	def depth_of_one_node(self, go_id):

		return len(GraphAlgo.shortest_path(self.go_digraph, self.root, go_id))
		
	def lowest_common_ancestor(self, go_id1, go_id2):
		subgraph1 = self.go_digraph.back_bfs_subgraph(go_id1)
		subgraph2 = self.go_digraph.back_bfs_subgraph(go_id2)
		set1 = Set(subgraph1.node_list())
		set2 = Set(subgraph2.node_list())
		intersection_set = set1 & set2
		lc_ancestors = []
		tuple_list = []
		for node in intersection_set:
			tuple_list.append((self.depth_of_one_node(node), node))
		#sort based on depth, first field in the tuple
		tuple_list.sort()
		#in ascending order
		tuple_list.reverse()
		max_depth = tuple_list[0][0]
		for (depth,node) in tuple_list:
			if depth == max_depth:
				lc_ancestors.append(node)
			elif depth < max_depth:
				break
		
		return lc_ancestors

	def raw_distance(self, go_id1, go_id2):
		distance = len(GraphAlgo.shortest_path(self.go_graph, go_id1, go_id2))-1
		return distance

	def lee_distance(self, lc_ancestors):
		#this distance is based on paper Lee2004a
		#all the nodes in the lc_ancestors list are of the same depth, so first one is enough
		depth = self.depth_of_one_node(lc_ancestors[0])
		return 15-depth
		
	def jasmine_distance(self, go_id1, go_id2, lc_ancestors):
		min_distance = 100
		for ancestor in lc_ancestors:
			distance1 = len(GraphAlgo.shortest_path(self.go_digraph, ancestor, go_id1))-1
			distance2 = len(GraphAlgo.shortest_path(self.go_digraph, ancestor, go_id2))-1
			if distance1 < min_distance:
				min_distance = distance1
			if distance2 < min_distance:
				min_distance = distance2
		return min_distance

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema=", "table=", \
		"branch=", "all", "new_table", "report", "commit", "log"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:b:anrcl", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = 'node_distance'
	branch = 1
	all = 0
	new_table = 0
	report = 0
	commit = 0
	log = 0
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
		elif opt in ("-b", "--branch"):
			branch = int(arg)
		elif opt in ("-a", "--all"):
			all = 1
		elif opt in ("-n", "--new_table"):
			new_table = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-l", "--log"):
			log = 1
			
	if schema:
		instance = go_node_distance(hostname, dbname, schema, table, branch, all, new_table, report, \
			commit, log)
		instance.dstruc_loadin()
		instance.run()
	else:
		print __doc__
		sys.exit(2)
