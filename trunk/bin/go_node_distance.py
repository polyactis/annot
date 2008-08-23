#!/usr/bin/env python
"""

Examples:
	#compute pairwise distance only for nodes in sc_54
	go_node_distance.py -k sc_54 -t node_dist -n -c
	
	#compute pairwise distance between all GO nodes
	go_node_distance.py  -t node_dist -n -c -a
	
	#compute the depth of all GO nodes
	go_node_distance.py -z localhost -d graphdb -e -r -c
	
Description:
	Program to compute three kinds of distances between two GO nodes and
	store them into a go table.
	branch illustration:
	0:	molecular_function
	1:	biological_process
	2:	cellular_component
	
	(02-17-05)
		an important bug appears in process_2indices(), if one index is parent of another index,
		the function can't correctly pick the lowest_common_ancestor. It picks the one just above
		that one.
		Now use the shorter index as fixed, and longer one as dynamic to be searched. Reverse
		check from the end of the fixed_index whether there's a common one in dynamic_index.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
sys.path.insert(0, os.path.expanduser('~/lib/python'))
sys.path.insert(0, os.path.join(os.path.expanduser('~/script')))
import sys, os, getopt
from graphlib import Graph, GraphAlgo
from sets import Set

class distance_of_2nodes:
	def __init__(self, raw_distance=10000, lee_distance=10000, jasmine_distance=10000):
		self.common_ancestor_set = Set()
		self.raw_distance = raw_distance
		self.lee_distance = lee_distance
		self.jasmine_distance = jasmine_distance

class go_node_distance:
	__doc__ = __doc__
	option_default_dict = {('drivername', 1,):['mysql', 'v', 1, 'which type of database? mysql or postgres', ],\
							('hostname', 1, ): ['papaya.usc.edu', 'z', 1, 'hostname of the db server', ],\
							('dbname', 1, ): ['stock_250k', 'd', 1, 'database name', ],\
							('schema', 0, ): [None, 'k', 1, 'database schema name', ],\
							('db_user', 1, ): [None, 'u', 1, 'database username', ],\
							('db_passwd', 1, ): [None, 'p', 1, 'database password', ],\
							("table", 0, ): ['node_distance', 't', 1, 'table to store the node distances'],\
							("depth", 0, ): [0, 'e', 0, 'compute the depth of nodes, default is pairwise distance'],\
							('all', 0, ): [0, 'a', 0, 'compute distances between all nodes in that branch,\
		instead of getting node_list from schema.go.'],\
							('new_table', 0, ): [0, 'n', 0, 'table is new. (create it first)'],\
							('log', 0, ): [0, 'l', 0, 'enable logging'],\
							('branch', 0, int):[2, '', 1, 'which GO branch, 1(molecular_function), 2(biological_process), 3(cellular_component)'],\
							("output_fname", 0, ): [None, 'o', 1, 'if you wanna results into a file'],\
							('commit', 0, int):[0, 'c', 0, 'commit the db operation.'],\
							('debug', 0, int):[0, 'b', 0, 'toggle debug mode'],\
							('report', 0, int):[0, 'r', 0, 'toggle report, more verbose stdout/stderr.']}
	def __init__(self, **keywords):
		"""
		2008-08-22
			use ProcessOptions
		"""
		from pymodule import ProcessOptions
		self.ad = ProcessOptions.process_function_arguments(keywords, self.option_default_dict, error_doc=self.__doc__, class_to_have_attr=self)

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
		
		self.table = table
		
		
		#debugging flags
		self.debug_process_2indices = int(self.log)
		#mapping for the branches
		self.branch_dict = {0:'molecular_function',
			1:'biological_process',
			2:'cellular_component'}
		self.branch = self.branch_dict[int(self.branch)]

		#mapping between go term id and its index list
		self.go_id2index = {}
		#mapping between go_id and go_acc
		self.go_id2acc = {}
		#GO DAG (directed)
		self.go_digraph = Graph.Graph()
		#GO undirected Graph, to compute distance between two nodes
		self.go_graph = Graph.Graph()
		#the list containing all indices
		self.index_list = Set()
		#the node_list to be pairwisely computed
		self.node_list = []
		
	def dstruc_loadin(self):
		"""
		03-04-05
			self.all decides where node_list comes.
			move the block of getting node_dist from schema.go to
			the case that self.all is not set.
		03-05-05
			self.all is deprecated
		03-06-05
			self.all is used again
		"""
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
			#setup the go_digraph and go_graph structure by term_id
			self.go_digraph.add_edge(row[0], row[1])
			self.go_graph.add_edge(row[0], row[1])
			self.go_graph.add_edge(row[1], row[0])

		#find the root based on different branches
		self.curs.execute("select id from go.term where name='%s'"%self.branch)
		rows = self.curs.fetchall()
		#it's also a term_id
		self.root = rows[0][0]
		
		if self.all==0:
			#get the node_list from schema.go
			self.curs.execute("select t.id from go g, go.term t where g.go_id=t.acc")
			rows = self.curs.fetchall()
			for row in rows:
				self.node_list.append(row[0])
		
		sys.stderr.write("Done\n")
		
	def go_index_setup(self):
		'''
		each index is in a tuple form
		the indices of one node is in a Set form
		'''
		sys.stderr.write("Computing graph indices of all GO nodes...")
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
		
		sys.stderr.write("Done\n")

	def node_depth(self):
		'''
		compute the depth of each node of GO tree from its index and put it into go.term.
		'''
		sys.stderr.write("Computing the depth of all GO nodes...")
		for (go_id, index_tuple_set) in self.go_id2index.iteritems():
			depth = 100
			index_list = []
			for index in index_tuple_set:
				index_acc_form = []
				for term_id in index:
					index_acc_form.append(self.go_id2acc[term_id])
				index_list.append('(%s)'%(','.join(index_acc_form)) )
				if len(index) < depth:
					depth = len(index)
			#Strings in a list will be wrapped by "'" after 'repr', which is not good for database submission
			#','.join(list) doesn't wrap "'" around the strings, which is good.
			self.curs.execute("update go.term set depth=%d, index='%s' where id=%d"%\
				(depth, ','.join(index_list), go_id) )
		sys.stderr.write("Done\n")

		
	def node_pairwise_distance(self):
		"""
		03-05-05
			lower the memory usage:
			instead of using the a data structure to hold all distances and submit at last.
			now compute one pair, submit one pair
		"""
		if self.new_table:
			sys.stderr.write("Creating new table...")
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
			sys.stderr.write("Done.\n")
			
		sys.stderr.write("Computing the pairwise distances of all GO nodes...")
	
		if self.all:
			go_id_list = self.go_id2index.keys()
		else:
			go_id_list = self.node_list
		no_of_nodes = len(go_id_list)
		for i in range(no_of_nodes):
			for j in range(i+1, no_of_nodes):
				pairwise_unit = distance_of_2nodes()
				go_id1 = go_id_list[i]
				go_id2 = go_id_list[j]
				if go_id1 < go_id2:
					key = (go_id1, go_id2)
				else:
					key = (go_id2, go_id1)
				index_set1 = self.go_id2index[go_id1]
				index_set2 = self.go_id2index[go_id2]
				for go_index1 in index_set1:
					for go_index2 in index_set2:
						self.process_2indices(go_index1, go_index2, pairwise_unit)
				
				if self.log:
					print "key is %s"%repr(key)
					print "raw_distance: %s"%pairwise_unit.raw_distance
					print "lee_distance: %s"%pairwise_unit.lee_distance
					print "jasmine_distance: %s"%pairwise_unit.jasmine_distance
					print "lca: %s"%repr(pairwise_unit.common_ancestor_set)
					raw_input("pause:")
				
				self.submit(key, pairwise_unit)
				'''
				#deprecated
				lc_ancestors = self.lowest_common_ancestor(go_id1, go_id2)
				raw_dist = self.raw_distance(go_id1, go_id2)
				lee_dist = self.lee_distance(lc_ancestors)
				jasmine_dist = self.jasmine_distance(go_id1, go_id2, lc_ancestors)
				'''
		sys.stderr.write("Done\n")
	
	def process_2indices(self, index1, index2, pairwise_unit):
		"""
		(02-17-05)
		an important bug appears in process_2indices(), if one index is parent of another index,
		the function can't correctly pick the lowest_common_ancestor. It picks the one just above that one.
		Now use the shorter index as fixed, and longer one as dynamic to be searched. Reverse
		check from the end of the fixed_index whether there's a common one in dynamic_index.
		03-05-05
			drop the global data structure, self.go_node_distance
			use the pairwise_unit instead
		
		"""
		if self.debug_process_2indices:
			print "\t\t##Enter process_2indices()"
		#compute three distances for these two indices
		go_id1 = index1[-1]
		go_id2 = index2[-1]
		if go_id1 == go_id2:
			#indices pointing to same go id, stop!
			return
		if self.debug_process_2indices:
			print "index of %d: %s\n"%(go_id1, repr(index1)) 
			print "index of %d: %s\n"%(go_id2, repr(index2))

		depth1 = len(index1)
		depth2 = len(index2)
		if depth1 <=  depth2:
			min_depth = depth1
			fixed_index = index1
			#convert to list because we want to use .index() method
			dynamic_index = list(index2)
			index_set_form = Set(index2)
		else:
			min_depth = depth2
			fixed_index = index2
			dynamic_index = list(index1)
			index_set_form = Set(index1)
		#find the length of the longest common head sequence
		lca_fixed_index = 0
		lca_dynamic_index = 0
		ran = range(min_depth)
		#start from the end
		ran.reverse()
		for i in ran:
			if fixed_index[i] in index_set_form:
				lca_fixed_index = i
				lca_dynamic_index = dynamic_index.index(fixed_index[i])
				break
		#i points to the first different column
		lca = fixed_index[lca_fixed_index]
		
		#add the lca to the common_ancestor set of this pair
		pairwise_unit.common_ancestor_set.add(lca)

		distance1_from_lca = len(fixed_index) - (lca_fixed_index+1)
		distance2_from_lca = len(dynamic_index) - (lca_dynamic_index+1)
		raw_distance = distance1_from_lca + distance2_from_lca
		lee_distance = 15 - (min(lca_fixed_index, lca_dynamic_index) +1)
		jasmine_distance = min(distance1_from_lca, distance2_from_lca)
		
		if raw_distance < pairwise_unit.raw_distance:
			pairwise_unit.raw_distance = raw_distance
			if self.debug_process_2indices:
				print "raw_distance replaced"
		if lee_distance < pairwise_unit.lee_distance:
			pairwise_unit.lee_distance = lee_distance
			if self.debug_process_2indices:
				print "lee_distance replaced"
		if jasmine_distance < pairwise_unit.jasmine_distance:
			pairwise_unit.jasmine_distance = jasmine_distance
			if self.debug_process_2indices:
				print "jasmine_distance replaced"

		if self.debug_process_2indices:
			print ">%d %d: %d(raw) %d(lee) %d(jasmine) %s(lca)\n"%(go_id1, go_id2, raw_distance, \
				lee_distance, jasmine_distance, lca)
			print "Corresponding GO acc: %s and %s:\n"%(self.go_id2acc[go_id1], self.go_id2acc[go_id2])
			print "\tOne lowest common ancestor: %s(%s)\n"%(lca, self.go_id2acc[lca])
			print "\t\t##Leave process_2indices()"
		
	def submit(self, key, pairwise_unit):
		"""
		03-02-05
			create two indices on go_id1 and go_id2 at the end of submit()
		03-05-05
			the index-creating part was split to another function, create_indices()
			now submit() is executed by one go_id pair after another.
		"""
		go_id1 = key[0]
		go_id2 = key[1]
		raw_distance = pairwise_unit.raw_distance
		lee_distance = pairwise_unit.lee_distance
		jasmine_distance = pairwise_unit.jasmine_distance
		common_ancestor_list = list(pairwise_unit.common_ancestor_set)
		self.curs.execute("insert into go.%s(go_id1, go_id2, raw_distance, lee_distance, \
			jasmine_distance, common_ancestor_list) values(%d, %d, %d, %d, %d, ARRAY%s)"%\
			(self.table, go_id1, go_id2, raw_distance, lee_distance, jasmine_distance,\
			repr(common_ancestor_list)) )
	
	def create_indices(self):
		"""
		03-05-05
			split from submit()
			add an index for (go_id1, go_id2)
		"""
		sys.stderr.write("Creating three indices...")
		#create indices
		self.curs.execute("set search_path to go")
		self.curs.execute("create index %s_go_id1_idx on go.%s(go_id1)"%\
			(self.table, self.table))
		self.curs.execute("create index %s_go_id2_idx on go.%s(go_id2)"%\
			(self.table, self.table))
		self.curs.execute("create index %s_go_id1_go_id2_idx on go.%s(go_id1, go_id2)"%\
			(self.table, self.table))
		sys.stderr.write("Done.\n")


	def run(self):
		"""
		03-05-05
		03-05-05
		run()
			--dstruc_loadin
			--go_index_setup
			
			
			--node_depth
			or
			--node_pairwise_distance
				--process_2indices
				--submit
			--create_indices
			
		In this class, go_id refers to the term id in tables of schema go.
		the index for one go node is its path in tuple form
		"""
		self.dstruc_loadin()
		#setup self.go_id2index
		self.go_index_setup()
		
		if self.depth:
			self.node_depth()
		else:
			self.node_pairwise_distance()
			self.create_indices()
		if self.needcommit:
			self.curs.execute("end")
		else:
			print "not commited into database"


	'''
	code below is deprecated. It's doing the same thing via graph algorithm. But much slower.
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
	from pymodule import ProcessOptions
	main_class = go_node_distance
	po = ProcessOptions(sys.argv, main_class.option_default_dict, error_doc=main_class.__doc__)
	
	instance = main_class(**po.long_option2value)
	instance.run()