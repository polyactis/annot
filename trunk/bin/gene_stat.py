#!/usr/bin/env python
"""
Usage: gene_stat.py -k SCHEMA [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	cluster_stat(default)
	-m ..., --mcl_table=...	mcl_result(default), mcl_result table corresponding to above table.
	-g ..., --gene_table=...	table to store the stat results, p_gene(default), needed if commit
	-e ..., --depth_cut_off=...	the minimum depth for a go node to be valid, 3(default)
	-f ..., --dir_files=...	the directory containing all the files outputed by cluster_stat.py
	-x ..., --recurrence_gap_size=...	2(default, IGNORE)
	-y ..., --connectivity_gap_size=...	2(default, IGNORE)
	-l, --leave_one_out	use the leave_one_out stat method, default is no leave_one_out
	-w, --wu	Wu's strategy(Default is Jasmine's strategy)
	-r, --report	report the progress(a number)
	-c, --commit	commit the database transaction, records in table gene.
	-q ..., --subgraph_cut_off=...	the cut_off for the subgraph to be valid in one dataset, 0(default)
		NOTICE: 0 means the binary conversion won't be used, just summing the floats.
	-b, --debug	enable debugging, no debug by default
	-h, --help              show this help

Examples:
	gene_stat.py -k sc_54 -t cluster_stat2 -m mcl_result2 -g p_gene_2 -e 5 -l -w 

Description:
	02-21-05
		a slim version of gene_stat_plot.py, only does _gene_stat_leave_one_out and
	submit(), another part goes to p_gene_analysis.py
	

"""

import sys, os, psycopg, getopt, csv, fileinput, math
from sets import Set
from codense.common import db_connect
from numarray import greater_equal

class gene_stat:
	"""
	run()
		--dstruc_loadin()
		--core_from_files()
		or
		--core()
			--_gene_stat_leave_one_out()
				--index_tuple()
				--direct_match()
				--L1_match()
				--common_ancestor_deep_enough()
		--submit()
	
	03-08-05
		add two more parameters, recurrence_gap_size and connectivity_gap_size (make them explicit)
	03-14-05
		distance loading on the fly
	03-15-05
		several functions were changed because the need to record down the lowest common ancestors 
		between the predicted function and known functions.
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, table=None, \
		mcl_table=None, leave_one_out=0, wu=0, report=0,\
		depth_cut_off =3, dir_files=None, needcommit=0, gene_table='p_gene',\
		subgraph_cut_off=0, debug=0, recurrence_gap_size=2, connectivity_gap_size=2):
		
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.table = table
		self.mcl_table = mcl_table
		self.leave_one_out = int(leave_one_out)
		self.wu = int(wu)
		self.report = int(report)
		self.depth_cut_off = int(depth_cut_off)
		self.dir_files = dir_files
		self.needcommit = int(needcommit)
		self.gene_table = gene_table
		self.subgraph_cut_off = float(subgraph_cut_off)
		#debugging flag
		self.debug = int(debug)
		
		#debug flags in several functions
		self.debug_L1_match = 0
		self.debug_common_ancestor_deep_enough = 0
		
		self.distance_table = 'go.node_dist'
		
		#the gap between two recurrences
		self.recurrence_gap_size = int(recurrence_gap_size)
		self.connectivity_gap_size = int(connectivity_gap_size)
		
		#mapping between gene_no and go_no set
		self.known_genes_dict = {}
		#the dictionary having (recurrence, connectivity) as a key, 
		# [[p_value, cluster_id, gene_no, go_no, is_correct, is_correct_L1, is_correct_lca, cluster_size, unknown_gene_ratio], [...], ... ] as a value
		self.prediction_tuple2list ={}
		#mapping between a pair of go_no's and its associated distances
		self.go_no2distance = {}
		#mapping each go_no to its depth
		self.go_no2depth = {}
		#mapping each go term id to its depth
		self.go_term_id2depth = {}
		self.go_term_id2go_no = {}
		self.go_no2go_id = {}
		self.go_no2term_id = {}
		self.no_of_records = 0

	def dstruc_loadin(self, curs):
		"""
		03-14-05
			remove the distance loading part
		"""
		sys.stderr.write("Loading Data STructure...\n")
		from codense.common import get_known_genes_dict, get_go_no2go_id,\
			get_go_no2term_id, get_go_no2depth, get_go_term_id2go_no, \
			get_go_term_id2depth
		
		self.known_genes_dict = get_known_genes_dict(curs)
		self.go_no2go_id = get_go_no2go_id(curs)
		self.go_no2term_id = get_go_no2term_id(curs)
		self.go_no2depth = get_go_no2depth(curs)
		self.go_term_id2go_no = get_go_term_id2go_no(curs)
		self.go_term_id2depth = get_go_term_id2depth(curs)
		
		sys.stderr.write("Done\n")

	def core_from_files(self):
		#following codes are attaching directory path to each file in the list
		file_list = os.listdir(self.dir_files)
		file_path_list = []
		for filename in file_list:
			file_path_list.append(os.path.join(self.dir_files, filename))
		#multiple files constitute the source of data
		self.files = fileinput.input(file_path_list)
		#wrap it with a reader
		self.reader = csv.reader(self.files, delimiter='\t')
		for row in self.reader:
			row[0] = int(row[0])
			row[1] = int(row[1])
			row[3] = float(row[3])
			curs.execute("select recurrence_array, vertex_set from %s where mcl_id=%d"%(self.mcl_table, int(row[0])) )
			rows = curs.fetchall()
			#first append the recurrence_array
			row.append(rows[0][0])
			#second append the vertex_set
			row.append(rows[0][1])
			#only leave_one_out
			self._gene_stat_leave_one_out(row)

			if self.report and self.no_of_records%2000==0:
				sys.stderr.write('%s%s'%('\x08'*20, self.no_of_records))
		if self.report:
			sys.stderr.write('%s%s'%('\x08'*20, self.no_of_records))

	def core(self, curs):
		"""
		03-14-05
			load go_no2distance on demand
		"""
		sys.stderr.write("Starting gene-stat...\n")
		from gene_p_map_redundancy import gene_p_map_redundancy
		node_distance_class = gene_p_map_redundancy()
		#the central function of the class
		if self.leave_one_out:
			#leave_one_out method gets data from both cluster_stat-like and mcl_result-like table
			curs.execute("DECLARE crs CURSOR FOR select c.mcl_id, c.leave_one_out, c.p_value_vector, \
				 c.connectivity, m.recurrence_array, m.vertex_set from %s c, %s m where c.mcl_id=m.mcl_id"\
				%(self.table, self.mcl_table))
		else:
			#no leave_one_out method gets data only from mcl_result-like table
			curs.execute("DECLARE crs CURSOR FOR select mcl_id, vertex_set, p_value_min, go_no_vector, unknown_gene_ratio, \
				recurrence_array from %s where connectivity>=%f and p_value_min notnull and array_upper(recurrence_array, 1)>=%d\
				and array_upper(vertex_set, 1)<=%d"%(self.mcl_table, self.connectivity_cut_off, self.recurrence_cut_off, self.cluster_size_cut_off))
		
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		while rows:
			for row in rows:
				if self.leave_one_out:
					#in leave_one_out, only one gene's function is predicted based on one row
					self._gene_stat_leave_one_out(row, node_distance_class, curs)
				else:
					#in no leave_one_out, function of all vertices in that cluster is predicted based on one row
					self._gene_stat_no_leave_one_out(row)
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, self.no_of_records))
			
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		sys.stderr.write("Done.\n")
	
	def _gene_stat_no_leave_one_out(self, row):
		"""
		later
		"""
		pass
	
	def _gene_stat_leave_one_out(self, row, node_distance_class=None, curs=None):
		"""
		03-08-05
			set a default(1.0) for min_p_value
			fix a bug, alter >self.depth_cut_off to >= self.depth_cut_off
			
		03-08-05
			don't take floor of recurrence and connectivity anymore, p_gene_analysis.py and p_gene_lm.py will take care of this.
		
		03-08-05
			fix another important bug
			when looking for other functions that have same min_p_value, the depth_cut_off requirement is forgotten.
		
		03-14-05
			getting go_no2distance via node_distance_class
		
		03-15-05
			get lca's for common_ancestor_deep_enough approach
		"""
		
		mcl_id = row[0]
		gene_no = row[1]
		p_value_vector = row[2][1:-1].split(',')
		connectivity = float(row[3])
		recurrence_array = row[4][1:-1].split(',')
		recurrence_array = map(float, recurrence_array)
		vertex_set = row[5][1:-1].split(',')
		vertex_set = map(int, vertex_set)
		
		#subgraph_cut_off is the cutoff for a cluster to be counted as occurred
		if self.subgraph_cut_off!=0:
			#0 means no cutoff
			recurrence_array = greater_equal(recurrence_array, self.subgraph_cut_off)
		"""
		#take the floor of the recurrence
		recurrence = int(math.floor(sum(recurrence_array)/self.recurrence_gap_size)*self.recurrence_gap_size)
		#take the floor of the connectivity *10
		connectivity = int(math.floor(connectivity*10/self.connectivity_gap_size)*self.connectivity_gap_size)
		"""
		recurrence = sum(recurrence_array)

		#setup in prediction_tuple2list
		prediction_tuple = (recurrence, connectivity)
		if prediction_tuple not in self.prediction_tuple2list:
			self.prediction_tuple2list[prediction_tuple] = []
		
		#default min_p_value is 1.0
		min_p_value = 1.0
		
		#transform into float type
		p_value_index_tuple_list = self.index_tuple(p_value_vector)
		for (p_value, index) in p_value_index_tuple_list:
			if self.wu:
				#index 0 corresponds to go_no 0.
				go_no = index
			else:
				#index 0 corresponds to go_no 1
				go_no = index+1
			if self.go_no2depth[go_no] >= self.depth_cut_off:
				min_p_value = p_value
				break

		#p-value 1.0 means the corresponding function has no associated genes in the cluster.
		#this situation means the cluster's associated functions are above the depth_cut_off
		if min_p_value >= 1.0:
			return
			
		if self.wu:
			unknown_gene_ratio = float(p_value_vector[0])
		else:
			unknown_gene_ratio = -1
		#The cluster is an eligible cluster. Passing all the cut_offs.
		#
		self.no_of_records += 1
		
		#looking for go_nos that have the same min_p_value
		for (p_value, index) in p_value_index_tuple_list:
			if p_value > min_p_value:
				break
			elif index == 0 or p_value==1.0:
				#0 is the unknown function, this is almost impossible because its depth = 2(see condition above)
				#1.0 is for function that has no associated genes
				continue
			elif p_value == min_p_value:
				if self.wu:
					#index 0 corresponds to go_no 0.
					go_no = index
				else:
					#index 0 corresponds to go_no 1
					go_no = index+1
				if self.go_no2depth[go_no] < self.depth_cut_off:
					#again we need the go_no to be deep enough
					continue
				
				if gene_no in self.known_genes_dict:
					k_functions_set = self.known_genes_dict[gene_no]
					is_correct = self.direct_match(go_no, k_functions_set)
					is_correct_L1 = self.L1_match(go_no, k_functions_set, node_distance_class, curs)
					is_correct_lca = self.common_ancestor_deep_enough(go_no, k_functions_set, node_distance_class, curs)
				else:
					#unknown gene
					is_correct = -1
					is_correct_L1 = -1
					is_correct_lca = -1
					#clear lca_list
					self.lca_list = []
					
				prediction_list = [p_value, mcl_id, gene_no, go_no, is_correct, is_correct_L1, \
					is_correct_lca, len(vertex_set), unknown_gene_ratio, self.lca_list]
				self.prediction_tuple2list[prediction_tuple].append(prediction_list)

	def index_tuple(self, list):
		"""
		03-08-05
			input: list
			output: a new list, each element is a tuple of (value, index)
			
			value is from list, index is value's index in the list. The new list is sorted.
		"""
		new_list = []
		for i in range(len(list)):
			#value is position 0, and index is position 1
			new_list.append((float(list[i]), i))
		#the sort is based on position 0
		new_list.sort()
		return new_list

	def direct_match(self, p_go_no, k_functions_set, node_distance_class=None, curs=None):
		if self.go_no2depth[p_go_no] < self.depth_cut_off:
			#first see if it's deep enough
			return 0
		else:
			if p_go_no in k_functions_set:
				return 1
			else:
				return 0
	
	def L1_match(self, p_go_no, k_functions_set, node_distance_class=None, curs=None):
		"""
		03-14-05
			getting go_no2distance via node_distance_class
		"""
		if self.debug_L1_match:
			print "\t\t ### In function L1_match() "
		#default not match
		flag = 0
		if self.go_no2depth[p_go_no] < self.depth_cut_off:
			#not good 
			return 0
		for k_go_no in k_functions_set:
			if self.go_no2depth[k_go_no] < self.depth_cut_off:
				#the known function is above the depth_cut_off, discard it
				continue
			elif k_go_no == p_go_no:
				flag = 1
				break
			elif k_go_no < p_go_no:
				key = (k_go_no, p_go_no)
			elif  k_go_no > p_go_no:
				key = (p_go_no, k_go_no)
			if key in self.go_no2distance:
				jasmine_distance = self.go_no2distance[key][2]
			else:
				jasmine_distance = node_distance_class.get_distance(curs, k_go_no, p_go_no, self.distance_table, \
					self.go_no2distance, self.go_no2term_id)
			if jasmine_distance == 1:
				#jasmine distance = 1
				if self.debug_L1_match:
					print 'One of %s and %s are one step away from their lowest common ancestor, \
						with depth_cut_off, %d'%(p_go_no, k_go_no, self.depth_cut_off)
					raw_input("Pause:")
				flag = 1
				break
			
		if self.debug_L1_match:
			print "\t\t ###leave function L1_match()"
		return flag
	
	def common_ancestor_deep_enough(self, p_go_no, k_functions_set, node_distance_class=None, curs=None):
		"""
		03-10-05
			add three lines in the beginning to check if p_go_no is deep enough.
			find a bug, in case of k_go_no==p_go_no, the function should return 1.
		03-14-05
			getting go_no2distance via node_distance_class
		03-15-05
			add a feature to get the lca's between p_go_no and k_functions_set
		"""
		if self.debug_common_ancestor_deep_enough:
			print "\t\t ### Enter common_ancestor_deep_enough() "
		#the lowest common ancestor is blank list by default, use global variable not to change the interface.
		self.lca_list = []
		
		if self.go_no2depth[p_go_no] < self.depth_cut_off:
			#not good
			if self.debug_common_ancestor_deep_enough:
				print "%s doesn't pass the depth_cut_off.\n"%(self.go_no2go_id[p_go_no])
			self.lca_list = []
			return 0
		#each element is a tuple, like (depth, go_no)
		ancestor_set = Set()
		for k_go_no in k_functions_set:
			if k_go_no == p_go_no:
				#bug here corrected.
				if self.debug_common_ancestor_deep_enough:
					print "%s hits itself.\n"%(self.go_no2go_id[p_go_no])
				self.lca_list = [p_go_no]
				return 1
			elif k_go_no < p_go_no:
				key = (k_go_no, p_go_no)
			elif k_go_no > p_go_no:
				key = (p_go_no, k_go_no)
			if key not in self.go_no2distance:
				#load in the distances and common_ancestors
				node_distance_class.get_distance(curs, k_go_no, p_go_no, self.distance_table, \
					self.go_no2distance, self.go_no2term_id)
				
			for ancestor in self.go_no2distance[key][3]:
				#NOTE: ancestor is in term_id. neither go_no nor go_id
				depth = self.go_term_id2depth[ancestor]
				if depth >= self.depth_cut_off:
					ancestor_set.add((depth, ancestor))
					if self.debug_common_ancestor_deep_enough:
						print "%s's common_ancestor %s passes the depth_cut_off, with depth %d.\n"%\
						(self.go_no2go_id[p_go_no], self.go_no2go_id[self.go_term_id2go_no[ancestor]], depth)
		#convert ancestor_set to a list to sort
		depth_lca_list = list(ancestor_set)

		if len(depth_lca_list)>0:
			#sort based on depth
			depth_lca_list.sort()
			#find the lca's(might have >1)
			max_depth = depth_lca_list[-1][0]
			if self.debug_common_ancestor_deep_enough:
				print "the maximum depth of the common_ancestors between %s and %s is %d"%\
					(p_go_no, repr(k_functions_set), max_depth)
			#only get those lca's whose depth ==max_depth
			for (depth,lca) in depth_lca_list:
				if depth==max_depth:
					if self.debug_common_ancestor_deep_enough:
						print "lca: %s has the max_depth"%self.go_no2go_id[self.go_term_id2go_no[lca]]
					#map term_id to go_no
					self.lca_list.append(self.go_term_id2go_no[lca])
			return 1
		else:
			if self.debug_common_ancestor_deep_enough:
				print "No common_ancestors deep enough for function %s"%self.go_no2go_id[p_go_no]
			self.lca_list = []
			return 0

	
	def submit(self, curs, gene_table):
		"""
		02-21-05
			Changes to table p_gene,
			1. one row means one gene, one cluster, one function. No merging of the clusters.
			2. avg_p_value is the real p-value.
			3. cluster_context and cluster_array loses its meaning. But I kept cluster_array
				because it's easy. And cluster_context is empty.
				context_specific.py and subgraph_visualize.py are going to be changed.
			4. p_value_cut_off is the real p_value(same as avg_p_value).
			5. recurrence_cut_off is the real recurrence of the cluster.
			6. connectivity_cut_off is the real connectivity of the cluster.
			7. cluster_size_cut_off is the real size of the cluster.
			8. all the predictions are kept, no any cutoff.
			9. add one field, mcl_id to the end of table p_gene to ease table linking.
			10. add two more integers, is_correct_L1, is_correct_lca
			11. unknown_cut_off is the real unknown_gene_ratio
			12. e_accuracy is deprecated
		03-15-05
			add another column lca_list to store the lowest common ancestors between
			the predicted functions and the assigned ones
		03-27-05
			use %s instead of %f when inserting p_value_cut_off and avg_p_value
			'%f' only preserves 6 float digits. '%s' preserves all and converts
			long float to scientific representation.
		"""
		sys.stderr.write("Creating table %s..."%gene_table)
		if gene_table!='p_gene':
			#create the table if it's not 'p_gene'
			curs.execute("create table %s(\
				p_gene_id       serial,\
				gene_no integer,\
				go_no   integer,\
				is_correct      integer,\
				is_correct_L1	integer,\
				is_correct_lca	integer,\
				avg_p_value     float,\
				e_accuracy      float,\
				no_of_clusters  integer,\
				cluster_context varchar,\
				cluster_array   integer[],\
				p_value_cut_off float,\
				recurrence_cut_off      float,\
				connectivity_cut_off    float,\
				cluster_size_cut_off    integer,\
				unknown_cut_off      float,\
				depth_cut_off integer,\
				mcl_id integer,\
				lca_list integer[]\
				)"%gene_table)
		sys.stderr.write("Done.\n")
		
		sys.stderr.write("Database submitting...")
		"""the value of self.prediction_tuple2list, [[p_value, cluster_id, gene_no, go_no, is_correct, \
			is_correct_L1, is_correct_lca, cluster_size, unknown_gene_ratio], [...],  ... ] """
		for (tuple, prediction_list)  in self.prediction_tuple2list.iteritems():
			recurrence = tuple[0]
			connectivity = tuple[1]
			for unit in prediction_list:
				p_value = unit[0]
				mcl_id = unit[1]
				gene_no = unit[2]
				go_no = unit[3]
				is_correct = unit[4]
				is_correct_L1 = unit[5]
				is_correct_lca = unit[6]
				cluster_size = unit[7]
				unknown_gene_ratio = unit[8]
				lca_list = unit[9]
				if len(lca_list)==0:
					curs.execute("insert into %s(gene_no, go_no, is_correct, is_correct_L1, is_correct_lca, \
						avg_p_value, no_of_clusters, cluster_array, p_value_cut_off, recurrence_cut_off,\
						connectivity_cut_off, cluster_size_cut_off, unknown_cut_off, depth_cut_off, mcl_id)\
						values(%d, %d, %d, %d, %d, %s, %s, ARRAY%s, %s, %s, %s, %s, %s, %s, %s)"%\
						(gene_table, gene_no, go_no, is_correct, is_correct_L1, is_correct_lca, \
						p_value, 1, repr([mcl_id]), p_value, recurrence, connectivity, \
						cluster_size, unknown_gene_ratio, self.depth_cut_off, mcl_id))

				else:
					lca_list_string = '{'+repr(lca_list)[1:-1] + '}'
					curs.execute("insert into %s(gene_no, go_no, is_correct, is_correct_L1, is_correct_lca, \
						avg_p_value, no_of_clusters, cluster_array, p_value_cut_off, recurrence_cut_off,\
						connectivity_cut_off, cluster_size_cut_off, unknown_cut_off, depth_cut_off, mcl_id, lca_list)\
						values(%d, %d, %d, %d, %d, %s, %s, ARRAY%s, %s, %s, %s, %s, %s, %s, %s, '%s')"%\
						(gene_table, gene_no, go_no, is_correct, is_correct_L1, is_correct_lca, \
						p_value, 1, repr([mcl_id]), p_value, recurrence, connectivity, \
						cluster_size, unknown_gene_ratio, self.depth_cut_off, mcl_id, lca_list_string))


		sys.stderr.write("done.\n")

	def run(self):
		"""
		03-14-05
			module reuse direction
			load go distance on demand
		"""
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		self.dstruc_loadin(curs)
		if self.dir_files and self.leave_one_out==0:
			sys.stderr.write("working on files of cluster_stat.py results, it must be leave_one_out.\n")
			sys.exit(2)
		if self.dir_files:
			self.core_from_files()
		else:
			self.core(curs)

		if self.needcommit and self.leave_one_out:
			#Database updating is too slow. Do it only if needcommit.
			self.submit(curs, self.gene_table)
		if self.needcommit:
			curs.execute("end")

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema=", "table=", "mcl_table=", \
		"depth_cut_off=", "dir_files=", "leave_one_out", "wu", "report", "commit", "gene_table=", \
		"subgraph_cut_off=", "debug", "recurrence_gap_size=", "connectivity_gap_size="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:m:e:f:x:y:lwrcg:q:b", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = 'cluster_stat'
	mcl_table = 'mcl_result'
	depth_cut_off = 3
	dir_files = None
	leave_one_out = 0
	wu = 0
	report = 0
	commit = 0
	gene_table = 'p_gene'
	subgraph_cut_off = 0
	debug = 0
	recurrence_gap_size = 2
	connectivity_gap_size = 2
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
		elif opt in ("-e", "--depth_cut_off"):
			depth_cut_off = int(arg)
		elif opt in ("-f", "--dir_files"):
			dir_files = arg
		elif opt in ("-l", "--leave_one_out"):
			leave_one_out = 1
		elif opt in ("-w", "--wu"):
			wu = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-g", "--gene_table"):
			gene_table = arg
		elif opt in ("-q", "--subgraph_cut_off="):
			subgraph_cut_off = float(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-x", "--recurrence_gap_size"):
			recurrence_gap_size = int(arg)
		elif opt in ("-y", "--connectivity_gap_size"):
			connectivity_gap_size = int(arg)

	if schema:
		instance = gene_stat(hostname, dbname, schema, table, mcl_table, \
			leave_one_out, wu, report, depth_cut_off, dir_files, commit, gene_table, \
			subgraph_cut_off, debug, recurrence_gap_size, connectivity_gap_size)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
