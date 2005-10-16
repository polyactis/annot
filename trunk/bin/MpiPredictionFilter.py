#!/usr/bin/env mpipython
"""
Usage: MpiPredictionFilter.py -k SCHEMA -i OLD_INPUT_FNAME -j NEW_INPUT_FNAME [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-i ...,	the old input_fname
	-j ...,	the new input_fname
	-m ...,	max cluster size, 100000 (default, no cut)
	-u ...,	unknown_gene_ratio, 1(default, no cut)
	-p ...,	p_value_cut_off, 1(default, no cut)
	-y ...,	is_correct type (2 lca, default)
	-a ...,	cluster accuracy cutoff, (0 default, no cut)
	-e ...,	exponent of the layer, 1(default)
	-s ...,	score for genes, 1,-1,0(default)
	-l ...,	maximum layer, 0(default, no gradient stuff, -e, -s, -l no effect)
	-x ...,	recurrence_x, either to do cutoff, or exponent, 0.8(default)
	-w ...,	recurrence_x's type, 0(nothing, default), 1(cutoff), 2(exponent)
	-v ...,	the size of each passing block, 2000(default)
	-c,	commit the database transaction
	-b,	debug version.
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	mpirun -np 4 -machinefile ~/hostfile /usr/bin/mpipython ~/script/annot/bin/MpiPredictionFilter.py
		-k hs_fim_40 -i hs_fim_40m4x40 -j hs_fim_40m4x40mstest -c -w 1
	
Description:
	Program to filter predictions to form another setting. It creates views for
	splat_table, mcl_table and pattern_table(possibly). The p_gene_table is
	real.
	10-09-05, score is a list, 1,-1,0 corresponds to gene with this function,
	gene without this function, unknown gene

"""

import sys, os, getopt, cPickle	#10-15-05 cPickle is used to serialize objects to pass between nodes
sys.path += [os.path.expanduser('~/script/annot/bin')]
from codense.common import db_connect, form_schema_tables, get_gene_no2go_no_set
from Scientific import MPI
from codense.common import mpi_synchronize
from gene_stat import gene_stat
from sets import Set
import math

class prediction_attributes:
	"""
	10-15-05
		add type
	"""
	def __init__(self, row, type=1):
		self.p_gene_id = row[0]
		self.gene_no = row[1]
		self.go_no = row[2]
		self.is_correct = row[3]
		self.is_correct_l1 = row[4]
		self.is_correct_lca = row[5]
		self.avg_p_value = row[6]
		self.no_of_clusters = row[7]
		self.cluster_array = row[8]
		self.p_value_cut_off = row[9]
		self.recurrence_cut_off = row[10]
		self.connectivity_cut_off = row[11]
		self.cluster_size_cut_off = row[12]
		self.unknown_cut_off = row[13]
		self.depth_cut_off = row[14]
		self.mcl_id = row[15]
		self.lca_list = row[16]
		if type==1:
			self.vertex_set = row[17]	#vertex_set and edge_set is used by SettingCmp.py
			self.edge_set = row[18]
			self.d_matrix = row[19]
		elif type==2:
			self.vertex_gradient = row[17]	#for output_node()
			self.edge_gradient = row[18]
		
		self.is_correct_dict = {0:self.is_correct, 1:self.is_correct_l1, 2:self.is_correct_lca}

class gradient_class:
	"""
	10-10-05
		class to calculate the gradients
	10-15-05
		change the interface to allow more efficient calling
	"""
	def __init__(self, gene_no2go, exponent, score_list, max_layer, debug=0):
		self.gene_no2go = gene_no2go
		self.exponent = float(exponent)
		self.score_list = score_list
		self.max_layer = int(max_layer)
		self.debug = int(debug)

	def get_vertex2no(self, vertex_set):
		vertex2no = {}
		for i in range(len(vertex_set)):
			vertex2no[vertex_set[i]] = i
		return vertex2no
	
	def get_d_row(self, d_matrix, index):
		"""
		10-09-05
			just the row corresponds to the gene_no
		"""
		d_matrix = d_matrix[2:-2].split('},{')
		d_matrix[index] = d_matrix[index].split(',')
		d_matrix[index] = map(int, d_matrix[index])
		return d_matrix[index]
	
	def return_score(self, go_no, v_go_no_set, score_list):
		if go_no in v_go_no_set:
			score = score_list[0]
		elif 0 in v_go_no_set:
			score = score_list[-1]
		else:
			score = score_list[1]
		return score
	
	def cal_vertex_gradient(self, gene_no, go_no, vertex_set, vertex2no, d_row, gene_no2go, exponent, score_list, max_layer):
		vertex_gradient = 0.0
		if self.debug:
			print "gene_no and go_no:",gene_no, go_no, vertex_set, d_row
		for vertex in vertex_set:
			v_no = vertex2no[vertex]
			layer = d_row[v_no]
			if layer<=max_layer and layer>=1:
				v_go_no_set = gene_no2go[vertex]
				score = self.return_score(go_no, v_go_no_set, score_list)
				vertex_gradient += 1/(math.pow(layer, exponent)*math.pow(layer, exponent))*score
				if self.debug:
					print "vertex, v_go_no, score, vertex_gradient:", vertex, v_go_no_set, score, vertex_gradient
		return vertex_gradient
	
	def cal_edge_gradient(self, gene_no, go_no, edge_set_string, vertex2no, d_row, gene_no2go, exponent, score_list, max_layer):
		"""
		10-12-05
		edge score changed from multiply to plus, like score1+score2
		"""
		edge_gradient = 0.0
		edge_set = edge_set_string[2:-2].split('},{')
		if self.debug:
			print "gene_no and go_no:",gene_no, go_no, edge_set, d_row
		for edge in edge_set:
			edge = edge.split(',')
			edge = map(int, edge)
			v1_no = vertex2no[edge[0]]
			v2_no = vertex2no[edge[1]]
			v1_layer = d_row[v1_no]
			v2_layer = d_row[v2_no]
			if v1_layer>0 and v1_layer<=max_layer and v2_layer>0 and v2_layer<=max_layer:
				v1_go_no_set = gene_no2go[edge[0]]
				v2_go_no_set = gene_no2go[edge[1]]
				score1 = self.return_score(go_no, v1_go_no_set, score_list)
				score2 = self.return_score(go_no, v2_go_no_set, score_list)
				edge_gradient += 1/(math.pow(v1_layer, exponent)*math.pow(v2_layer, exponent))*(score1+score2)
				if self.debug:
					print "edge,v1_go_no, v2_go_no, score1, score2, edge_gradient:",edge,v1_go_no_set, v2_go_no_set, score1, score2, edge_gradient
		edge_gradient /= float(len(vertex2no))*(len(vertex2no)-1)/2
		if self.debug:
			print "edge_gradient:", edge_gradient
		return edge_gradient
	
	def cal_gradient(self, gene_no, go_no, vertex_set_string, edge_set_string, d_matrix_string):
		vertex_set = vertex_set_string[1:-1].split(',')
		vertex_set = map(int, vertex_set)
		vertex2no = self.get_vertex2no(vertex_set)
		d_row = self.get_d_row(d_matrix_string, vertex2no[gene_no])
		vertex_gradient = self.cal_vertex_gradient(gene_no,go_no, vertex_set, vertex2no, d_row, \
			self.gene_no2go, self.exponent, self.score_list, self.max_layer)
		edge_gradient = self.cal_edge_gradient(gene_no, go_no, edge_set_string, vertex2no, d_row, \
			self.gene_no2go, self.exponent, self.score_list, self.max_layer)
		return (vertex_gradient, edge_gradient)


class MpiPredictionFilter:
	"""
	09-26-05
		upgrade the function-form to class, 
		p_gene_view is not a view anymore, a real table derived from p_gene_table
			because runtime select cluster_size_cut_off<=max_size blows the memory
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, input_fname=None,\
		jnput_fname=None, max_size=100000, unknown_gene_ratio=1, p_value_cut_off=1,\
		is_correct_type=2, acc_cut_off=0, exponent=1, score_list='1,-1,0', max_layer=0, recurrence_x=0.8,\
		recurrence_x_type=0, size=2000, commit=0, debug=0, report=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_fname = input_fname
		self.jnput_fname = jnput_fname
		self.max_size = int(max_size)
		self.unknown_gene_ratio = float(unknown_gene_ratio)
		self.p_value_cut_off = float(p_value_cut_off)
		self.is_correct_type = int(is_correct_type)
		self.acc_cut_off = float(acc_cut_off)
		self.exponent = float(exponent)
		score_list = score_list.split(',')
		self.score_list = map(float, score_list)
		self.max_layer = int(max_layer)
		self.recurrence_x = float(recurrence_x)
		self.recurrence_x_type = int(recurrence_x_type)
		self.size = int(size)
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)

	
	def get_mcl_id2accuracy(self, curs, p_gene_table, crs_sentence, is_correct_type):
		"""
		10-05-05
		10-12-05
			correct a bug in accuracy calculation
		"""
		sys.stderr.write("Getting mcl_id2accuracy...\n")
		mcl_id2accuracy = {}
		curs.execute("%s"%crs_sentence)
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		counter = 0
		while rows:
			for row in rows:
				prediction_attr_instance = prediction_attributes(row)
				if self.is_good_prediction(prediction_attr_instance):
					if prediction_attr_instance.mcl_id not in mcl_id2accuracy:
						mcl_id2accuracy[prediction_attr_instance.mcl_id]  = []
					mcl_id2accuracy[prediction_attr_instance.mcl_id].append(prediction_attr_instance.is_correct_dict[is_correct_type])
				
				counter += 1
			if self.report:
				sys.stderr.write("%s%s"%('\x08'*20, counter))
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		curs.execute("close crs")
		for mcl_id, is_correct_ls in mcl_id2accuracy.iteritems():
			only_one_ls = map((lambda x: int(x>=1)), is_correct_ls)	#10-12-05 only correct(1) predictions are counted as 1
			one_or_zero_ls = map((lambda x: int(x>=0)), is_correct_ls)	#10-12-05 only correct(1) or wrong(0) are counted as 1, -1 is discarded
			accuracy = sum(only_one_ls)/float(sum(one_or_zero_ls))	#10-12-05 sum the one_or_zero_ls
			mcl_id2accuracy[mcl_id] = accuracy
		sys.stderr.write(" %s clusters. Done.\n"%len(mcl_id2accuracy))
		return mcl_id2accuracy
	
	def fetch_predictions(self, curs, size):
		"""
		10-15-05
		"""
		if self.report:
			sys.stderr.write("Fetching predictions...\n")
		curs.execute("fetch %s from crs"%size)
		rows = curs.fetchall()
		prediction_ls = []
		for row in rows:
			prediction_ls.append(list(row))
		if self.report:
			sys.stderr.write("Fetching done.\n")
		return prediction_ls
	
	def input_node(self, communicator, curs, schema_instance, crs_sentence, size):
		"""
		10-15-05
		10-16-05
			ask output_node for a free_computing_node
		"""
		sys.stderr.write("Reading predictions...\n")
		node_rank = communicator.rank
		curs.execute("%s"%(crs_sentence))
		prediction_ls = self.fetch_predictions(curs, size)
		counter = 0
		while prediction_ls:
			communicator.send("1", communicator.size-1, 1)	#WATCH: tag is 1, to the output_node.
			free_computing_node, source, tag = communicator.receiveString(communicator.size-1, 2)
				#WATCH: tag is 2, from the output_node
			prediction_ls_pickle = cPickle.dumps(prediction_ls, -1)
			communicator.send(prediction_ls_pickle, int(free_computing_node), 0)	#WATCH: int()
			if self.debug:
				sys.stderr.write("block %s sent to %s.\n"%(counter, node_to_receive_block))
			prediction_ls = self.fetch_predictions(curs, size)
			counter += 1
		#tell computing_node to exit the loop
		for node in range(1, communicator.size-1):	#send it to the computing_node
			communicator.send("-1", node, 0)
		sys.stderr.write("Node no.%s reading done\n"%(node_rank))	

	
	def is_good_prediction(self, prediction_attr_instance, mcl_id2accuracy=None):
		"""
		10-05-05
			judge whether this prediction should be taken into consideration
		10-06-05
			add judgement by acc_cut_off
		"""
		if prediction_attr_instance.p_value_cut_off>self.p_value_cut_off or \
			prediction_attr_instance.cluster_size_cut_off>self.max_size or \
			prediction_attr_instance.unknown_cut_off>self.unknown_gene_ratio:
			return 0
		else:
			if mcl_id2accuracy:
				if mcl_id2accuracy[prediction_attr_instance.mcl_id]>=self.acc_cut_off:
					return 1
				else:
					return 0
			else:
				return 1
	
	def node_fire(self, communicator, data, gene_no2go, exponent, score_list, max_layer, mcl_id2accuracy, \
		acc_cut_off, gradient_class_instance, functor):
		"""
		10-15-05
			passing row  should be faster than a p_attr_instance
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		prediction_ls = cPickle.loads(data)
		counter = 0
		no_of_good_predictions = 0
		good_prediction_ls = []
		for row in prediction_ls:
			p_attr_instance = prediction_attributes(row)
			if functor:	#if functor is None, keep the old recurrence_cut_off
				recurrence_array = row[-1][1:-1].split(',')
				recurrence_array = map(float, recurrence_array)
				recurrence_array = map(functor, recurrence_array)
				row[10] = sum(recurrence_array)	#replace the old recurrence_cut_off
			
			if self.is_good_prediction(p_attr_instance, mcl_id2accuracy):
				no_of_good_predictions += 1
				row = row[:-2]	#discard d_matrix_string and recurrence_array
				if max_layer:	#not 0
					row[17], row[18] = gradient_class_instance.cal_gradient(\
						p_attr_instance.gene_no, p_attr_instance.go_no, p_attr_instance.vertex_set,\
						p_attr_instance.edge_set, p_attr_instance.d_matrix)
				else:
					row[17], row[18] = 0,0	#vertex_gradient and edge_gradient equal to 0
				good_prediction_ls.append(row)
			counter += 1
			if self.report:
				sys.stderr.write("%s%s:%s"%('\x08'*20, counter, no_of_good_predictions))
		sys.stderr.write("Node no.%s done with %s good predictions.\n"%(node_rank, no_of_good_predictions))
		return good_prediction_ls
		
	def computing_node(self, communicator,  gene_no2go, exponent, score_list, \
				max_layer, mcl_id2accuracy, acc_cut_off, functor):
		"""
		10-15-05
			
		"""
		node_rank = communicator.rank
		data, source, tag = communicator.receiveString(0, 0)	#get data from node 0
		gradient_class_instance = gradient_class(gene_no2go, exponent, score_list, max_layer, self.debug)
		while 1:
			if data=="-1":
				if self.debug:
					sys.stderr.write("node %s breaked.\n"%node_rank)
				break
			else:
				prediction_ls = self.node_fire(communicator, data, gene_no2go, exponent, score_list, max_layer, \
					mcl_id2accuracy, acc_cut_off, gradient_class_instance, functor)
				prediction_ls_pickle = cPickle.dumps(prediction_ls, -1)
				communicator.send(prediction_ls_pickle, communicator.size-1, 1)
			data, source, tag = communicator.receiveString(0, 0)	#get data from node 0
		#tell the last node to stop
		communicator.send("-1", communicator.size-1, 1)
	
	def view_from_table(self, curs, table, view):
		curs.execute("CREATE OR REPLACE VIEW %s AS SELECT * FROM %s"%(view, table))
	
	def createGeneTable(self, curs, gene_table):
		"""
		10-09-05
			copied from gene_stat.py
			add vertex_gradient and edge_gradient
		"""
		sys.stderr.write("Creating table %s..."%gene_table)
		if gene_table!='p_gene':
			#create the table if it's not 'p_gene'
			curs.execute("create table %s(\
				p_gene_id       serial primary key,\
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
				lca_list integer[],\
				vertex_gradient	float,\
				edge_gradient	float\
				)"%gene_table)
		sys.stderr.write("Done.\n")	

	def submit_to_p_gene_table(self, curs, p_gene_table, p_attr_instance):
		"""
		10-05-05
		10-09-05
			add vertex_gradient and edge_gradient
		"""
		if p_attr_instance.lca_list:
			curs.execute("insert into %s(p_gene_id, gene_no, go_no, is_correct, is_correct_l1, \
			is_correct_lca, avg_p_value, no_of_clusters, cluster_array, p_value_cut_off, recurrence_cut_off, \
			connectivity_cut_off, cluster_size_cut_off, unknown_cut_off, depth_cut_off, mcl_id, lca_list, \
			vertex_gradient, edge_gradient)\
			values(%s, %s, %s, %s, %s,\
			%s, %s, %s, '%s', %s, %s,\
			%s, %s, %s, %s, %s, '%s', \
			%s, %s)"%\
			(p_gene_table, \
			p_attr_instance.p_gene_id, p_attr_instance.gene_no, p_attr_instance.go_no, p_attr_instance.is_correct, p_attr_instance.is_correct_l1,\
			p_attr_instance.is_correct_lca, p_attr_instance.avg_p_value, p_attr_instance.no_of_clusters, p_attr_instance.cluster_array, p_attr_instance.p_value_cut_off, p_attr_instance.recurrence_cut_off,\
			p_attr_instance.connectivity_cut_off, p_attr_instance.cluster_size_cut_off, p_attr_instance.unknown_cut_off, p_attr_instance.depth_cut_off, p_attr_instance.mcl_id, p_attr_instance.lca_list,\
			p_attr_instance.vertex_gradient, p_attr_instance.edge_gradient))
		else:
			curs.execute("insert into %s(p_gene_id, gene_no, go_no, is_correct, is_correct_l1, \
			is_correct_lca, avg_p_value, no_of_clusters, cluster_array, p_value_cut_off, recurrence_cut_off, \
			connectivity_cut_off, cluster_size_cut_off, unknown_cut_off, depth_cut_off, mcl_id,\
			vertex_gradient, edge_gradient)\
			values(%s, %s, %s, %s, %s,\
			%s, %s, %s, '%s', %s, %s,\
			%s, %s, %s, %s, %s,\
			%s, %s)"%\
			(p_gene_table, \
			p_attr_instance.p_gene_id, p_attr_instance.gene_no, p_attr_instance.go_no, p_attr_instance.is_correct, p_attr_instance.is_correct_l1,\
			p_attr_instance.is_correct_lca, p_attr_instance.avg_p_value, p_attr_instance.no_of_clusters, p_attr_instance.cluster_array, p_attr_instance.p_value_cut_off, p_attr_instance.recurrence_cut_off,\
			p_attr_instance.connectivity_cut_off, p_attr_instance.cluster_size_cut_off, p_attr_instance.unknown_cut_off, p_attr_instance.depth_cut_off, p_attr_instance.mcl_id,\
			p_attr_instance.vertex_gradient, p_attr_instance.edge_gradient))
	
	def output_node(self, communicator, curs, output_table):
		"""
		10-15-05
		10-16-05
			reserve a pool of free_computing_nodes for input_node to choose
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s ready to accept output...\n"%node_rank)
		data, source, tag = communicator.receiveString(None, 1)
		no_of_resting_nodes = 0	#to keep track how many computing_nodes have rested
		free_computing_nodes = range(1,communicator.size-1)
		while 1:
			if data=="-1":
				no_of_resting_nodes += 1
				if self.debug:
					sys.stderr.write("node %s(%s-th) rested.\n"%(source, no_of_resting_nodes))
				if no_of_resting_nodes==communicator.size-2:	#all computing_nodes have stopped, i'm done.
					break
					if self.debug:
						sys.stderr.write("node %s output finished.\n"%node_rank)
			elif data=="1":	#the input_node is asking me for free computing_node
				free_computing_node = free_computing_nodes.pop(0)
				communicator.send(str(free_computing_node), source, 2)	#WATCH tag is 2.
			else:
				free_computing_nodes.append(source)	#append the free computing_node
				prediction_ls = cPickle.loads(data)
				for row in prediction_ls:
					p_attr_instance = prediction_attributes(row, type=2)
					self.submit_to_p_gene_table(curs, output_table, p_attr_instance)
			data, source, tag = communicator.receiveString(None, 1)
		sys.stderr.write("Node no.%s output done.\n"%node_rank)
	
	def run(self):
		"""
		10-05-05
		10-12-05
			use max_layer to control whether to turn on the gradient or not
		10-16-05
			transformed to MPI version
		
			if node_rank==0
				--db_connect()
				--form_schema_tables()
				--form_schema_tables()
				--get_gene_no2go_no_set()
				--get_mcl_id2accuracy()
			elif computing_node:
				(prepare data)
			elif output_node:
				--db_connect()
				--form_schema_tables()
				--form_schema_tables()
				--view_from_table()
				--view_from_table()
				--view_from_table()
				--createGeneTable()
			
			--mpi_synchronize()
			
			if input_node:
				--input_node()
					--fetch_predictions()
			elif computing_node:
				--computing_node()
					--node_fire()
						--gradient_class()
			elif output_node:
				--output_node()
					--submit_to_p_gene_table()
		"""		
		communicator = MPI.world.duplicate()
		node_rank = communicator.rank
		if node_rank == 0:
			(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
			old_schema_instance = form_schema_tables(self.input_fname)
			new_schema_instance = form_schema_tables(self.jnput_fname)
			gene_no2go = get_gene_no2go_no_set(curs)
			gene_no2go_pickle = cPickle.dumps(gene_no2go, -1)	#-1 means use the highest protocol
			if self.acc_cut_off:
				mcl_id2accuracy = self.get_mcl_id2accuracy(curs, old_schema_instance.p_gene_table, crs_sentence, self.is_correct_type)
			else:
				mcl_id2accuracy = None
			mcl_id2accuracy_pickle = cPickle.dumps(mcl_id2accuracy, -1)	#-1 means use the highest protocol
			for node in range(1, communicator.size-1):	#send it to the computing_node
				communicator.send(gene_no2go_pickle, node, 0)
			for node in range(1, communicator.size-1):	#send it to the computing_node
				communicator.send(mcl_id2accuracy_pickle, node, 0)
		elif node_rank<=communicator.size-2:	#exclude the last node
			data, source, tag = communicator.receiveString(0, 0)
			gene_no2go = cPickle.loads(data)	#take the data
			data, source, tag = communicator.receiveString(0, 0)
			mcl_id2accuracy = cPickle.loads(data)	#take the data
			#choose a functor for recurrence_array
			functor_dict = {0: None,
				1: lambda x: int(x>=self.recurrence_x),
				2: lambda x: math.pow(x, self.recurrence_x)}
			functor = functor_dict[self.recurrence_x_type]
		elif node_rank==communicator.size-1:
			(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
			old_schema_instance = form_schema_tables(self.input_fname)
			new_schema_instance = form_schema_tables(self.jnput_fname)
			self.view_from_table(curs, old_schema_instance.splat_table, new_schema_instance.splat_table)
			self.view_from_table(curs, old_schema_instance.mcl_table, new_schema_instance.mcl_table)
			self.view_from_table(curs, old_schema_instance.pattern_table, new_schema_instance.pattern_table)
			self.createGeneTable(curs, new_schema_instance.p_gene_table)
		
		mpi_synchronize(communicator)
		
		if node_rank == 0:
			if self.max_layer:	#not 0
				crs_sentence = 'DECLARE crs CURSOR FOR SELECT p.p_gene_id, p.gene_no, p.go_no, p.is_correct, p.is_correct_l1, \
				p.is_correct_lca, p.avg_p_value, p.no_of_clusters, p.cluster_array, p.p_value_cut_off, p.recurrence_cut_off, \
				p.connectivity_cut_off, p.cluster_size_cut_off, p.unknown_cut_off, p.depth_cut_off, p.mcl_id, p.lca_list, \
				d.vertex_set, s.edge_set, d.d_matrix, m.recurrence_array from %s p, %s s, %s d, %s m where \
				p.mcl_id=s.splat_id and p.mcl_id=d.splat_id and p.mcl_id=m.mcl_id'%(old_schema_instance.p_gene_table, \
				old_schema_instance.splat_table, old_schema_instance.d_matrix_table, old_schema_instance.mcl_table)
			else:
				crs_sentence = "DECLARE crs CURSOR FOR SELECT p.p_gene_id, p.gene_no, p.go_no, p.is_correct, p.is_correct_l1, \
				p.is_correct_lca, p.avg_p_value, p.no_of_clusters, p.cluster_array, p.p_value_cut_off, p.recurrence_cut_off, \
				p.connectivity_cut_off, p.cluster_size_cut_off, p.unknown_cut_off, p.depth_cut_off, p.mcl_id, p.lca_list, 'vertex_set',\
				'edge_set', 'd_matrix', m.recurrence_array\
				from %s p, %s m where p.mcl_id=m.mcl_id"%(old_schema_instance.p_gene_table, old_schema_instance.mcl_table)
				
				#some placeholders 'vertex_set', 'edge_set', 'd_matrix' for prediction_attributes()
			self.input_node(communicator, curs, old_schema_instance, crs_sentence, self.size)
		elif node_rank<=communicator.size-2:	#exclude the last node
			self.computing_node(communicator, gene_no2go, self.exponent, self.score_list, \
				self.max_layer, mcl_id2accuracy, self.acc_cut_off, functor)
		elif node_rank==communicator.size-1:
			self.output_node(communicator, curs, new_schema_instance.p_gene_table)
			if self.commit:
				curs.execute("end")

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:j:m:u:p:y:a:e:s:l:x:w:v:cbr", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	input_fname = None
	jnput_fname = None
	max_size = 100000
	unknown_gene_ratio = 1
	p_value_cut_off = 1
	is_correct_type = 2
	acc_cut_off = 0
	exponent = 1
	score_list = '1,-1,0'
	max_layer = 0
	recurrence_x = 0.8
	recurrence_x_type = 0
	size = 2000
	commit = 0
	debug = 0
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
		elif opt in ("-i"):
			input_fname = arg
		elif opt in ("-j"):
			jnput_fname = arg
		elif opt in ("-m"):
			max_size = int(arg)
		elif opt in ("-u"):
			unknown_gene_ratio = float(arg)
		elif opt in ("-p"):
			p_value_cut_off = float(arg)
		elif opt in ("-y"):
			is_correct_type = int(arg)
		elif opt in ("-a"):
			acc_cut_off = float(arg)
		elif opt in ("-e"):
			exponent = float(arg)
		elif opt in ("-s"):
			score_list = arg
		elif opt in ("-l"):
			max_layer = int(arg)
		elif opt in ("-x"):
			recurrence_x = float(arg)
		elif opt in ("-w"):
			recurrence_x_type = int(arg)
		elif opt in ("-v"):
			size = int(arg)
		elif opt in ("-c"):
			commit = 1
		elif opt in ("-b"):
			debug = 1
		elif opt in ("-r"):
			report = 1
	if schema and input_fname and jnput_fname:
		instance = MpiPredictionFilter(hostname, dbname, schema, input_fname, jnput_fname, max_size, \
			unknown_gene_ratio, p_value_cut_off, is_correct_type, acc_cut_off, exponent, score_list, \
			max_layer, recurrence_x, recurrence_x_type, size, commit, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
