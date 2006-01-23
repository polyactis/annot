#!/usr/bin/env mpipython
"""
Usage: MpiStatCluster.py -k SCHEMA -i -j [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-i ...,	input_fname, from MpiBFSCluster.py
	-j ..., jnput_fname(output file)
	--de=...,	the minimum GO tree depth(root is 0) to consider, 5(default)
	-u ...,	min #genes associated with the candidate go in layer 1, 2(default)
	-g ...,	min ratio of #associated genes on layer 1, 0.5(default)
	-o ...,	min #associated genes on layer 2, 2(default)
	-p ...,	min ratio of #associated genes on layer 2, 0.1(default)
	-e ...,	exponent of the layer, 2(default)
	-s ...,	score for genes, 3,-1,0(default)
	-l ...,	maximum layer, 5(default)
	-q ...,	exponent of the normalization factor in edge-gradient's denominator, 0.8(default)
	-t ...,	edge_gradient type, 1 (edge_gradient), 0(edge gradient_score),2(vertex_gradient, default)
	-x ...,	recurrence_x, either to do cutoff, or exponent, 5(default)	(IGNORE)
	-w ...,	recurrence_x's type, 0(nothing), 1(cutoff), 2(exponent, default)	(IGNORE)
	-v ...,	the size of message(by string_length of each prediction), 10000000(default)
	-f ...,	file contains the go_no2edge_counter_list data structure
	-n,	need to createGeneTable()	(IGNORE)
	-c,	commit the database transaction (IGNORE)
	-b,	debug version.
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	mpirun -np 20 -machinefile ~/hostfile /usr/bin/mpipython ~/script/annot/bin/MpiStatCluster.py
	-k hs_fim_40 -i -j -r
	
Description:
	Program to do function prediction statistics based on gradient.
	01-02-06 program starts to work on plain files
	output format:
	id, gene_no, go_no, GradientScorePrediction_instance.depth, gradient_score, edge_gradient, \
	is_correct, is_correct_L1, is_correct_lca, one_dim_list2string(lca_list)
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import sys, os, getopt, csv, math, Numeric, cPickle
from Scientific import MPI
from codense.common import mpi_synchronize, db_connect, get_gene_no2go_no_set, get_go_no2depth,\
	get_go_no2edge_counter_list, combine_numerator_a_denominator_dict, one_dim_list2string, \
	get_no_of_unknown_genes, get_go_no2gene_no_set, form_schema_tables, output_node, input_node,\
	computing_node
from sets import Set
from MpiPredictionFilter import prediction_attributes, gradient_class, MpiPredictionFilter
from gene_stat import gene_stat
from heapq import heappush, heappop
if sys.version_info[:2] < (2, 3):       #python2.2 or lower needs some extra
	from python2_3 import *
	
class MpiSetupGo_no2Edges:
	"""
	10-20-05
		not used because cPickle approach is easier.
	"""
	def __init__(self,communicator=None, hostname='zhoudb', dbname='graphdb', schema=None, \
		edge_table='edge_cor_vector', gene_no2go_no_set=None, debug=0, report=0):
		self.communicator = communicator
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.edge_table = edge_table
		self.gene_no2go_no_set = gene_no2go_no_set
		self.debug = int(debug)
		self.report = int(report)
		
		self.go_no2edge_counter_list = {}	#value is a six integer list, each corresponds to below edge_type
		self.edge_type2index = {(0,0):0, (0,1):1, (1,0):1, (0,2):2, (2,0):2, (1,1):3, (1,2):4, (2,1):4, (2,2):5 }
			#this is used to map to the index of the value of go_no2edge_counter_list,
			#1st digit is gene1, 2nd digit is gene2,(order doesn't matter), 0 means gene having this function,
			#1 means gene having the other function, 2 means unknown gene
	
	def fetch_edge_block(self, curs, size):
		"""
		10-10-05
			format:
			[
			[g1, g2],
			[g3, g4],
			...
			]
		"""
		curs.execute("fetch %s from crs"%(size))
		rows = curs.fetchall()
		block = []
		for row in rows:
			edge = row[0][1:-1].split(',')
			edge = map(int, edge)
			block.append(edge)
		block = Numeric.array(block)
		return block
	
	def input_node(self, communicator, curs, edge_table, size=10000):
		sys.stderr.write("Reading clusters from tables...\n")
		node_rank = communicator.rank
		curs.execute("DECLARE crs CURSOR FOR select edge_name\
			from %s "%(edge_table))
		cluster_block = self.fetch_edge_block(curs, size)
		which_block = 0
		while cluster_block:
			for node in range(1, communicator.size-1):	#send it to the computing_node
				communicator.send(cluster_block, node, 0)
			if self.report:
				sys.stderr.write("block %s sent to everyone.\n"%(which_block))
			cluster_block = self.fetch_edge_block(curs, size)
			which_block += 1
		#tell computing_node to exit the loop
		stop_signal = Numeric.zeros((1,1), Numeric.Int)
		stop_signal[0,0] = -1
		for node in range(1, communicator.size-1):	#send it to the computing_node
			communicator.send(stop_signal, node, 0)
		sys.stderr.write("Node no.%s reading done\n"%(node_rank))
	
	def node_fire(self, communicator, data, gene_no2go_no_set):
		"""
		10-10-05
		"""
		node_rank = communicator.rank
		if self.debug:
			sys.stderr.write("Node no.%s processing...\n"%node_rank)
		edge_type = [0,0]	#to map to the index in the value of go_no2edge_counter_list
		for no in range(0, len(data), 2):
			gene_no1 = data[no]
			gene_no2 = data[no+1]
			go_no_set1 = gene_no2go_no_set[gene_no1]
			go_no_set2 = gene_no2go_no_set[gene_no2]
			go_no_set = go_no_set1 | go_no_set2
			for go_no in go_no_set:
				if go_no not in self.go_no2edge_counter_list:
					self.go_no2edge_counter_list[go_no] = [0]*6
				if go_no==0:
					if go_no in go_no_set1:
						edge_type[0] = 2
					else:
						edge_type[0] = 1
					if go_no in go_no_set2:
						edge_type[1] = 2
					else:
						edge_type[1] = 1
				else:
					if go_no in go_no_set1:
						edge_type[0] = 0
					else:
						edge_type[0] = 1
					if go_no in go_no_set2:
						edge_type[1] = 0
					else:
						edge_type[1] = 1
				self.go_no2edge_counter_list[self.edge_type2index[tuple(edge_type)]] += 1	#see doc in __init__()
				
		if self.debug:
			sys.stderr.write("Node no.%s done.\n"%node_rank)
	
	def computing_node(self, communicator, gene_no2go_no_set):
		"""
		10-10-05
			
		"""
		node_rank = communicator.rank
		data, source, tag, count = communicator.receive(Numeric.Int, 0, 0)	#get data from node 0
		while 1:
			if data[0]==-1:
				if self.debug:
					sys.stderr.write("node %s breaked.\n"%node_rank)
				break
			else:
				self.node_fire(communicator, data, gene_no2go_no_set)
			data, source, tag, count = communicator.receive(Numeric.Int, 0, 0)	#get data from node 0
		
	def run(self):
		"""
		10-10-05
			--db_connect()
			--input_node()
				--fetch_edge_block()
			--computing_node()
				--node_fire()
		"""
		node_rank = self.communicator.rank
		if node_rank == 0:
			(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
			self.input_node(self.communicator, curs, self.edge_table_table)
		elif node_rank<=self.communicator.size-3:	#exclude the last two nodes
			self.computing_node(self.communicator, self.gene_no2go_no_set)
		mpi_synchronize(self.communicator)
		del conn, curs

class GradientScorePrediction(gradient_class):
	"""
	10-16-05
		inherit from gradient_class
		calculate the gradient score for each go_no and return the one with highest gradient-score
		
	"""
	def __init__(self, gene_no2go, go_no2gene_no_set, go_no2depth, go_no2edge_counter_list, no_of_unknown_genes, \
		depth, min_layer1_associated_genes, min_layer1_ratio, min_layer2_associated_genes, min_layer2_ratio, exponent,\
		score_list, max_layer, norm_exp, eg_d_type, debug):
		gradient_class.__init__(self, gene_no2go, exponent, score_list, max_layer, norm_exp, eg_d_type, debug)
		self.go_no2gene_no_set = go_no2gene_no_set
		self.go_no2depth = go_no2depth
		self.go_no2edge_counter_list = go_no2edge_counter_list
		self.no_of_unknown_genes = int(no_of_unknown_genes)
		self.depth = int(depth)
		self.min_layer1_associated_genes = int(min_layer1_associated_genes)
		self.min_layer1_ratio = float(min_layer1_ratio)
		self.min_layer2_associated_genes = int(min_layer2_associated_genes)
		self.min_layer2_ratio = float(min_layer2_ratio)
		
		#get the no_of_total_edges from 0's value
		if go_no2edge_counter_list:	#01-22-06
			self.no_of_total_edges = go_no2edge_counter_list[0][0]+go_no2edge_counter_list[0][1]+go_no2edge_counter_list[0][3]
		else:
			self.no_of_total_edges = None
		self.no_of_total_genes = len(gene_no2go)	#get the no_of_total_genes from gene_no2go
		self.vertex_set = None
		self.edge_set = None
		self.d_matrix = None
		self.vertex2no = None
		self.edge_type2index = {(0,0):0, (0,1):1, (1,0):1, (0,2):2, (2,0):2, (1,1):3, (1,2):4, (2,1):4, (2,2):5 }
			#this is used to map to the index of the value of go_no2edge_counter_list,
			#1st digit is gene1, 2nd digit is gene2,(order doesn't matter), 0 means gene having this function,
			#1 means gene having the other function, 2 means unknown gene
		self.cal_func_dict = {0: self.cal_score_edge_gradient_1,
			1: self.cal_score_edge_gradient_1,
			2: self.cal_score_edge_gradient_2
			}
		self.cal_denominator_dict = {0: self.cal_gradient_denominator_1,
			1: self.cal_gradient_denominator_1,
			2: self.cal_gradient_denominator_2}
	
	def get_d_matrix(self, d_matrix_string):
		"""
		01-02-06
			d_matrix_string is from the output file of MpiBFSCluster.py
		"""
		d_matrix = d_matrix_string[2:-2].split('], [')
		for index in range(len(d_matrix)):
			d_matrix[index] = d_matrix[index].split(',')
			d_matrix[index] = map(int, d_matrix[index])
		return d_matrix
	
	def get_candidate_go_nos(self, gene_no2go, go_no2depth, vertex2no, layer2no_of_vertices, d_row, depth,\
		min_layer1_associated_genes, min_layer1_ratio, min_layer2_associated_genes, min_layer2_ratio):
		"""
		10-20-05
			several criteria
			1.  associated genes on layer 1
			2. >= depth
			3. >=min_layer1_associated_genes associated
		10-25-05 add more requirements
			1. layer1>=2 genes, >=50%
			2. layer2>=2 genes, >=10%
			
			The 0 function shall be excluded because of depth.
		"""
		if self.debug:
			sys.stderr.write("Getting candidate go_no...\n")
		go_no2layer2associated_genes = {}
		for vertex, no in vertex2no.iteritems():
			layer = d_row[no]
			if layer>=1 and layer<=2:	#only layer 1
				for go_no in gene_no2go[vertex]:
					if go_no2depth[go_no]>=depth:	#only go_no under depth_cut_off
						if go_no not in go_no2layer2associated_genes:
							go_no2layer2associated_genes[go_no] = {}
						if layer not in go_no2layer2associated_genes[go_no]:
							go_no2layer2associated_genes[go_no][layer] = 0
						go_no2layer2associated_genes[go_no][layer] += 1
		candidate_go_no_ls = []
		for go_no, layer2associated_genes in go_no2layer2associated_genes.iteritems():
			if 1 in layer2associated_genes:
				layer1_associated_genes = layer2associated_genes[1]
			else:
				layer1_associated_genes = 0.0
			if 2 in layer2associated_genes:
				layer2_associated_genes = layer2associated_genes[2]
			else:
				layer2_associated_genes = 0.0
			if 1 in layer2no_of_vertices:
				layer1_ratio = float(layer1_associated_genes)/layer2no_of_vertices[1]
			else:
				layer1_ratio = 0.0
			if 2 in layer2no_of_vertices:
				layer2_ratio = float(layer2_associated_genes)/layer2no_of_vertices[2]
			else:
				layer2_ratio = 0.0
			if layer1_associated_genes>=min_layer1_associated_genes and layer1_ratio>=min_layer1_ratio \
				and layer2_associated_genes>=min_layer2_associated_genes and layer2_ratio>=min_layer2_ratio:
				candidate_go_no_ls.append(go_no)
		if self.debug:
			sys.stderr.write("candidate go_no done.\n")
		return candidate_go_no_ls
	
	def setup_cluster_info(self, vertex_set, edge_set, d_matrix_string):
		"""
		10-20-05
			Must be run before get_prediction()
			
			--get_d_matrix()
		"""
		self.vertex_set = vertex_set
		self.edge_set = edge_set
		self.d_matrix = self.get_d_matrix(d_matrix_string)
		self.vertex2no = self.get_vertex2no(self.vertex_set)
	
	def cal_layer2no_of_vertices(self, d_row):
		"""
		10-25-05 calculate the layer2no_of_vertices
		"""
		layer2no_of_vertices = {}
		for i in range(len(d_row)):
			layer = d_row[i]
			if layer not in layer2no_of_vertices:
				layer2no_of_vertices[layer] = 0
			layer2no_of_vertices[layer] += 1
		return layer2no_of_vertices
	
	def cal_gradient_denominator_1(self, d_row, exponent, max_layer, norm_exp, layer2no_of_vertices):
		"""
		10-21-05
			return layer2edge_gradient_denominator
		10-25-05 get layer2no_of_vertices directly from arguments
		"""
		layer2gradient_score_denominator = {}
		layer2edge_gradient_denominator = {}
		for i in range(1, max_layer+1):
			if i not in layer2no_of_vertices:
				break
			#1st, possible edges between i-1 and i
			n_i_1 = layer2no_of_vertices[i-1]
			n_i = layer2no_of_vertices[i]
			if i==1:	#10-19-05 layer 0 and 1 uses a different a exponent
				edge_gradient_denominator_index = math.pow(n_i_1*n_i, norm_exp-0.2)
			else:
				edge_gradient_denominator_index = math.pow(n_i_1*n_i, norm_exp)
			key_layer = 2*i-1
			denominator = math.pow(key_layer, exponent)
			layer2edge_gradient_denominator[key_layer] = edge_gradient_denominator_index*denominator
			layer2gradient_score_denominator[key_layer] = denominator
			#2nd possible edges within i
			#edge_gradient_denominator += n_i*(n_i-1)/(4*math.pow(i,exponent))
			key_layer = 2*i
			edge_gradient_denominator_index = math.pow(n_i*(n_i-1)/2, norm_exp)
			denominator = math.pow(key_layer,exponent)
			layer2edge_gradient_denominator[key_layer] = edge_gradient_denominator_index*denominator
			layer2gradient_score_denominator[key_layer] = denominator
		if self.debug:
			print "d_row",d_row, "layer2gradient_score_denominator", layer2gradient_score_denominator, \
				"layer2edge_gradient_denominator", layer2edge_gradient_denominator
		return layer2gradient_score_denominator, layer2edge_gradient_denominator
	
	def cal_gradient_denominator_2(self, d_row, exponent, max_layer, norm_exp, layer2no_of_vertices):
		"""
		10-21-05
			type 2, no consideration between layer i and i-1
		10-25-05 get layer2no_of_vertices directly from arguments
			
			return layer2edge_gradient_denominator
		"""
		layer2gradient_score_denominator = {}
		layer2edge_gradient_denominator = {}
		for i in range(1, max_layer+1):
			if i not in layer2no_of_vertices:
				break
			key_layer = 2*i
			n_i = layer2no_of_vertices[i]
			edge_gradient_denominator_index = math.pow(n_i, norm_exp)
			denominator = math.pow(key_layer,exponent)
			layer2edge_gradient_denominator[key_layer] = edge_gradient_denominator_index*denominator
			layer2gradient_score_denominator[key_layer] = denominator
		if self.debug:
			print "d_row",d_row, "layer2gradient_score_denominator", layer2gradient_score_denominator, \
				"layer2edge_gradient_denominator", layer2edge_gradient_denominator
		return layer2gradient_score_denominator, layer2edge_gradient_denominator
	
	def cal_score_edge_gradient_1(self, gene_no, go_no, edge_set, vertex2no, d_row, gene_no2go, go_no2gene_no_set, exponent, \
		score_list, max_layer, no_of_total_edges, no_of_total_genes, no_of_unknown_genes):
		"""
		10-17-05
			edges between layer 0 and 1, we use vertex-based score or probability
			
			[return] layer2gradient_score, layer2edge_gradient
			[global structures]: go_no2edge_counter_list, edge_type2index
		01-02-06
			math.log is one-argument function (=ln) in python2.2
		"""
		if self.debug:
			print "gene_no",gene_no, "go_no", go_no, "edge_set", edge_set, "d_row", d_row
		layer2edge_gradient = {}
		layer2gradient_score = {}
		for edge in edge_set:
			v1_no = vertex2no[edge[0]]
			v2_no = vertex2no[edge[1]]
			v1_layer = d_row[v1_no]
			v2_layer = d_row[v2_no]
			if v1_layer>0 and v1_layer<=max_layer and v2_layer>0 and v2_layer<=max_layer:
				v1_go_no_set = gene_no2go[edge[0]]
				v2_go_no_set = gene_no2go[edge[1]]
				v1_type = self.return_vertex_type(go_no, v1_go_no_set)
				v2_type = self.return_vertex_type(go_no, v2_go_no_set)
				score1 = score_list[v1_type]
				score2 = score_list[v2_type]
				key_layer = v1_layer + v2_layer
				edge_gradient_part = score1+score2
				if key_layer not in layer2edge_gradient:
					layer2edge_gradient[key_layer] = 0.0
					layer2gradient_score[key_layer] = 0.0
				layer2edge_gradient[key_layer] += edge_gradient_part
				edge_type = (v1_type, v2_type)
				p_e = self.go_no2edge_counter_list[go_no][self.edge_type2index[edge_type]]/float(no_of_total_edges)
				layer2gradient_score[key_layer] += edge_gradient_part*(-math.log(p_e)/math.log(10))	#01-02-06
				if self.debug:
					print "edge,v1_go_no, v2_go_no, score1, score2, edge_type, p_e, layer2gradient_score, layer2edge_gradient",\
					edge,v1_go_no_set, v2_go_no_set, score1, score2, edge_type, p_e,layer2gradient_score, layer2edge_gradient
			elif v1_layer==0 and v2_layer<=max_layer:
				v2_go_no_set = gene_no2go[edge[1]]
				v2_type = self.return_vertex_type(go_no, v2_go_no_set)
				score2 = score_list[v2_type]
				edge_gradient_part = score2
				key_layer = v2_layer*2 -1
				if key_layer not in layer2edge_gradient:
					layer2edge_gradient[key_layer] = 0.0
					layer2gradient_score[key_layer] = 0.0
				layer2edge_gradient[key_layer] += edge_gradient_part
				if v2_type ==0:
					p_e = len(go_no2gene_no_set[go_no])/float(no_of_total_genes)
				elif v2_type == 1:	#genes with other function
					p_e = (no_of_total_genes-len(go_no2gene_no_set[go_no])-no_of_unknown_genes)/float(no_of_total_genes)
				elif v2_type == 2:	#unknown genes
					p_e = no_of_unknown_genes/float(no_of_total_genes)
				layer2gradient_score[key_layer] += edge_gradient_part*(-math.log(p_e)/math.log(10))
				if self.debug:
					print "edge, v2_go_no, score2, p_e, layer2gradient_score, layer2edge_gradient",\
					edge, v2_go_no_set, score2, p_e,layer2gradient_score, layer2edge_gradient
			elif v2_layer==0 and v1_layer<=max_layer:
				v1_go_no_set = gene_no2go[edge[0]]
				v1_type = self.return_vertex_type(go_no, v1_go_no_set)
				score1 = score_list[v1_type]
				edge_gradient_part = score1
				key_layer = v1_layer*2 - 1
				if key_layer not in layer2edge_gradient:
					layer2edge_gradient[key_layer] = 0.0
					layer2gradient_score[key_layer] = 0.0
				layer2edge_gradient[key_layer] += edge_gradient_part
				if v1_type ==0:
					p_e = len(go_no2gene_no_set[go_no])/float(no_of_total_genes)
				elif v1_type == 1:	#genes with other function
					p_e = (no_of_total_genes-len(go_no2gene_no_set[go_no])-no_of_unknown_genes)/float(no_of_total_genes)
				elif v1_type == 2:	#unknown genes
					p_e = no_of_unknown_genes/float(no_of_total_genes)
				layer2gradient_score[key_layer] += edge_gradient_part*(-math.log(p_e)/math.log(10))
				if self.debug:
					print "edge, v1_go_no, score1, p_e,layer2gradient_score, layer2edge_gradient",\
					edge, v1_go_no_set, score1, p_e,layer2gradient_score, layer2edge_gradient
		return layer2gradient_score, layer2edge_gradient
	
	def cal_score_edge_gradient_2(self, gene_no, go_no, edge_set, vertex2no, d_row, gene_no2go, go_no2gene_no_set, exponent, \
		score_list, max_layer, no_of_total_edges, no_of_total_genes, no_of_unknown_genes):
		"""
		10-20-05 type 2: vertex based gradient
			
			[return]  layer2gradient_score, layer2edge_gradient
		01-02-06
			math.log is one-argument function (=ln) in python2.2
		"""
		if self.debug:
			print "gene_no",gene_no, "go_no", go_no, "edge_set", edge_set, "d_row", d_row
		layer2edge_gradient = {}
		layer2gradient_score = {}
		for vertex in vertex2no:
			v_no = vertex2no[vertex]
			layer = d_row[v_no]
			if layer<=max_layer and layer>=1:
				v_go_no_set = gene_no2go[vertex]
				v_type = self.return_vertex_type(go_no, v_go_no_set)
				score = score_list[v_type]
				key_layer = layer*2
				if key_layer not in layer2edge_gradient:
					layer2edge_gradient[key_layer] = 0.0
					layer2gradient_score[key_layer] = 0.0
				layer2edge_gradient[key_layer] += score
				if v_type ==0:
					p_e = len(go_no2gene_no_set[go_no])/float(no_of_total_genes)
				elif v_type == 1:	#genes with other function
					p_e = (no_of_total_genes-len(go_no2gene_no_set[go_no])-no_of_unknown_genes)/float(no_of_total_genes)
				elif v_type == 2:	#unknown genes
					p_e = no_of_unknown_genes/float(no_of_total_genes)
				layer2gradient_score[key_layer] += score*(-math.log(p_e)/math.log(10))
				if self.debug:
					print "vertex, v_go_no, score, p_e, layer2gradient_score, layer2edge_gradient", \
						vertex, v_go_no_set, score, p_e, layer2gradient_score, layer2edge_gradient
		return layer2gradient_score, layer2edge_gradient
		
	def get_prediction(self, gene_no):
		"""
		10-20-05
			--setup_cluster_info() must be run ahead.
		10-22-05
			use heap to keep track of predictions
			no-prob to probability conversion degree changed to 10
		10-23-05 type 0 refers to edge gradient_score
			
			return list of [gradient_score, edge_gradient, go_no]
			
			--cal_denominator_dict()
			--get_candidate_go_nos()
			(loop)
				--cal_func_dict()
					--return_vertex_type()
				--combine_numerator_a_denominator_dict()
		"""
		if self.vertex_set==None or self.edge_set==None or self.d_matrix == None or self.vertex2no == None:
			sys.stderr.write("In GradientScorePrediction, fill in vertex_set, edge_set, d_matrix first by setup_cluster_info().\n")
			sys.exit(2)
		result_heap = []
		d_row = self.d_matrix[self.vertex2no[gene_no]]
		layer2no_of_vertices = self.cal_layer2no_of_vertices(d_row)
		layer2gradient_score_denominator, layer2edge_gradient_denominator = self.cal_denominator_dict[self.eg_d_type](\
			d_row, self.exponent, self.max_layer, self.norm_exp, layer2no_of_vertices)
		
		#10-23-05 comment out 
		#degree = sum([int(layer==1) for layer in d_row])	#to judge whether consider probability or not
		#check each go_no
		candidate_go_no_ls = self.get_candidate_go_nos(self.gene_no2go, self.go_no2depth, self.vertex2no, \
			layer2no_of_vertices, d_row, self.depth, self.min_layer1_associated_genes, self.min_layer1_ratio, \
			self.min_layer2_associated_genes, self.min_layer2_ratio)
		for go_no in candidate_go_no_ls:
			layer2gradient_score, layer2edge_gradient = self.cal_func_dict[self.eg_d_type](gene_no, go_no, self.edge_set, \
				self.vertex2no, d_row, self.gene_no2go, self.go_no2gene_no_set, self.exponent, self.score_list, self.max_layer, \
				self.no_of_total_edges, self.no_of_total_genes, self.no_of_unknown_genes)
			if self.eg_d_type==0:	#10-23-05 type 0 refers to edge gradient_score
				layer2gradient_score = combine_numerator_a_denominator_dict(layer2gradient_score, layer2gradient_score_denominator)
				gradient_score = sum(layer2gradient_score.values())
			else:	#use edge_gradient
				layer2gradient_score = combine_numerator_a_denominator_dict(layer2edge_gradient, layer2edge_gradient_denominator)
				gradient_score = sum(layer2gradient_score.values())
			heappush(result_heap, [-gradient_score, layer2edge_gradient, go_no])	#10-22-05, -gradient_score due to min heap
		if result_heap:	#not empty 10-22-05 revamped
			max_prediction = heappop(result_heap)
			#reverse the sign
			max_prediction[0] = -max_prediction[0]
			#replace the layer2edge_gradient with edge_gradient
			max_prediction[1] = sum(combine_numerator_a_denominator_dict(max_prediction[1], layer2edge_gradient_denominator).values())
			list_to_return = [max_prediction]
			for i in range(len(result_heap)):
				one_prediction = heappop(result_heap)
				if one_prediction[0] == -max_prediction[0]:
					one_prediction[0] = -one_prediction[0]
					one_prediction[1] = sum(combine_numerator_a_denominator_dict(one_prediction[1], layer2edge_gradient_denominator).values())
					list_to_return.append(one_prediction)
				else:	#other have lower gradient_score
					break
			return list_to_return
		else:
			return []
	


class MpiStatCluster:
	"""
	10-20-05
		largely copied from MpiPredictionFilter.py
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, input_fname=None,\
		jnput_fname=None, depth=5, min_layer1_associated_genes=2, min_layer1_ratio=0.5, min_layer2_associated_genes=2,\
		min_layer2_ratio=0.1, exponent=2, score_list='3,-1,0', \
		max_layer=5, norm_exp=0.8, eg_d_type=2, recurrence_x=5,recurrence_x_type=2, message_size=10000000, \
		go_no2edge_counter_list_fname=None, new_table=0, commit=0, debug=0, report=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_fname = input_fname
		self.jnput_fname = jnput_fname
		self.depth = int(depth)
		self.min_layer1_associated_genes = int(min_layer1_associated_genes)
		self.min_layer1_ratio = float(min_layer1_ratio)
		self.min_layer2_associated_genes = int(min_layer2_associated_genes)
		self.min_layer2_ratio = float(min_layer2_ratio)
		self.exponent = float(exponent)
		score_list = score_list.split(',')
		self.score_list = map(float, score_list)
		self.max_layer = int(max_layer)
		self.norm_exp = float(norm_exp)
		self.eg_d_type = int(eg_d_type)
		self.recurrence_x = float(recurrence_x)
		self.recurrence_x_type = int(recurrence_x_type)
		self.message_size = int(message_size)
		self.go_no2edge_counter_list_fname = go_no2edge_counter_list_fname
		self.new_table = int(new_table)
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
		
		self.edge_type2index = {(0,0):0, (0,1):1, (1,0):1, (0,2):2, (2,0):2, (1,1):3, (1,2):4, (2,1):4, (2,2):5 }
			#this is used to map to the index of the value of go_no2edge_counter_list,
			#1st digit is gene1, 2nd digit is gene2,(order doesn't matter), 0 means gene having this function,
			#1 means gene having the other function, 2 means unknown gene
	
	def input_handler(self, parameter_list, message_size, report=0):
		"""
		01-02-06
			input is from MpiBFSCluster.py's output
		"""
		if report:
			sys.stderr.write("Fetching stuff...\n")
		reader = parameter_list[0]
		block = []
		string_length = 0
		for row in reader:
			self.counter += 1
			row.insert(0, self.counter)	#counter is used as id = line no
			block.append(row)
			string_length += len(repr(row))	#the length to control MPI message size
			if string_length>=message_size:
				break
		if report:
			sys.stderr.write("Fetching done.\n")
		return block
	
	def node_fire_handler(self, communicator, data, parameter_list):
		"""
		10-21-05
			called by common.computing_node()
		10-22-05
			GradientScorePrediction changed its return value of get_prediction()
		01-02-06
			the stuff from input_handler() changed due to the input of the program
			starts to come from plain file
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		GradientScorePrediction_instance, functor = parameter_list
		data = cPickle.loads(data)
		no_of_clusters = 0
		no_of_predictions = 0
		prediction_ls = []
		for pattern in data:
			#01-02-06
			id, vertex_set_string, edge_set_string, d_matrix_string = pattern
			edge_set = edge_set_string[2:-2].split('], [')
			for i in range(len(edge_set)):
				edge_set[i] = edge_set[i].split(',')
				edge_set[i] = map(int, edge_set[i])
			no_of_edges = len(edge_set)
			vertex_set = vertex_set_string[1:-1].split(',')
			vertex_set = map(int, vertex_set)
			
			#Before call its get_prediction(), setup the vertex_set, edge_set, d_matrix
			GradientScorePrediction_instance.setup_cluster_info(vertex_set, edge_set, d_matrix_string)
			for gene_no in vertex_set:
				prediction_list = GradientScorePrediction_instance.get_prediction(gene_no)	#10-22-05
				for gradient_score, edge_gradient, go_no in prediction_list:	#not None, None means no prediction.
					#ask the judge_node if its' correct
					communicator.send(cPickle.dumps([gene_no, go_no], -1), communicator.size-2, 3)	#tag is 3
					data, source, tag = communicator.receiveString(communicator.size-2, 3)
					is_correct, is_correct_L1, is_correct_lca, lca_list = cPickle.loads(data)
					#01-02-06 reshuffle the row
					row = [id, gene_no, go_no, GradientScorePrediction_instance.depth, gradient_score, edge_gradient, \
						is_correct, is_correct_L1, is_correct_lca, one_dim_list2string(lca_list)]
					prediction_ls.append(row)
					no_of_predictions += 1
			no_of_clusters += 1
		sys.stderr.write("Node no.%s done with %s clusters, %s predictions.\n"%(node_rank, no_of_clusters, no_of_predictions))
		return prediction_ls
	
	def cleanup_handler(self, communicator):
		#tell the size-2 node to stop
		communicator.send("-1", communicator.size-2, 3)	#tag is 3
		#tell the last node to stop
		communicator.send("-1", communicator.size-1, 1)	#tag is 1
	
	def judge_node(self, communicator, curs, gene_stat_instance, node_distance_class):
		"""
		10-20-05
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s ready to judge prediction...\n"%node_rank)
		data, source, tag = communicator.receiveString(None, 3)	#tag is 3 for judge_node
		no_of_resting_nodes = 0	#to keep track how many computing_nodes have rested
		while 1:
			if data=="-1":
				no_of_resting_nodes += 1
				if self.debug:
					sys.stderr.write("From No.%s: node %s(%s-th) rested.\n"%(node_rank, source, no_of_resting_nodes))
				if no_of_resting_nodes==communicator.size-3:	#WATCH: its' size-3
					break
					if self.report:
						sys.stderr.write("node %s done, exit.\n"%node_rank)
			else:
				gene_no, go_no = cPickle.loads(data)
				if gene_no in gene_stat_instance.known_genes_dict:
					k_functions_set = gene_stat_instance.known_genes_dict[gene_no]
					is_correct = gene_stat_instance.direct_match(go_no, k_functions_set)
					is_correct_L1 = -2
					#gene_stat_instance.L1_match(go_no, k_functions_set, node_distance_class, curs)
					#L1 is a waste of time. 10-20-05
					is_correct_lca = gene_stat_instance.common_ancestor_deep_enough(go_no, k_functions_set, node_distance_class, curs)
				else:
					#unknown gene
					is_correct = -1
					is_correct_L1 = -1
					is_correct_lca = -1
					#clear lca_list
					gene_stat_instance.lca_list = []
				is_correct_ls = [is_correct, is_correct_L1, is_correct_lca, gene_stat_instance.lca_list]
				is_correct_ls_pickle = cPickle.dumps(is_correct_ls, -1)
				communicator.send(is_correct_ls_pickle, source, 3)
			data, source, tag = communicator.receiveString(None, 3)
		sys.stderr.write("Node no.%s done.\n"%node_rank)
	
	def output_node_handler(self, communicator, parameter_list, data):
		"""
		10-21-05
			called by common.output_node()
		01-02-06
			output goes to plain file, not database
		"""
		writer = parameter_list[0]
		prediction_ls = cPickle.loads(data)
		for row in prediction_ls:
			writer.writerow(row)

	
	def run(self):
		"""
		09-05-05
		10-23-05
			create views from old schema
			result goes to the new schema's p_gene_table
		
			(input_node)
				--db_connect()
				--form_schema_tables()
				--form_schema_tables()
				--get_gene_no2go_no_set()
				--get_go_no2depth()
				(pass data to computing_node)
			(computing_node)
				(take data from other nodes, 0 and size-1)
			(judge_node)
				--gene_stat()
				--db_connect()
				--gene_p_map_redundancy()
			(output_node)
				--db_connect()
				--form_schema_tables()
				--form_schema_tables()
				--MpiPredictionFilter()
				--MpiPredictionFilter_instance.createGeneTable()
				--get_go_no2edge_counter_list()(if necessary)
				(pass go_no2edge_counter_list to computing_node)
			
			(input_node)
				--fetch_cluster_block()
			(computing_node)
				--get_no_of_unknown_genes()
				--node_fire_handler()
				--cleanup_handler()
			--judge_node()
				--gene_stat_instance.(match functions)
			--output_node()
				--output_node_handler()
					--MpiPredictionFilter_instance.submit_to_p_gene_table()
		"""
		communicator = MPI.world.duplicate()
		node_rank = communicator.rank
		if node_rank == 0:
			(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
			"""
			#01-02-06
			old_schema_instance = form_schema_tables(self.input_fname)
			new_schema_instance = form_schema_tables(self.jnput_fname)
			"""
			gene_no2go_no = get_gene_no2go_no_set(curs)
			gene_no2go_no_pickle = cPickle.dumps(gene_no2go_no, -1)	#-1 means use the highest protocol
			go_no2depth = get_go_no2depth(curs)
			go_no2depth_pickle = cPickle.dumps(go_no2depth, -1)
			go_no2gene_no_set = get_go_no2gene_no_set(curs)
			go_no2gene_no_set_pickle = cPickle.dumps(go_no2gene_no_set, -1)
			for node in range(1, communicator.size-2):	#send it to the computing_node
				communicator.send(gene_no2go_no_pickle, node, 0)
				communicator.send(go_no2depth_pickle, node, 0)
				communicator.send(go_no2gene_no_set_pickle, node, 0)
		elif node_rank<=communicator.size-3:	#WATCH: last 2 nodes are not here.
			data, source, tag = communicator.receiveString(0, 0)
			gene_no2go_no = cPickle.loads(data)	#take the data
			data, source, tag = communicator.receiveString(0, 0)
			go_no2depth = cPickle.loads(data)
			data, source, tag = communicator.receiveString(0, 0)
			go_no2gene_no_set = cPickle.loads(data)
			data, source, tag = communicator.receiveString(communicator.size-1, 0)	#from the last node
			go_no2edge_counter_list = cPickle.loads(data)
			#choose a functor for recurrence_array
			functor_dict = {0: None,
				1: lambda x: int(x>=self.recurrence_x),
				2: lambda x: math.pow(x, self.recurrence_x)}
			functor = functor_dict[self.recurrence_x_type]
		elif node_rank == communicator.size-2:	#judge node
			gene_stat_instance = gene_stat(depth_cut_off=self.depth)
			(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
			gene_stat_instance.dstruc_loadin(curs)
			from gene_p_map_redundancy import gene_p_map_redundancy
			node_distance_class = gene_p_map_redundancy()			
		elif node_rank==communicator.size-1:	#establish connection before pursuing
			(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
			"""
			#01-02-06, input and output are all directed to files
			old_schema_instance = form_schema_tables(self.input_fname)
			new_schema_instance = form_schema_tables(self.jnput_fname)
			MpiPredictionFilter_instance = MpiPredictionFilter()
			MpiPredictionFilter_instance.view_from_table(curs, old_schema_instance.splat_table, new_schema_instance.splat_table)
			MpiPredictionFilter_instance.view_from_table(curs, old_schema_instance.mcl_table, new_schema_instance.mcl_table)
			MpiPredictionFilter_instance.view_from_table(curs, old_schema_instance.pattern_table, new_schema_instance.pattern_table)
			if self.new_table:
				MpiPredictionFilter_instance.createGeneTable(curs, new_schema_instance.p_gene_table)
			"""
			if self.go_no2edge_counter_list_fname:
				go_no2edge_counter_list = cPickle.load(open(self.go_no2edge_counter_list_fname,'r'))
			else:
				if self.eg_d_type==2:
					go_no2edge_counter_list = None
				else:
					gene_no2go_no = get_gene_no2go_no_set(curs)
					go_no2edge_counter_list = get_go_no2edge_counter_list(curs, gene_no2go_no, self.edge_type2index)
			go_no2edge_counter_list_pickle = cPickle.dumps(go_no2edge_counter_list, -1)
			for node in range(1, communicator.size-2):	#send it to the computing_node
				communicator.send(go_no2edge_counter_list_pickle, node, 0)
		
		mpi_synchronize(communicator)
		
		free_computing_nodes = range(1,communicator.size-2)	#exclude the last node
		if node_rank == 0:
			"""
			curs.execute("DECLARE crs CURSOR FOR SELECT id, vertex_set, edge_set, no_of_edges,\
			connectivity, unknown_gene_ratio, recurrence_array, d_matrix from %s"%(old_schema_instance.pattern_table))
			"""
			self.counter = 0	#01-02-06 counter is used as id
			reader = csv.reader(open(self.input_fname, 'r'), delimiter='\t')
			parameter_list = [reader]
			input_node(communicator, parameter_list, free_computing_nodes, self.message_size, \
				self.report, input_handler=self.input_handler)
			del reader
		elif node_rank in free_computing_nodes:
			no_of_unknown_genes = get_no_of_unknown_genes(gene_no2go_no)
			GradientScorePrediction_instance = GradientScorePrediction(gene_no2go_no, go_no2gene_no_set, go_no2depth, \
				go_no2edge_counter_list, no_of_unknown_genes, self.depth, self.min_layer1_associated_genes, \
				self.min_layer1_ratio, self.min_layer2_associated_genes, self.min_layer2_ratio, self.exponent, \
				self.score_list, self.max_layer, self.norm_exp, self.eg_d_type, self.debug)
			parameter_list = [GradientScorePrediction_instance, functor]
			computing_node(communicator, parameter_list, self.node_fire_handler, self.cleanup_handler, self.report)
		elif node_rank == communicator.size-2:
			self.judge_node(communicator, curs, gene_stat_instance, node_distance_class)
		elif node_rank==communicator.size-1:
			#01-02-06 output goes to plain file, not database
			writer = csv.writer(open(self.jnput_fname, 'w'), delimiter='\t')
			parameter_list = [writer]
			output_node(communicator, free_computing_nodes, parameter_list, self.output_node_handler, self.report)
			del writer

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:j:u:g:o:p:e:s:l:q:t:x:w:v:f:ncbr", ["help", \
			"hostname=", "dbname=", "schema=", "de="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	input_fname = None
	jnput_fname = None
	depth = 5
	min_layer1_associated_genes = 2
	min_layer1_ratio = 0.5
	min_layer2_associated_genes = 2
	min_layer2_ratio = 0.1
	exponent = 2
	score_list = '3,-1,0'
	max_layer = 5
	norm_exp = 0.8
	eg_d_type = 2
	recurrence_x = 5
	recurrence_x_type = 2
	message_size = 10000000
	go_no2edge_counter_list_fname = None
	new_table = 0
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
		elif opt in ("-i",):
			input_fname = arg
		elif opt in ("-j",):
			jnput_fname = arg
		elif opt in ("--de",):
			depth = int(arg)
		elif opt in ("-u",):
			min_layer1_associated_genes = int(arg)
		elif opt in ("-g",):
			min_layer1_ratio = float(arg)
		elif opt in ("-o",):
			min_layer2_associated_genes = int(arg)
		elif opt in ("-p",):
			min_layer2_ratio = float(arg)
		elif opt in ("-e",):
			exponent = int(arg)
		elif opt in ("-s",):
			score_list = arg
		elif opt in ("-l",):
			max_layer = int(arg)
		elif opt in ("-q",):
			norm_exp = float(arg)
		elif opt in ("-t",):
			eg_d_type = int(arg)
		elif opt in ("-x",):
			recurrence_x = float(arg)
		elif opt in ("-w",):
			recurrence_x_type = int(arg)
		elif opt in ("-v",):
			message_size = int(arg)
		elif opt in ("-f",):
			go_no2edge_counter_list_fname = arg
		elif opt in ("-n",):
			new_table = 1
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-b",):
			debug = 1
		elif opt in ("-r",):
			report = 1
	if schema and input_fname and jnput_fname:
		instance = MpiStatCluster(hostname, dbname, schema, input_fname, \
			jnput_fname, depth, min_layer1_associated_genes, \
			min_layer1_ratio, min_layer2_associated_genes, min_layer2_ratio, \
			exponent, score_list, max_layer, norm_exp, eg_d_type, recurrence_x, recurrence_x_type, message_size, \
			go_no2edge_counter_list_fname, new_table, commit, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
