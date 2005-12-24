#!/usr/bin/env mpipython
"""
Usage: MpiClusterBsStat.py -k SCHEMA -g GOOD_CLUSTER_TABLE
	-l CLUSTER_BS_TABLE [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-g ...,	GOOD_CLUSTER_TABLE
	-l ...,	CLUSTER_BS_TABLE
	-s ...,	no of clusters per transmission, 2000(default)
	-a ...,	min ratio of #associated genes for one bs_no vs cluster size, 1/3(default)
	-t ...,	top_number of scores to be kept, 5(default)
	-e ...,	degree_cut_off, 0.3(default), for fuzzyDense
	-p ...,	p_value_cut_off, 0.001(default)
	-o ...,	output_file, fuzzyDense.out(default)
	-x ...,	tax_id, 4932(default, yeast), for get_gene_id2gene_symbol
	-n,	CLUSTER_BS_TABLE is new
	-c,	commit the database transaction
	-b,	debug version.
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	mpirun -np 20 -machinefile ~/hostfile /usr/bin/mpipython ~/script/annot/bin/MpiClusterBsStat.py
	-k mm_fim_97 -g good_clusters -l cluster_bs
	
Description:
	Program to do binding site enrichment analysis for a module.
	
"""

import sys, os, getopt, csv, math, Numeric, cPickle
sys.path += [os.path.expanduser('~/script/annot/bin')]
from Scientific import MPI
from codense.common import mpi_synchronize, db_connect, output_node, get_gene_id2gene_symbol, \
	dict_map, get_dataset_no2desc, input_node, computing_node
from sets import Set
from TF_functions import cluster_bs_analysis
from MpiCrackSplat import MpiCrackSplat	#12-18-05	for fill_edge2encodedOccurrence()
from fuzzyDense import fuzzyDense	#12-18-05	for fuzzyDense()

class MpiClusterBsStat:
	"""
	09-19-05
	"""
	def __init__(self,hostname='zhoudb', dbname='graphdb', schema=None, good_cluster_table=None,\
		cluster_bs_table=None,  size=2000, ratio_cutoff=1/3.0, top_number=5, degree_cut_off=0.3, p_value_cut_off=0.001, \
		output_file='fuzzyDense.out', tax_id=4932, new_table=0, commit=0, debug=0, report=0):
		"""
		12-15-05
			add p_value_cut_off
		12-20-05
			add output_file and tax_id, degree_cut_off
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.good_cluster_table = good_cluster_table
		self.cluster_bs_table = cluster_bs_table
		self.size = int(size)
		self.ratio_cutoff = float(ratio_cutoff)
		self.top_number = int(top_number)
		self.degree_cut_off = float(degree_cut_off)
		self.p_value_cut_off = float(p_value_cut_off)
		self.output_file = output_file
		self.tax_id = int(tax_id)
		self.new_table = int(new_table)
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
	
	def get_gene_no2bs_no_block(self, curs):
		"""
		09-19-05
			input_node takes care of this
			format:
				[gene_no1, bs_no11, bs_no12, ..., -1, gene_no2, gene_no21, ..., -1]
			-1 is the dilimiter
		"""
		sys.stderr.write("Getting gene_no2bs_no_block...\n")
		gene_no2bs_no_block = []
		curs.execute("DECLARE crs0 CURSOR FOR select gene_no, trans_fac from gene")
		curs.execute("fetch 5000 from crs0")
		rows = curs.fetchall()
		while rows:
			for row in rows:
				gene_no, trans_fac_ls = row
				if trans_fac_ls:	#throw away gene_nos with no trans_fac
					trans_fac_ls = trans_fac_ls[1:-1].split(',')
					trans_fac_ls = map(int, trans_fac_ls)
					gene_no2bs_no_block.append(gene_no)
					for bs_no in trans_fac_ls:
						gene_no2bs_no_block.append(bs_no)
					gene_no2bs_no_block.append(-1)	#delimiter
			curs.execute("fetch 5000 from crs0")
			rows = curs.fetchall()
		gene_no2bs_no_block = Numeric.array(gene_no2bs_no_block)
		sys.stderr.write("Done.\n")
		return gene_no2bs_no_block
	
	def construct_two_dicts(self, node_rank, gene_no2bs_no_block):
		"""
		09-19-05
			computing nodes need to have two dicts
		"""
		sys.stderr.write("Node no.%s constructing two dicts...\n"%(node_rank))
		gene_no_and_bs_no_list = []
		gene_no2bs_no_set = {}
		bs_no2gene_no_set = {}
		for no in gene_no2bs_no_block:
			if no==-1:
				gene_no = gene_no_and_bs_no_list[0]
				for bs_no in gene_no_and_bs_no_list[1:]:
					if gene_no not in gene_no2bs_no_set:
						gene_no2bs_no_set[gene_no] = Set()
					gene_no2bs_no_set[gene_no].add(bs_no)
					if bs_no not in bs_no2gene_no_set:
						bs_no2gene_no_set[bs_no] = Set()
					bs_no2gene_no_set[bs_no].add(gene_no)
				gene_no_and_bs_no_list = []
			else:
				gene_no_and_bs_no_list.append(no)
		sys.stderr.write("Node no.%s done with constructing two dicts.\n"%(node_rank))
		return gene_no2bs_no_set, bs_no2gene_no_set
	
	def fetch_cluster_block(self, curs, size):
		"""
		09-19-05
		12-18-05
			add recurrence_array
			
			format:
				[mcl_id1, gene_no, gene_no, ..., -1, recurrence, recurrence, ..., -2, mcl_id2, gene_no, gene_no, ...]
		"""
		curs.execute("fetch %s from crs"%(size))
		rows = curs.fetchall()
		cluster_block = []
		for row in rows:
			mcl_id, vertex_list, recurrence_array = row
			vertex_list = vertex_list[1:-1].split(',')
			vertex_list = map(int, vertex_list)
			recurrence_array = recurrence_array[1:-1].split(',')
			recurrence_array = map(int, recurrence_array)
			
			cluster_block.append(mcl_id)
			for vertex in vertex_list:
				cluster_block.append(vertex)
			cluster_block.append(-1)
			
			for recurrence in recurrence_array:
				cluster_block.append(recurrence)
			cluster_block.append(-2)
			
		cluster_block = Numeric.array(cluster_block)
		return cluster_block
	
	def input_node(self, communicator, curs, good_cluster_table, size):
		"""
		10-19-05
			ask output_node for a free_computing_node
		"""
		sys.stderr.write("Reading clusters from tables...\n")
		node_rank = communicator.rank
		curs.execute("DECLARE crs CURSOR FOR select distinct mcl_id, vertex_set, recurrence_array\
			from %s "%(good_cluster_table))
		cluster_block = self.fetch_cluster_block(curs, size)
		counter = 0	#10-19-05
		while cluster_block:
			communicator.send(Numeric.array([-2.0]), communicator.size-1, 1)	#10-19-05 WATCH: tag is 1, to the output_node. Message is array,Float.
			free_computing_node, source, tag = communicator.receiveString(communicator.size-1, 2)	#10-19-05
				#WATCH: tag is 2, from the output_node
			#node_to_receive_block = counter%(communicator.size-2)+1	#it's -2, work like this. 10-19-05 bad scheduling mechanism
			#	#"counter%(communicator.size-2)" is in range (0,size-3). Regard node 1 to size-2 as 0 to size-3. So need +1.
			communicator.send(cluster_block, int(free_computing_node), 0)	#10-19-05	#WATCH: int()
			if self.report:
				sys.stderr.write("block %s sent to %s.\n"%(counter, free_computing_node))	#10-19-05
			cluster_block = self.fetch_cluster_block(curs, size)
			counter += 1	#10-19-05
		#tell computing_node to exit the loop
		stop_signal = Numeric.zeros((1,1), Numeric.Int)
		stop_signal[0,0] = -1
		for node in range(1, communicator.size-1):	#send it to the computing_node
			communicator.send(stop_signal, node, 0)
		sys.stderr.write("Node no.%s reading done\n"%(node_rank))
	
	def encode_result_block(self, mcl_id, ls_to_return):
		block = [mcl_id]
		for ls in ls_to_return:
			score, score_type, bs_no_list, gene_no_list, global_ratio, local_ratio, expected_ratio, unknown_ratio = ls
			no_of_bs_nos = len(bs_no_list)
			no_of_genes = len(gene_no_list)
			row = [no_of_bs_nos, no_of_genes] + bs_no_list + gene_no_list + \
				[global_ratio, local_ratio, expected_ratio, unknown_ratio, score, score_type, -1]
			block += row
		return Numeric.array(block, Numeric.Float)
	
	def node_fire(self, communicator, data, gene_no2bs_no_set, bs_no2gene_no_set, ratio_cutoff, \
		top_number, p_value_cut_off, fuzzyDense_instance, degree_cut_off):
		"""
		09-19-05
		12-15-05
			add p_value_cut_off
		12-18-05
			fuzzyDense
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		tmp_list =[]
		for no in data:
			if no==-1:
				mcl_id = tmp_list[0]
				gene_no_list = tmp_list[1:]
				tmp_list = []
			elif no==-2:
				recurrence_array = tmp_list
				core_vertex_list, on_dataset_index_ls =  fuzzyDense_instance.get_core_vertex_set(gene_no_list, recurrence_array, degree_cut_off)
				if core_vertex_list:
					ls_to_return = cluster_bs_analysis(core_vertex_list, gene_no2bs_no_set, bs_no2gene_no_set, ratio_cutoff, \
						top_number, p_value_cut_off)
				else:
					ls_to_return = []
				result_block = self.encode_result_block(mcl_id, ls_to_return)
				communicator.send(result_block, communicator.size-1, 1)
				tmp_list = []
			else:
				tmp_list.append(no)
		sys.stderr.write("Node no.%s done.\n"%node_rank)
	
	def computing_node_handler(self, communicator, data, parameter_list):
		"""
		12-20-05
			called common's computing_node()
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		gene_no2bs_no_set, bs_no2gene_no_set, ratio_cutoff, \
		top_number, p_value_cut_off, fuzzyDense_instance, degree_cut_off = parameter_list
		result = []
		for row in data:
			id, vertex_set, recurrence_array = row
			vertex_set = vertex_set[1:-1].split(',')
			vertex_set = map(int, vertex_set)
			recurrence_array = recurrence_array[1:-1].split(',')
			recurrence_array = map(float, recurrence_array)
			core_vertex_list, on_dataset_index_ls =  fuzzyDense_instance.get_core_vertex_set(vertex_set, recurrence_array, degree_cut_off)
			if core_vertex_list:
				ls_to_return = cluster_bs_analysis(core_vertex_list, gene_no2bs_no_set, bs_no2gene_no_set, ratio_cutoff, \
					top_number, p_value_cut_off)
				if ls_to_return:
					result.append([id, core_vertex_list, on_dataset_index_ls, ls_to_return])
		sys.stderr.write("Node no.%s done.\n"%node_rank)
		return result
	
	def computing_node(self, communicator, gene_no2bs_no_set, bs_no2gene_no_set, ratio_cutoff, \
		top_number, p_value_cut_off, fuzzyDense_instance, degree_cut_off):
		"""
		09-19-05
		10-21-05
			stop_signal's dimension is 1. In fact, dimension (1,1) is same as (1).
		12-15-05
			add p_value_cut_off
		12-18-05
		"""
		node_rank = communicator.rank
		data, source, tag, count = communicator.receive(Numeric.Int, 0, 0)	#get data from node 0
		while 1:
			if data[0]==-1:
				if self.debug:
					sys.stderr.write("node %s breaked.\n"%node_rank)
				break
			else:
				self.node_fire(communicator, data, gene_no2bs_no_set, bs_no2gene_no_set, ratio_cutoff, \
					top_number, p_value_cut_off, fuzzyDense_instance, degree_cut_off)
			data, source, tag, count = communicator.receive(Numeric.Int, 0, 0)	#get data from node 0
		#tell the last node to stop
		stop_signal = Numeric.zeros((1), Numeric.Float)
		stop_signal[0] = -1.0
		communicator.send(stop_signal, communicator.size-1, 1)
	
	def create_cluster_bs_table(self, curs, cluster_bs_table):
		sys.stderr.write("Creating %s...\n"%cluster_bs_table)
		curs.execute("create table %s(\
			id	serial primary key,\
			mcl_id	integer,\
			bs_no_list	integer[],\
			gene_no_list	integer[],\
			global_ratio	float,\
			local_ratio	float,\
			expected_ratio	float,\
			unknown_ratio	float,\
			score	float,\
			score_type	integer)"%(cluster_bs_table))
		sys.stderr.write("%s created.\n"%cluster_bs_table)
	
	def submit_cluster_bs_table(self, communicator, parameter_list, data):
		"""
		09-20-05
			data is Float
		"""
		curs, cluster_bs_table = parameter_list
		mcl_id = int(data[0])
		row = []
		for no in data[1:]:
			if no==-1.0:
				no_of_bs_nos = int(row[0])
				no_of_genes = int(row[1])
				offset1 = no_of_bs_nos + 2
				offset2 = no_of_bs_nos + no_of_genes + 2
				bs_no_list = row[2:offset1]
				bs_no_list = map(int, bs_no_list)
				gene_no_list = row[offset1: offset2]
				gene_no_list = map(int, gene_no_list)
				global_ratio, local_ratio, expected_ratio, unknown_ratio, score, score_type = row[offset2:]
				score_type = int(score_type)
				curs.execute("insert into %s(mcl_id, bs_no_list, gene_no_list, global_ratio, local_ratio, \
					expected_ratio, unknown_ratio, score, score_type) values(%s, '{%s}', '{%s}', %s, %s,\
					%s, %s, %s, %s)"%(cluster_bs_table, mcl_id, repr(bs_no_list)[1:-1], repr(gene_no_list)[1:-1],\
					global_ratio, local_ratio, expected_ratio, unknown_ratio, score, score_type))
				row = []
			else:
				row.append(no)
		if len(data)==1:	#no cluster_bs_analysis result
			curs.execute("insert into %s(mcl_id) values(%s)"%(cluster_bs_table, mcl_id))
	
	def output_cluster_bs_data(self, communicator, parameter_list, data):
		"""
		12-20-05 for darwin output
		
			out:=[
			[id, [fuzzyDense gene list], { [ [TF gene list1], [TF target gene list1], p-value], [ [TF gene list2], [TF target gene list2], p-value], ...}, \
				{[dataset_no1, description], [dataset_no2, description], ... }  ],
			[...],
			...
			[]]:
			
		"""
		outf, gene_id2symbol, dataset_no2desc = parameter_list
		data = cPickle.loads(data)
		for row in data:
			id, core_vertex_ls, on_dataset_index_ls, ls_to_return = row
			#prepare the dataset_no_desc_ls
			dataset_no_desc_ls = []
			for dataset_index in on_dataset_index_ls:
				dataset_no = dataset_index +1
				dataset_no_desc_ls.append([dataset_no, dataset_no2desc[dataset_no]])
			#prepare the tfbs_row_darwin_ls
			tfbs_row_darwin_ls = []
			for tfbs_row in ls_to_return:
				score, score_type, bs_no_list, target_gene_no_list, global_ratio, local_ratio, expected_ratio, unknown_ratio = tfbs_row
				bs_no_symbol_list = dict_map(gene_id2symbol, bs_no_list)
				target_gene_no_symbol_list = dict_map(gene_id2symbol, target_gene_no_list)
				tfbs_row_darwin_ls.append([bs_no_symbol_list, target_gene_no_symbol_list, score])
			#translate the core_vertex_ls
			core_vertex_symbol_ls = dict_map(gene_id2symbol, core_vertex_ls)
			
			#output them all
			outf.write('[%s, %s,{%s},{%s}],\n'%(id, repr(core_vertex_symbol_ls), repr(tfbs_row_darwin_ls)[1:-1], repr(dataset_no_desc_ls)[1:-1]))
	
	def run(self):
		"""
		09-05-05
			
			--db_connect()
			--get_gene_no2bs_no_block()
			--construct_two_dicts()
			
			--input_node()
				--fetch_cluster_block()
			--computing_node()
				--node_fire()
					--cluster_bs_analysis()
			--create_cluster_bs_table()
			--output_node()
				--submit_cluster_bs_table()
		"""
		communicator = MPI.world.duplicate()
		node_rank = communicator.rank
		free_computing_nodes = range(1,communicator.size-1)
		
		if node_rank == 0:
			(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
			gene_no2bs_no_block = self.get_gene_no2bs_no_block(curs)
			for node in range(1, communicator.size-1):	#send it to the computing_node
				communicator.send(gene_no2bs_no_block, node, 0)
			#12-18-05 get edge2encodedOccurrence
			MpiCrackSplat_instance = MpiCrackSplat()
			edge2encodedOccurrence = {}
			min_sup = 3	#need to expose them
			max_sup = 30
			no_of_datasets = MpiCrackSplat_instance.fill_edge2encodedOccurrence(self.hostname, self.dbname, \
				self.schema, edge2encodedOccurrence, min_sup, max_sup)
			edge2encodedOccurrence_pickle = cPickle.dumps(edge2encodedOccurrence, -1)
			for node in free_computing_nodes:	#send it to the computing_node
				communicator.send(edge2encodedOccurrence_pickle, node, 0)
			
		elif node_rank>0 and node_rank<communicator.size-1:
			data, source, tag, count = communicator.receive(Numeric.Int, 0, 0)
			gene_no2bs_no_set, bs_no2gene_no_set = self.construct_two_dicts(node_rank, data)
			#12-18-05
			data, source, tag = communicator.receiveString(0, 0)
			edge2encodedOccurrence = cPickle.loads(data)
			
		elif node_rank==communicator.size-1:	#establish connection before pursuing
			(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
			
			#12-20-05 for darwin output
			gene_id2symbol = get_gene_id2gene_symbol(curs, self.tax_id)
			dataset_no2desc = get_dataset_no2desc(curs)
			
			
		mpi_synchronize(communicator)
		
		if node_rank == 0:
			#self.input_node(communicator, curs, self.good_cluster_table, self.size)
			
			#12-20-05
			curs.execute("DECLARE crs CURSOR FOR select distinct id, vertex_set, recurrence_array\
				from %s "%(self.good_cluster_table))
			input_node(communicator, curs, free_computing_nodes, self.size, self.report)
			curs.execute("close crs")
			
		elif node_rank<=communicator.size-2:	#exclude the last node
			#12-18-05
			#fuzzyDense_instance = fuzzyDense(edge2encodedOccurrence)
			#degree_cut_off = 0.3
			
			#self.computing_node(communicator, gene_no2bs_no_set, bs_no2gene_no_set, self.ratio_cutoff, \
				#self.top_number, self.p_value_cut_off, fuzzyDense_instance, degree_cut_off)				#12-15-05 add p_value_cut_off
			
			#12-20-05
			fuzzyDense_instance = fuzzyDense(edge2encodedOccurrence)
			parameter_list = [gene_no2bs_no_set, bs_no2gene_no_set, self.ratio_cutoff, \
				self.top_number, self.p_value_cut_off, fuzzyDense_instance, self.degree_cut_off]
			computing_node(communicator, parameter_list, self.computing_node_handler, report=self.report)
			
		elif node_rank==communicator.size-1:
			"""
			#12-20-05 comment out
			if self.new_table:
				self.create_cluster_bs_table(curs, self.cluster_bs_table)
			parameter_list = [curs, self.cluster_bs_table]
			output_node(communicator, free_computing_nodes, parameter_list, self.submit_cluster_bs_table, report=self.report, type=Numeric.Float)
			if self.commit:
				curs.execute("end")
			"""
			outf = open(self.output_file, 'w')
			outf.write("out:=[\n")
			parameter_list = [outf, gene_id2symbol, dataset_no2desc]
			output_node(communicator, free_computing_nodes, parameter_list, self.output_cluster_bs_data, report=self.report)
			
			outf.write('[]]:\n')
			outf.close()


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:g:l:s:a:t:e:p:o:x:ncbr", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	good_cluster_table = None
	cluster_bs_table = None
	size = 2000
	ratio_cutoff = 1.0/3
	top_number = 5
	degree_cut_off = 0.3
	p_value_cut_off = 0.001
	output_file = 'fuzzyDense.out'
	tax_id = 4932
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
		elif opt in ("-g",):
			good_cluster_table = arg
		elif opt in ("-l",):
			cluster_bs_table = arg
		elif opt in ("-s",):
			size = int(arg)
		elif opt in ("-a",):
			ratio_cutoff = float(arg)
		elif opt in ("-t",):
			top_number = int(arg)
		elif opt in ("-e",):
			degree_cut_off = float(arg)
		elif opt in ("-p",):
			p_value_cut_off = float(arg)
		elif opt in ("-n",):
			new_table = 1
		elif opt in ("-o",):
			output_file = arg
		elif opt in ("-x",):
			tax_id = int(arg)
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-b",):
			debug = 1
		elif opt in ("-r",):
			report = 1
	if schema and good_cluster_table:
		if cluster_bs_table or output_file:
			instance = MpiClusterBsStat(hostname, dbname, schema, good_cluster_table,\
				cluster_bs_table, size, ratio_cutoff, top_number, degree_cut_off, p_value_cut_off, output_file, \
				tax_id, new_table, commit, debug, report)
			instance.run()
	else:
		print __doc__
		sys.exit(2)
