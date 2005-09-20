#!/usr/bin/env mpipython
"""
Usage: MpiClusterBsStat.py -k SCHEMA -n GENE_P_TABLE -p P_GENE_TABLE
	-m MCL_TABLE -b CLUSTER_BS_TABLE [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-g ...,	GENE_P_TABLE
	-p ...,	P_GENE_TABLE
	-m ...,	MCL_TABLE
	-l ...,	CLUSTER_BS_TABLE
	-s ...,	no of clusters per transmission, 2000(default)
	-a ...,	min ratio of #associated genes for one bs_no vs cluster size, 1/3(default)
	-t ...,	top_number of scores to be kept, 5(default)
	-n,	CLUSTER_BS_TABLE is new
	-c,	commit the database transaction
	-b,	debug version.
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	mpirun -np 20 -machinefile ~/hostfile /usr/bin/mpipython ~/script/annot/bin/MpiClusterBsStat.py
	-k mm_fim_97 -g gene_p_mm_fim_97m4x200_e5_a60
	-p p_gene_mm_fim_97m4x200_e5 -m mcl_mm_fim_97m4x200 -l cluster_bs
	
Description:
	Program to do binding site enrichment analysis for a module.
	
"""

import sys, os, getopt, csv, math, Numeric
sys.path += [os.path.expanduser('~/script/annot/bin')]
from Scientific import MPI
from codense.common import mpi_synchronize, db_connect
from sets import Set
from TF_functions import cluster_bs_analysis

class MpiClusterBsStat:
	"""
	09-19-05
	"""
	def __init__(self,hostname='zhoudb', dbname='graphdb', schema=None, gene_p_table=None,\
		p_gene_table=None, mcl_table=None, cluster_bs_table=None,  size=2000, ratio_cutoff=1/3.0, \
		top_number=5, new_table=0, commit=0, debug=0, report=0):
		"""
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.gene_p_table = gene_p_table
		self.p_gene_table = p_gene_table
		self.mcl_table = mcl_table
		self.cluster_bs_table = cluster_bs_table
		self.size = int(size)
		self.ratio_cutoff = float(ratio_cutoff)
		self.top_number = int(top_number)
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
			format:
				[mcl_id1, gene_no, gene_no, ..., -1, mcl_id2, gene_no, gene_no, ...]
		"""
		curs.execute("fetch %s from crs"%(size))
		rows = curs.fetchall()
		cluster_block = []
		for row in rows:
			mcl_id, vertex_list = row
			vertex_list = vertex_list[1:-1].split(',')
			vertex_list = map(int, vertex_list)
			cluster_block.append(mcl_id)
			for vertex in vertex_list:
				cluster_block.append(vertex)
			cluster_block.append(-1)
		cluster_block = Numeric.array(cluster_block)
		return cluster_block
	
	def input_node(self, communicator, curs, gene_p_table, p_gene_table, mcl_table, size):
		sys.stderr.write("Reading clusters from tables...\n")
		node_rank = communicator.rank
		curs.execute("DECLARE crs CURSOR FOR select m.mcl_id, m.vertex_set from %s m limit 8000"%(mcl_table))
			#,\
			#%s p, %s m where n.p_gene_id=p.p_gene_id and m.mcl_id=p.mcl_id"%(gene_p_table,
			#p_gene_table, mcl_table))
		cluster_block = self.fetch_cluster_block(curs, size)
		which_node = 0
		while cluster_block:
			node_to_receive_block = which_node%(communicator.size-2)+1	#it's -2, work like this. 
				#"which_node%(communicator.size-2)" is in range (0,size-3). Regard node 1 to size-2 as 0 to size-3. So need +1.
			communicator.send(cluster_block, node_to_receive_block, 0)
			if self.debug:
				sys.stderr.write("block %s sent to %s.\n"%(which_node, node_to_receive_block))
			cluster_block = self.fetch_cluster_block(curs, size)
			which_node += 1
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
	
	def node_fire(self, communicator, data, gene_no2bs_no_set, bs_no2gene_no_set, ratio_cutoff, top_number):
		"""
		09-19-05
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		gene_no_list =[]
		for no in data:
			if no==-1:
				mcl_id = gene_no_list[0]
				gene_no_list = gene_no_list[1:]
				ls_to_return = cluster_bs_analysis(gene_no_list, gene_no2bs_no_set, bs_no2gene_no_set, ratio_cutoff, top_number)
				result_block = self.encode_result_block(mcl_id, ls_to_return)
				communicator.send(result_block, communicator.size-1, 1)
				gene_no_list = []
			else:
				gene_no_list.append(no)
		sys.stderr.write("Node no.%s done.\n"%node_rank)
	
	def computing_node(self, communicator, gene_no2bs_no_set, bs_no2gene_no_set, ratio_cutoff, top_number):
		"""
		09-19-05
			
		"""
		node_rank = communicator.rank
		data, source, tag, count = communicator.receive(Numeric.Int, 0, 0)	#get data from node 0
		while data:
			if data[0]==-1:
				if self.debug:
					sys.stderr.write("node %s breaked.\n"%node_rank)
				break
			else:
				self.node_fire(communicator, data, gene_no2bs_no_set, bs_no2gene_no_set, ratio_cutoff, top_number)
			data, source, tag, count = communicator.receive(Numeric.Int, 0, 0)	#get data from node 0
		#tell the last node to stop
		stop_signal = Numeric.zeros((1,1), Numeric.Float)
		stop_signal[0,0] = -1.0
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
	
	def submit_cluster_bs_table(self, curs, cluster_bs_table, data):
		"""
		09-20-05
			data is Float
		"""
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
				gene_no_list = row[offset1, offset2]
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
			
	
	def output_node(self, communicator, curs, cluster_bs_table):
		"""
		09-19-05
			
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s ready to accept output...\n"%node_rank)
		data, source, tag, count = communicator.receive(Numeric.Float, None, 1)
		no_of_resting_nodes = 0	#to keep track how many computing_nodes have rested
		sys.std .write("The data is %s.\n"%repr(data))
		while data:
			if data[0]==-1.0:
				no_of_resting_nodes += 1
				if self.debug:
					sys.stderr.write("node %s(%s-th) rested.\n"%(source, no_of_resting_nodes))
				if no_of_resting_nodes==communicator.size-2:	#all computing_nodes have stopped, i'm done.
					break
					if self.debug:
						sys.stderr.write("node %s output finished.\n"%node_rank)
			else:
				self.submit_cluster_bs_table(curs, cluster_bs_table, data)
			data, source, tag, count = communicator.receive(Numeric.Float, None, 1)
		sys.stderr.write("Node no.%s output done.\n"%node_rank)
	
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
		if node_rank == 0:
			(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
			gene_no2bs_no_block = self.get_gene_no2bs_no_block(curs)
			for node in range(1, communicator.size-1):	#send it to the computing_node
				communicator.send(gene_no2bs_no_block, node, 0)
		elif node_rank>0 and node_rank<communicator.size-1:
			data, source, tag, count = communicator.receive(Numeric.Int, 0, 0)
			gene_no2bs_no_set, bs_no2gene_no_set = self.construct_two_dicts(node_rank, data)
		elif node_rank==communicator.size-1:	#establish connection before pursuing
			(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		mpi_synchronize(communicator)
		
		if node_rank == 0:
			self.input_node(communicator, curs, self.gene_p_table, self.p_gene_table, self.mcl_table, self.size)
		elif node_rank<=communicator.size-2:	#exclude the last node
			self.computing_node(communicator, gene_no2bs_no_set, bs_no2gene_no_set, self.ratio_cutoff, self.top_number)
		elif node_rank==communicator.size-1:
			if self.new_table:
				self.create_cluster_bs_table(curs, self.cluster_bs_table)
			self.output_node(communicator, curs, self.cluster_bs_table)
			if self.commit:
				curs.execute("end")

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:g:p:m:l:s:a:t:ncbr", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	gene_p_table = None
	p_gene_table = None
	mcl_table = None
	cluster_bs_table = None
	size = 2000
	ratio_cutoff = 1.0/3
	top_number = 5
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
		elif opt in ("-g"):
			gene_p_table = arg
		elif opt in ("-p"):
			p_gene_table = arg
		elif opt in ("-m"):
			mcl_table = arg
		elif opt in ("-l"):
			cluster_bs_table = arg
		elif opt in ("-s"):
			size = int(arg)
		elif opt in ("-a"):
			ratio_cutoff = float(arg)
		elif opt in ("-t"):
			top_number = int(arg)
		elif opt in ("-n"):
			new_table = 1
		elif opt in ("-c"):
			commit = 1
		elif opt in ("-b"):
			debug = 1
		elif opt in ("-r"):
			report = 1
	if schema and gene_p_table and p_gene_table and mcl_table and cluster_bs_table:
		instance = MpiClusterBsStat(hostname, dbname, schema, gene_p_table, p_gene_table,\
			mcl_table, cluster_bs_table, size, ratio_cutoff, top_number, new_table, commit, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
