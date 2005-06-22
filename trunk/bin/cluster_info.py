#!/usr/bin/env python
"""
Usage: cluster_info.py -k SCHEMA -t SPLAT_TABLE -m MCL_TABLE [OPTION] OUTPUT_FILE

Option:
	R_FILE is the file to store the R code.
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	splat_result table
	-m ..., --mcl_table=...	mcl_result corresponding to table above
	-h, --help              show this help

Examples:
	cluster_info.py -k sc_54_6661 -t splat_result -m mcl_result /tmp/darwin.input
	
Description:
	This program outputs cluster information to be read by Jasmine's darwin.
	
	Also it's a backend for GuiAnalyzer.py
"""


import sys, os, psycopg, getopt, csv, math
from sets import Set
from rpy import r
from numarray import array
from numarray import reshape
from numarray import greater_equal, greater
from codense.common import db_connect
from codense.codense2db import codense2db
from codense.codense2db import cluster_dstructure
from codense.common import parse_splat_table_edge_set
from codense.common import get_no_of_total_genes
from graphlib import Graph
from CoexprFromCooccu  import CoexprFromCooccu

class cluster_info:
	"""
	04-06-05
		This class investigates the information about a cluster, such
		as recurrence vector, connectivity vector, edge_set's cor_vector
		or sig_vector.
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, table=None, mcl_table=None, \
			gene_p_table=None, gene_table=None, function=0, functioncolor='green', centralnode=1, mcl_id=1, \
			type=1, output_fname=None, plot_type="dot", label=1):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.table = table
		self.mcl_table = mcl_table
		self.gene_p_table = gene_p_table
		self.gene_table = gene_table
		self.function = int(function)
		self.functioncolor = functioncolor
		self.centralnode = int(centralnode)
		self.mcl_id = int(mcl_id)
		self.type = int(type)
		self.output_fname = output_fname
		self.plot_type = plot_type
		self.label = int(label)
	
		"""
		04-06-05
			other initializations
		"""
		#the table for edge_correlation_vector
		self.edge_table = 'edge_cor_vector'
		#mapping between go_no and go_id
		self.go_no2go_id = {}
		#mapping between go_no and go's name
		self.go_no2go_name = {}
		#mapping between gene_no and gene_id
		self.gene_no2gene_id = {}
		self.gene_id2gene_no = {}
		self.global_gene_to_go_dict = {}
		self.label_dict = {}

		self.order_1st_id2all_clusters = {}
		self.codense2db_instance = codense2db()
	
	def get_basic_cluster_dstructure(self, curs, mcl_id, splat_table, mcl_table):
		"""
		04-06-05
		
		"""
		sys.stderr.write("Getting the basic information of cluster no.%s..."%mcl_id)
		unit = cluster_dstructure()
		curs.execute("select m.mcl_id, m.vertex_set, m.connectivity, 0,\
			m.recurrence_array, s.edge_set, s.connectivity from %s m, %s s where m.splat_id=s.splat_id and \
			m.mcl_id=%s"\
			%(mcl_table, splat_table, mcl_id))	#06-20-05	connectivity_original faked to be 0
		rows = curs.fetchall()
		if len(rows)>0:
			for row in rows:
				unit.cluster_id = row[0]
				vertex_set = row[1][1:-1].split(',')
				unit.vertex_set = map(int, vertex_set)
				unit.connectivity = row[2]
				unit.connectivity_original = row[3]
				recurrence_array = row[4][1:-1].split(',')
				unit.recurrence_array = map(float, recurrence_array)
				unit.edge_set = parse_splat_table_edge_set(row[5])
				unit.splat_connectivity = row[6]
			sys.stderr.write("Done.\n")
		else:
			unit = None
			sys.stderr.write("Cluster: %s not found.\n"%mcl_id)
		return unit

	def get_go_functions_of_this_gene_set(self, curs, gene_set, gene_table='gene'):
		"""
		04-06-05
			input: gene_set
			output: go_no2association_genes
		"""
		sys.stderr.write("Getting go_functions of the gene set...")
		go_no2association_genes = {}
		for gene_no in gene_set:
			curs.execute("select gene_no,go_functions from %s where gene_no=%s"%(gene_table, gene_no))
			rows = curs.fetchall()
			for row in rows:
				go_functions_list = row[1][1:-1].split(',')
				go_functions_list = map(int, go_functions_list)
				for go_no in go_functions_list:
					if go_no not in go_no2association_genes:
						go_no2association_genes[go_no] = []
					go_no2association_genes[go_no].append(gene_no)
		sys.stderr.write("Done.\n")
		return go_no2association_genes
	
	def get_information_of_go_functions(self, curs, go_no2association_genes, cluster_size, \
		no_of_total_genes, p_value_cut_off=0, go_table='go'):
		"""
		04-06-05
			input: go_no_list
			output: go_no2information
			
			information includes go_no, go_id, name, depth, no_of_associated genes
		"""
		sys.stderr.write("Getting information about a list of go_nos...")
		go_no2information = {}
		for go_no,association_genes in go_no2association_genes.iteritems():
			no_of_associated_genes = len(association_genes)
			curs.execute("select go_no, go_id, name, depth, array_upper(gene_array,1) from %s \
				where go_no=%s"%(go_table, go_no))
			rows = curs.fetchall()
			for row in rows:
				p_value = r.phyper(no_of_associated_genes-1, row[-1],no_of_total_genes-row[-1], cluster_size,lower_tail = r.FALSE)
				if p_value_cut_off:	#non zero, needs cut some p-values
					if p_value>p_value_cut_off:
						continue
				go_no2information[go_no] = list(row) + [no_of_associated_genes, p_value]	#go_no, go_id, name, depth, population size, local size, p_value
					
		sys.stderr.write("Done.\n")
		return go_no2information
	
	def get_cor_sig_2d_list(self, curs, edge_set, edge_table='edge_cor_vector'):
		"""
		04-06-05
			input: edge_set
			output: (cor_2d_list, sig_2d_list)
			
			--one_dim2two_dim()
		"""
		sys.stderr.write("Getting cor_vector and sig_vector ...")
		no_of_edges = len(edge_set)
		(combined_cor_vector, combined_sig_vector) = self.codense2db_instance.get_combined_cor_vector(curs, edge_set, edge_table=edge_table)
		
		cor_2d_list = self.one_dim2two_dim(combined_cor_vector, no_of_edges)
		sig_2d_list = self.one_dim2two_dim(combined_sig_vector, no_of_edges)

		sys.stderr.write("Done.\n")
		return (cor_2d_list, sig_2d_list)
	
	def one_dim2two_dim(self, list_1d, x_dimension):
		list_2d = []
		y_dimension = len(list_1d)/x_dimension
		for i in range(x_dimension):
			list_2d.append(list_1d[i*y_dimension:(i+1)*y_dimension])
		return list_2d
	
	def graph_from_node_edge_set(self, node_set, edge_set, weight_list=None):
		"""
		04-06-05
			generate a graph from node_set and edge_set.
			edge_set and weight_list correspond to each other.
		"""
		sys.stderr.write("Generating graph from node_set and edge_set...")
		graph = Graph.Graph()
		for node in node_set:
			graph.add_node(node)
		for i in range(len(edge_set)):
			edge = edge_set[i]
			if weight_list:
				graph.add_edge(edge[0], edge[1], weight_list[i])
			else:
				graph.add_edge(edge[0], edge[1])
		sys.stderr.write("Done.\n")
		return graph
	
	def column_output(self, output_fname, list_2d):
		"""
		04-06-05
			a convinenient output function
		"""
		outf = open(output_fname, 'w')
		writer = csv.writer(outf, delimiter='\t')
		for list in list_2d:
			writer.writerow(list)
		del writer
		outf.close()

	def get_cluster_accuracy(self, curs, p_gene_table, mcl_id_list, p_value_cut_off=0.01):
		"""
		04-07-05
			get the accuracy, no_of_corrected predictions, known predictions, total predictions for each cluster
		"""
		accuracy2cluster = []
		for mcl_id in mcl_id_list:
			curs.execute("select is_correct_lca from %s where mcl_id=%s and p_value_cut_off<=%s"%(p_gene_table,\
				mcl_id, p_value_cut_off))
			rows = curs.fetchall()
			if rows:
				is_correct_lca_array = array(rows)
				correct_array = greater(is_correct_lca_array[:,0],0)
				known_array = greater_equal(is_correct_lca_array[:,0],0)
				accuracy = float(sum(correct_array))/float(sum(known_array))
				accuracy2cluster.append([accuracy, sum(correct_array), sum(known_array), len(correct_array), mcl_id]) 
		return accuracy2cluster

	
	def data_fetch(self, curs, splat_table, mcl_table, crs_no=0):
		"""
		04-17-05
			fetch cluster_dstructures for all clusters(Jasmine's request)	
		04-19-05
			1. return a mcl_id2cluster_dstructure
			2. crs_no
		"""
		mcl_id2cluster_dstructure = {}
		
		no_of_total_genes = get_no_of_total_genes(curs)
		sys.stderr.write("Getting the basic information for all clusters...\n")
		curs.execute("DECLARE crs%s CURSOR FOR select m.mcl_id, m.vertex_set, m.connectivity, 0,\
			m.recurrence_array, s.edge_set, s.connectivity, m.cooccurrent_cluster_id from %s m, %s s where \
			m.splat_id=s.splat_id"\
			%(crs_no, mcl_table, splat_table))	#06-20-05	connectivity_original faked to be 0
		curs.execute("fetch 5000 from crs%s"%crs_no)
		rows = curs.fetchall()
		while rows:
			for row in rows:
				unit = cluster_dstructure()
				unit.cluster_id = row[0]
				vertex_set = row[1][1:-1].split(',')
				unit.vertex_set = map(int, vertex_set)
				unit.connectivity = row[2]
				unit.connectivity_original = row[3]
				recurrence_array = row[4][1:-1].split(',')
				unit.recurrence_array = map(float, recurrence_array)
				unit.edge_set = parse_splat_table_edge_set(row[5])
				unit.splat_connectivity = row[6]
				unit.cooccurrent_cluster_id = row[7]
				unit.go_no2association_genes = self.get_go_functions_of_this_gene_set(curs, unit.vertex_set)
				unit.go_no2information = self.get_information_of_go_functions(curs, \
					unit.go_no2association_genes, len(unit.vertex_set), no_of_total_genes, p_value_cut_off=0.05)	#jasmine wants to cut some go-nos.
				unit.edge_cor_2d_list, unit.edge_sig_2d_list = self.get_cor_sig_2d_list(curs, unit.edge_set)
				
				mcl_id2cluster_dstructure[unit.cluster_id] = unit
				"""
				order_1st_id, order_2nd_id = map(int, unit.cooccurrent_cluster_id.split('.'))
				if order_1st_id not in self.order_1st_id2all_clusters:
					self.order_1st_id2all_clusters[order_1st_id] = {}
				if order_2nd_id not in self.order_1st_id2all_clusters[order_1st_id]:
					self.order_1st_id2all_clusters[order_1st_id][order_2nd_id] = []
				self.order_1st_id2all_clusters[order_1st_id][order_2nd_id].append(unit)
				"""
			curs.execute("fetch 5000 from crs%s"%crs_no)
			rows = curs.fetchall()
		sys.stderr.write("Done.\n")
		return mcl_id2cluster_dstructure
		
	def cluster_dstructure_output(self, curs, output_fname, order_1st_id2all_clusters):
		"""
		04-17-05
			output it in the format Jasmine's Darwin can read
		"""
		from codense.common import get_gene_no2gene_id
		gene_no2gene_id = get_gene_no2gene_id(curs)
		sys.stderr.write("Outputting cluster information...")
		outf = open(output_fname, 'w')
		str_tmp_list0 = []	#hold the 1st-order clusters
		for order_1st_id,all_2nd_order_clusters in order_1st_id2all_clusters.iteritems():
			str_tmp_list1 = []	#hold the 2nd-order clusters
			for order_2nd_id,cluster_list in all_2nd_order_clusters.iteritems():
				str_tmp_list2 = []	#hold the connected components
				for cluster in cluster_list:
					str_tmp = self.return_string_form_of_cluster_dstructure(cluster, gene_no2gene_id)
					str_tmp_list2.append(str_tmp)
				str_tmp_list1.append("[%s]"%','.join(str_tmp_list2))
			str_tmp_list0.append("[%s]"%",".join(str_tmp_list1))
		#'r:=' is for directly read in as an array
		outf.write("r:=[%s]:"%",".join(str_tmp_list0))
		outf.close()
		sys.stderr.write("Done.\n")
					
	def return_string_form_of_cluster_dstructure(self, cluster, gene_no2gene_id):
		"""
		04-19-05
			wrap the information jasmine wants out of a cluster_dstructure into a string form, which
			is in array-form of darwin.
		"""
		str_tmp ="["
		str_tmp+="[%s],\n"%cluster.cluster_id
		
		str_tmp_list = []
		for edge in cluster.edge_set:
			str_tmp_list.append("{'%s','%s'}"%(gene_no2gene_id[edge[0]], gene_no2gene_id[edge[1]]))
		str_tmp += "[%s],\n"%(','.join(str_tmp_list))
		
		str_tmp_list = []
		for go_no, information in cluster.go_no2information.iteritems():
			str_tmp_list.append("['%s','%s',%s,%s,%s,%s]\n"%(information[1], information[2], \
				information[3],information[4],information[5],information[6]))
		str_tmp += "[%s],"%(",".join(str_tmp_list))
		
		edge_sig_2d_array = array(cluster.edge_sig_2d_list)
		str_tmp_list = []
		for i in range(edge_sig_2d_array.shape[1]):
			str_tmp_list.append("%s\n"%repr(list(edge_sig_2d_array[:,i])))
		str_tmp += "[%s]]"%(",".join(str_tmp_list))
		
		return str_tmp
	
	def cluster_dstructure_output_with_both_hierarchy(self, curs, output_fname, \
		pre_2nd_cc_hierarchy, mcl_id2cluster_dstructure, mcl_id_2nd_order2cluster_dstructure):
		"""
		04-19-05
			jasmine wants to put 2nd-order clusters and its connected components into one file.
		"""
		from codense.common import get_gene_no2gene_id
		gene_no2gene_id = get_gene_no2gene_id(curs)
		sys.stderr.write("Outputting cluster information...")
		outf = open(output_fname, 'w')
		str_tmp_list0 = []	#hold the 1st-order clusters
		for pregraph_id,mcl_id_2nd_order_dict in pre_2nd_cc_hierarchy.iteritems():
			str_tmp_list1 = []	#hold the 2nd-order clusters
			for mcl_id_2nd_order,mcl_id_set in mcl_id_2nd_order_dict.iteritems():
				str_tmp_list2 = []	#hold the connected components
				#first one is the 2nd-order cluster
				str_tmp = self.return_string_form_of_cluster_dstructure(mcl_id_2nd_order2cluster_dstructure[mcl_id_2nd_order],\
					gene_no2gene_id)
				str_tmp_list2.append(str_tmp)
				for mcl_id in mcl_id_set:
					str_tmp = self.return_string_form_of_cluster_dstructure(mcl_id2cluster_dstructure[mcl_id],\
						gene_no2gene_id)
					str_tmp_list2.append(str_tmp)
				str_tmp_list1.append("[%s]"%','.join(str_tmp_list2))
			str_tmp_list0.append("[%s]"%",".join(str_tmp_list1))
		#'r:=' is for directly read in as an array
		outf.write("r:=[%s]:"%",".join(str_tmp_list0))
		outf.close()
		sys.stderr.write("Done.\n")
	
	def get_cluster_dstructure(self, curs, mcl_id, splat_table, mcl_table):
		"""
		04-18-05
			called by GuiAnalyzer.py
			
			--get_basic_cluster_dstructure()
			--get_go_functions_of_this_gene_set()
			--get_information_of_go_functions()
			--get_cor_sig_2d_list()
			--graph_from_node_edge_set()
			--column_output()
		"""
		no_of_total_genes = get_no_of_total_genes(curs)
		cluster  = self.get_basic_cluster_dstructure(curs, mcl_id, splat_table, mcl_table)
		if cluster:	#not None
			cluster.go_no2association_genes = self.get_go_functions_of_this_gene_set(curs, cluster.vertex_set)
			cluster.go_no2information = self.get_information_of_go_functions(curs, cluster.go_no2association_genes, \
				len(cluster.vertex_set), no_of_total_genes)
			cluster.edge_cor_2d_list, cluster.edge_sig_2d_list = self.get_cor_sig_2d_list(curs, cluster.edge_set)
			#graph = self.graph_from_node_edge_set(cluster.vertex_set, cluster.edge_set)
		return cluster
		
		"""
		print "vertex_set"
		print cluster.vertex_set
		print "edge_set"
		print cluster.edge_set
		recurrence_list_2d = ['recurrence_array']+cluster.recurrence_array
		recurrence_list_2d_1 =  ['recurrence_array_1']+cluster.recurrence_array
		recurrence_list_2d = [recurrence_list_2d, recurrence_list_2d_1]
		self.column_output('/tmp/yh/recurrence_array',recurrence_list_2d)

		print cluster.splat_connectivity
		print "connectivity"
		print cluster.connectivity
		print "connectivity_original"
		print cluster.connectivity_original
		cor_list_2d = []
		sig_list_2d = []
		for i in range(len(cluster.edge_set)):
			cor_list_2d.append([repr(cluster.edge_set[i])]+cluster.edge_cor_2d_list[i])
			sig_list_2d.append([repr(cluster.edge_set[i])]+cluster.edge_sig_2d_list[i])
		self.column_output('/tmp/yh/edge_cor_2d_list', cor_list_2d)
		self.column_output('/tmp/yh/edge_sig_2d_list', sig_list_2d)

		go_no_list_2d = []
		for go_no,information in cluster.go_no2information.iteritems():
			go_no_list_2d.append(list(information)+[len(cluster.go_no2association_genes[go_no])])
		#self.column_output('/tmp/yh/go_no_list_2d', go_no_list_2d)
		"""
	
	def run(self):
		"""
		04-18-05
			Serve for jasmine's darwin input.
		04-19-05
			changed to put 2nd-order clusters and its connected components into one file.
			
			--db_connect()
			--CoexprFromCooccu.data_fetch()
			--data_fetch()
			--data_fetch()
			--cluster_dstructure_output_with_both_hierarchy()

		"""
		(conn, curs) = db_connect(self.hostname, self.dbname, self.schema)
		
		e_splat_table = self.table+'e'
		e_mcl_table = self.mcl_table+'e'
		CoexprFromCooccu_instance = CoexprFromCooccu()
		pre_2nd_cc_hierarchy = CoexprFromCooccu_instance.data_fetch(curs, self.mcl_table, e_mcl_table)
		mcl_id2cluster_dstructure = self.data_fetch(curs, self.table,  self.mcl_table, crs_no=1)
		mcl_id_2nd_order2cluster_dstructure = self.data_fetch(curs, e_splat_table, e_mcl_table, crs_no=2)
		self.cluster_dstructure_output_with_both_hierarchy(curs, self.output_fname, pre_2nd_cc_hierarchy,\
			mcl_id2cluster_dstructure, mcl_id_2nd_order2cluster_dstructure)
		#self.cluster_dstructure_output(curs, self.output_fname, self.order_1st_id2all_clusters)
		

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema=", "table=", "mcl_table=", "gene_p_table=", \
		"gene_table=", "function=", "functioncolor=", "centralnode=", "mcl_id=", "type=", "plot_type=", "label="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:m:n:g:f:o:c:l:y:p:b:", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = None
	table = None
	mcl_table = None
	gene_p_table = None
	gene_table = None
	function = 0
	functioncolor = 'green'
	centralnode = 0
	mcl_id = 0
	type = 0
	plot_type = "dot"
	label = 1
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
		elif opt in ("-n", "--gene_p_table="):
			gene_p_table = arg
		elif opt in ("-g", "--gene_table"):
			gene_table = arg
		elif opt in ("-f", "--function"):
			function = int(arg)
		elif opt in ("-o", "--functioncolor"):
			functioncolor = arg
		elif opt in ("-c", "--centralnode"):
			centralnode = int(arg)
		elif opt in ("-l", "--mcl_id"):
			mcl_id = int(arg)
		elif opt in ("-y", "--type"):
			type = int(arg)
		elif opt in ("-p", "--plot_type"):
			plot_type = arg
		elif opt in ("-b", "--label"):
			label = int(arg)
			
	if schema and table and mcl_table:
		instance = cluster_info(hostname, dbname, schema, table, mcl_table, gene_p_table,\
			gene_table, function, functioncolor, centralnode, mcl_id, type, args[0], plot_type, label)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
