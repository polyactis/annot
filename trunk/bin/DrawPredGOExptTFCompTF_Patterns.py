#!/usr/bin/env python
"""
Usage: DrawPredGOExptTFCompTF_Patterns.py -k SCHEMA -i xxx -g xx -c xxx -e xxx -p xxx
	-o xxx [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-i ..., p_gene_table
	-g ...,	gene_p_table
	-n ...,	prot_interaction_table, 'mrinal_pi.intact_interaction'(default)
	-c ...,	comp_cluster_bs_table
	-e ...,	expt_cluster_bs_table
	-p ...,	pattern_table = 
	-o ...,	output_dir = 
	-m ...,	comp_tf_mapping_table = 'graph.gene_id2mt_no'
	-t ...,	expt_tf_mapping_table = 'graph.tf_mapping'
	-x ...,	tax_id, 9606(default)
	-b,	debug version.
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	python2.3 ./DrawPredGOExptTFCompTF_Patterns.py -k hs_fim_65 -c cluster_bs_hs_fim_65_m5x65s4l5e0p001geneid \
		-e cluster_bs_hs_fim_65_m5x65s4l5e0p001expt2 -p pattern_hs_fim_65_m5x65s4l5 -o /tmp/tf_go_patterns/ -i p_gene_hs_fim_65_n2s175_m5x65s4l5_ft2_e5
		-g gene_p_hs_fim_65_n2s175_m5x65s4l5_ft2_e5_000001a60
	
Description:
	Draw GO function graph, Experimental TF graph, Computational TF graph.
		(depends on whether they exist or not)
	color notation:
	 	green: standout but not associated
		yellow: standout and associated
		red: associated
		blue: not associated and not associated
	Edges colored red in GO function graph denote protein interaction.
		They might not be part of the pattern (compare it with TF graph).
	
	For the augmented PI graph:
		nodes and edges from the coexpr pattern are marked in 'green' and 'red' color, respectively
		nodes and edges added by protein interaction are marked in 'red' and 'gray' color, respectively
		overlapping edges are widened with 'red' color
"""


import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import sys, os, getopt, csv, math
from codense.common import db_connect, get_gene_id2gene_symbol, get_go_id2name, pg_1d_array2python_ls
from sets import Set
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import pylab
import Numeric

class DrawPredGOExptTFCompTF_Patterns:
	def __init__(self,hostname='zhoudb', dbname='graphdb', schema=None,  \
		input_fname=None, gene_p_table=None, prot_interaction_table=None, comp_cluster_bs_table=None, \
		expt_cluster_bs_table=None, gene_table='gene', go_table='go',\
		pattern_table=None, output_dir=None, \
		comp_tf_mapping_table='graph.gene_id2mt_no', expt_tf_mapping_table='graph.tf_mapping',\
		tax_id=9606, debug=0, report=0):
		"""
		2006-10-14
		2006-11-20
			add gene_p_table
			
			add prot_interaction_table
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_fname = input_fname
		self.gene_p_table = gene_p_table
		self.prot_interaction_table = prot_interaction_table
		self.comp_cluster_bs_table = comp_cluster_bs_table
		self.expt_cluster_bs_table = expt_cluster_bs_table
		self.gene_table = gene_table
		self.go_table = go_table
		self.pattern_table = pattern_table
		self.output_dir = output_dir
		self.comp_tf_mapping_table = comp_tf_mapping_table
		self.expt_tf_mapping_table = expt_tf_mapping_table
		self.tax_id = int(tax_id)
		self.debug = int(debug)
		self.report = int(report)
	

	def get_mt_no2gene_id_set(self, curs, tf_mapping_table, gene_id_set):
		sys.stderr.write("Getting mt_no2gene_id_set from %s..."%tf_mapping_table)
		curs.execute("select gene_id, mt_no from %s"%tf_mapping_table)
		rows = curs.fetchall()
		mt_no2gene_id_set = {}
		for row in rows:
			gene_id, mt_no = row
			if gene_id in gene_id_set:
				if mt_no not in mt_no2gene_id_set:
					mt_no2gene_id_set[mt_no] = Set()
				mt_no2gene_id_set[mt_no].add(int(gene_id))	#gene_id in these tf_mapping_table are char type
		sys.stderr.write("Done.\n")
		return mt_no2gene_id_set
	
	def get_mcl_id2mt_no_set(self, curs, cluster_bs_table):
		sys.stderr.write("Getting mcl_id2mt_no_set from %s ...\n"%cluster_bs_table)
		curs.execute("DECLARE crs CURSOR FOR select mcl_id, bs_no_list from %s"%cluster_bs_table )
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		mcl_id2mt_no_set = {}
		counter = 0
		while rows:
			for row in rows:
				mcl_id, bs_no_list = row
				bs_no_list = bs_no_list[1:-1].split(',')
				bs_no_list = map(int, bs_no_list)
				if mcl_id not in mcl_id2mt_no_set:
					mcl_id2mt_no_set[mcl_id] = Set()
				mcl_id2mt_no_set[mcl_id] |= Set(bs_no_list)
				counter += 1
			sys.stderr.write("%s%s"%('\x08'*10, counter))
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		curs.execute("close crs")
		sys.stderr.write("done\n")
		return mcl_id2mt_no_set
	
	def get_mcl_id2pred_go_id2gene_id_set(self, input_fname):
		#input_fname is the filtered one of haifeng's prediction output
		sys.stderr.write("Getting mcl_id2pred_go_id2gene_id_set from %s ..."%input_fname)
		reader = csv.reader(open(input_fname), delimiter='\t')
		mcl_id2pred_go_id2gene_id_set = {}
		for row in reader:
			vertex_id, pattern_id, go_id, is_known, is_correct, m1_score, p_value = row
			pattern_id = int(pattern_id)
			if pattern_id not in mcl_id2pred_go_id2gene_id_set:
				mcl_id2pred_go_id2gene_id_set[pattern_id] = {}
			if go_id not in mcl_id2pred_go_id2gene_id_set[pattern_id]:
				mcl_id2pred_go_id2gene_id_set[pattern_id][go_id] = Set()
			mcl_id2pred_go_id2gene_id_set[pattern_id][go_id].add(int(vertex_id))
		del reader
		sys.stderr.write("done\n")
		return mcl_id2pred_go_id2gene_id_set

	def get_mcl_id2pred_go_id2gene_id_set_from_db(self, curs, p_gene_table, gene_p_table):
		"""
		2006-11-20
			from p_gene_table
		"""
		sys.stderr.write("Getting mcl_id2pred_go_id2gene_id_set from %s ..."%p_gene_table)
		curs.execute("DECLARE crs1 CURSOR for select p.gene_no, p.mcl_id, go.go_id from %s p, %s g, go\
		where go.go_no=p.go_no and p.p_gene_id=g.p_gene_id"%(p_gene_table, gene_p_table))
		counter = 0
		curs.execute("fetch 1000 from crs1")
		rows = curs.fetchall()
		mcl_id2pred_go_id2gene_id_set = {}
		while rows:
			for row in rows:
				gene_no, pattern_id, go_id = row
				if pattern_id not in mcl_id2pred_go_id2gene_id_set:
					mcl_id2pred_go_id2gene_id_set[pattern_id] = {}
				if go_id not in mcl_id2pred_go_id2gene_id_set[pattern_id]:
					mcl_id2pred_go_id2gene_id_set[pattern_id][go_id] = Set()
				mcl_id2pred_go_id2gene_id_set[pattern_id][go_id].add(gene_no)
				counter += 1
			sys.stderr.write("%s%s"%('\x08'*30, counter))
			curs.execute("fetch 1000 from crs1")
			rows = curs.fetchall()
		curs.execute("close crs1")
		sys.stderr.write("done\n")
		return mcl_id2pred_go_id2gene_id_set
	
	def get_go_id2gene_set_from_db(self, curs, go_table):
		sys.stderr.write("Getting go_id2gene_set from %s..."%go_table)
		go_id2gene_set = {}
		curs.execute("select go_id, gene_array from %s"%go_table)
		rows = curs.fetchall()
		for row in rows:
			go_id, gene_array = row
			if gene_array:
				gene_array = gene_array[1:-1].split(',')
				gene_array = map(int, gene_array)
				go_id2gene_set[go_id] = Set(gene_array)
		sys.stderr.write("done\n")
		return go_id2gene_set
	
	def get_prot_interaction_graph(self, curs, prot_interaction_table, tax_id):
		"""
		2006-11-20
		"""
		sys.stderr.write("Getting protein interaction graph ...\n")
		curs.execute("DECLARE crs0 CURSOR for select gene_id_array, interaction_type_id\
			from %s where tax_id=%s"%(prot_interaction_table, tax_id))
		curs.execute("fetch 3000 from crs0")
		rows = curs.fetchall()
		prot_interaction_graph = nx.XGraph()
		counter = 0
		while rows:
			for row in rows:
				gene_id_array, interaction_type_id = row
				gene_id_array = pg_1d_array2python_ls(gene_id_array)
				if len(gene_id_array)>1:
					for i in range(len(gene_id_array)):
						for j in range(i+1, len(gene_id_array)):
							prot_interaction_graph.add_edge(gene_id_array[i], gene_id_array[j], interaction_type_id)
				counter += 1
			if self.report:
				sys.stderr.write("%s%s"%('\x08'*20, counter))
			curs.execute("fetch 3000 from crs0")
			rows = curs.fetchall()
		curs.execute("close crs0")
		sys.stderr.write("%s nodes and %s edges, Done.\n"%(nx.number_of_nodes(prot_interaction_graph),\
			nx.number_of_edges(prot_interaction_graph)))
		return prot_interaction_graph
	
	def draw_augmented_PI_graph(self, old_g, prot_interaction_graph, old_sub_label_map, gene_id2gene_symbol, output_fname_prefix):
		"""
		2006-12-16
			draw a graph based on old_g but augmented by prot_interaction_graph
				augment the old graph by node-pairwise searching against the interaction graph
			old nodes and edges are marked in 'green' and 'red' color, respectively
			overlapping edges are widened with 'red' color
		"""
		g = nx.XGraph()
		
		node_list = old_g.nodes()
		no_of_nodes = len(node_list)
		overlapping_edge_list = []
		non_standout_edge_list = []
		for m in range(no_of_nodes):
			for n in range(m+1, no_of_nodes):
				u = node_list[m]
				v = node_list[n]
				if prot_interaction_graph.has_node(u) and prot_interaction_graph.has_node(v):
					shortest_path_list = nx.shortest_path(prot_interaction_graph, u, v)
					if shortest_path_list:	#check the whole shortest path
						for i in range(len(shortest_path_list)-1):
							j = i+1
							if not g.has_edge(shortest_path_list[i], shortest_path_list[j]):
								g.add_edge(shortest_path_list[i], shortest_path_list[j], 1)
								if old_g.has_edge(shortest_path_list[i], shortest_path_list[j]):
									overlapping_edge_list.append((shortest_path_list[i], shortest_path_list[j]))
								else:
									non_standout_edge_list.append((shortest_path_list[i], shortest_path_list[j]))		
		standout_edge_list = []
		for (u,v) in old_g.edges():
			if not g.has_edge(u,v):
				standout_edge_list.append((u,v))
		
		sub_label_map = old_sub_label_map.copy()
		standout_node_list = []
		non_standout_node_list = []
		for v in g:
			if v not in old_g:
				sub_label_map[v] = gene_id2gene_symbol[v]
				non_standout_node_list.append(v)
		for v in old_g:
			standout_node_list.append(v)
			if v not in g:
				g.add_node(v)
		pos = nx.spring_layout(g)	#position is determined by the interaction graph
		
		if standout_edge_list:
			nx.draw_networkx_edges(g, pos, alpha=0.4, edge_color='r', edgelist=standout_edge_list)
		if non_standout_edge_list:
			nx.draw_networkx_edges(g, pos, alpha=0.4, edgelist=non_standout_edge_list)
		if overlapping_edge_list:
			nx.draw_networkx_edges(g, pos, alpha=0.4, edge_color='r', width=5, edgelist=overlapping_edge_list)
		if standout_node_list:
			nx.draw_networkx_nodes(g, pos, nodelist= standout_node_list, node_color='g', alpha=0.4)
		if non_standout_node_list:
			nx.draw_networkx_nodes(g, pos, nodelist= non_standout_node_list, alpha=0.4)
		nx.draw_networkx_labels(g, pos, labels=sub_label_map)
		pylab.savefig('%s.png'%(output_fname_prefix), dpi=300)
		pylab.clf()
		
	def draw_pattern(self, figure_no, old_g, pos, sub_label_map, title_map, go_id_or_mt_no_struct, go_id_or_mt_no2gene_id_set, \
		output_fname_prefix, is_go_function=0, prot_interaction_graph=None):
		"""
		2006-11-20
			add prot_interaction_graph
		"""
		g = old_g.copy()
		for key in go_id_or_mt_no_struct:
			figure_no += 1
			pylab.figure()
			pylab.title(title_map[key])
			standout_gene_id_list = []
			standout_and_associated_gene_id_list = []
			associated_gene_id_list = []
			other_gene_id_list = []
			for v in g:
				if is_go_function:
					if v in go_id_or_mt_no_struct[key] and v in go_id_or_mt_no2gene_id_set[key]:
						standout_and_associated_gene_id_list.append(v)
					elif v in go_id_or_mt_no_struct[key] and v not in go_id_or_mt_no2gene_id_set[key]:
						standout_gene_id_list.append(v)
					elif v in go_id_or_mt_no2gene_id_set[key]:
						associated_gene_id_list.append(v)
					else:
						other_gene_id_list.append(v)
				else:
					if v in go_id_or_mt_no2gene_id_set[key]:
						associated_gene_id_list.append(v)
					else:
						other_gene_id_list.append(v)
			if prot_interaction_graph:
				nodes_of_g = g.nodes()
				sub_prot_graph = prot_interaction_graph.subgraph(nodes_of_g)
				interaction_edge_list = []
				non_interaction_edge_list = []
				for (u, v) in g.edges():
					if not sub_prot_graph.has_edge(u,v):
						non_interaction_edge_list.append((u,v))
				for (u, v, interaction_type_id) in sub_prot_graph.edges():
					interaction_edge_list.append((u,v))
					#if not g.has_edge(u,v):	#expand g, not necesary
					#	print 'added'
					#	g.add_edge(u,v)
				nx.draw_networkx_edges(g, pos, alpha=0.4, edge_color='r', edgelist=interaction_edge_list)
				nx.draw_networkx_edges(g, pos, alpha=0.4, edgelist=non_interaction_edge_list)
			else:
				nx.draw_networkx_edges(g, pos, alpha=0.4)
			if standout_gene_id_list:
				nx.draw_networkx_nodes(g, pos, nodelist= standout_gene_id_list, node_color='g', alpha=0.4)
			if standout_and_associated_gene_id_list:
				nx.draw_networkx_nodes(g, pos, nodelist= standout_and_associated_gene_id_list, node_color='y', alpha=0.4)
			if associated_gene_id_list:
				nx.draw_networkx_nodes(g, pos, nodelist= associated_gene_id_list, node_color='r', alpha=0.4)
			if other_gene_id_list:
				nx.draw_networkx_nodes(g, pos, nodelist= other_gene_id_list, node_color='b', alpha=0.4)
			nx.draw_networkx_labels(g, pos, labels=sub_label_map)
			#nx.draw(g, pos, node_color=pylab.array(color_gene_id_list), labels=sub_label_map, alpha=0.4)
			pylab.savefig('%s_%s.png'%(output_fname_prefix, key), dpi=300)
			pylab.clf()
		return figure_no
		
	
	def draw_all_patterns(self, curs, mcl_id2pred_go_id2gene_id_set, comp_cluster_bs_table, expt_cluster_bs_table, gene_table, go_table,\
		pattern_table, output_dir, gene_id2gene_symbol, go_id2name,\
		comp_tf_mapping_table='graph.gene_id2mt_no', expt_tf_mapping_table='graph.tf_mapping',\
		prot_interaction_graph=None):
		"""
		2006-11-20
			make it one by one, from user input
			input_fname becomes mcl_id2pred_go_id2gene_id_set
		2006-11-21
			add prot_interaction_graph
		2006-12-16
			add --draw_augmented_PI_graph()
		"""
		sys.stderr.write("Getting gene_id set from %s..."%gene_table)
		curs.execute("select gene_id from %s"%gene_table)
		gene_id_set = Set()
		rows = curs.fetchall()
		for row in rows:
			gene_id_set.add(row[0])
		sys.stderr.write("Done.\n")
		
		comp_mt_no2gene_id_set = self.get_mt_no2gene_id_set(curs, comp_tf_mapping_table, gene_id_set)
		expt_mt_no2gene_id_set = self.get_mt_no2gene_id_set(curs, expt_tf_mapping_table, gene_id_set)
		go_id2gene_set = self.get_go_id2gene_set_from_db(curs, go_table)
		
		comp_mcl_id2mt_no_set = self.get_mcl_id2mt_no_set(curs, comp_cluster_bs_table)
		expt_mcl_id2mt_no_set = self.get_mcl_id2mt_no_set(curs, expt_cluster_bs_table)
		
		sys.stderr.write("Start to draw patterns....\n")
		pattern_id = raw_input("Please input a pattern id:")
		curs.execute("select id, edge_set from %s where id=%s"%(pattern_table, pattern_id))
		rows = curs.fetchall()
		counter = 0
		while rows:
			for row in rows:
				mcl_id, edge_set = row
				g = nx.Graph()
				edge_set = edge_set[2:-2].split('},{')
				for edge in edge_set:
					edge = edge.split(',')
					edge = map(int, edge)
					g.add_edge(edge[0], edge[1])
				pos = nx.spring_layout(g)
				#pos = nx.graphviz_layout(g)
				sub_label_map = {}
				for v in g:
					sub_label_map[v] = gene_id2gene_symbol[v]
				
				#2006-12-16
				output_fname_prefix = os.path.join(output_dir, 'id_%s_aug_prot_inter'%mcl_id)
				self.draw_augmented_PI_graph(g, prot_interaction_graph, sub_label_map, gene_id2gene_symbol, output_fname_prefix)
				
				if mcl_id in mcl_id2pred_go_id2gene_id_set:
					#1st, draw the GO association
					output_fname_prefix = os.path.join(output_dir, 'id_%s_go_func'%mcl_id)
					counter = self.draw_pattern(counter, g, pos, sub_label_map, go_id2name, mcl_id2pred_go_id2gene_id_set[mcl_id], go_id2gene_set, \
						output_fname_prefix, 1, prot_interaction_graph)
						# add prot_interaction_graph
				if mcl_id in comp_mcl_id2mt_no_set:
					#2nd, draw the comp mt_no association
					output_fname_prefix = os.path.join(output_dir, 'id_%s_comp_mt_no'%mcl_id)
					counter = self.draw_pattern(counter, g, pos, sub_label_map, gene_id2gene_symbol, comp_mcl_id2mt_no_set[mcl_id], comp_mt_no2gene_id_set, \
						output_fname_prefix, is_go_function=0)
				if mcl_id in expt_mcl_id2mt_no_set:
					#3rd, draw the expt mt_no association
					output_fname_prefix = os.path.join(output_dir, 	'id_%s_expt_mt_no'%mcl_id)
					counter = self.draw_pattern(counter, g, pos, sub_label_map, gene_id2gene_symbol, expt_mcl_id2mt_no_set[mcl_id], expt_mt_no2gene_id_set, \
						output_fname_prefix, is_go_function=0)
				if self.debug:
					sys.stderr.write("Exit on Debug\n")
					sys.exit(1)
			pattern_id = raw_input("Please input a pattern id:")
			curs.execute("select id, edge_set from %s where id=%s"%(pattern_table, pattern_id))
			rows = curs.fetchall()
		curs.execute("close crs")
		sys.stderr.write("done\n")
	
	def run(self):
		"""
		
		--db_connect()
		--get_gene_id2gene_symbol()
		--get_go_id2name()
		--get_mcl_id2pred_go_id2gene_id_set_from_db()
		--get_prot_interaction_graph()
		--draw_all_patterns()
			--get_mt_no2gene_id_set()
			--get_go_id2gene_set_from_db()
			--get_mcl_id2mt_no_set()
			--draw_augmented_PI_graph()
			--draw_pattern()
		
		2006-11-21
			add prot_interaction_graph
		"""
		conn, curs = db_connect(self.hostname, self.dbname, self.schema)
		gene_id2gene_symbol = get_gene_id2gene_symbol(curs, self.tax_id)
		go_id2name = get_go_id2name(curs)
		if not os.path.isdir(self.output_dir):
			os.makedirs(self.output_dir)
		#mcl_id2pred_go_id2gene_id_set = self.get_mcl_id2pred_go_id2gene_id_set(input_fname)
		mcl_id2pred_go_id2gene_id_set = self.get_mcl_id2pred_go_id2gene_id_set_from_db(curs, self.input_fname, self.gene_p_table)
		prot_interaction_graph = self.get_prot_interaction_graph(curs, self.prot_interaction_table, self.tax_id)
		self.draw_all_patterns(curs, mcl_id2pred_go_id2gene_id_set, self.comp_cluster_bs_table, self.expt_cluster_bs_table, \
			self.gene_table, self.go_table,\
			self.pattern_table, self.output_dir, gene_id2gene_symbol, go_id2name,\
			self.comp_tf_mapping_table, self.expt_tf_mapping_table, prot_interaction_graph)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:g:n:c:e:p:o:m:t:x:br", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	input_fname = None
	gene_p_table = None
	prot_interaction_table = 'mrinal_pi.intact_interaction'
	comp_cluster_bs_table = None
	expt_cluster_bs_table = None
	gene_table = 'gene'
	go_table = 'go'
	pattern_table = None
	output_dir = None
	comp_tf_mapping_table = 'graph.gene_id2mt_no'
	expt_tf_mapping_table = 'graph.tf_mapping'
	tax_id = 9606
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
		elif opt in ("-g",):
			gene_p_table = arg
		elif opt in ("-n",):
			prot_interaction_table = arg
		elif opt in ("-c",):
			comp_cluster_bs_table = arg
		elif opt in ("-e",):
			expt_cluster_bs_table = arg
		elif opt in ("-p",):
			pattern_table = arg
		elif opt in ("-o",):
			output_dir = arg
		elif opt in ("-m",):
			comp_tf_mapping_table = arg
		elif opt in ("-t",):
			expt_tf_mapping_table = arg
		elif opt in ("-x",):
			tax_id = int(arg)
		elif opt in ("-b",):
			debug = 1
		elif opt in ("-r",):
			report = 1
	if schema and input_fname and gene_p_table and comp_cluster_bs_table and expt_cluster_bs_table and pattern_table and output_dir:
			instance = DrawPredGOExptTFCompTF_Patterns(hostname, dbname, schema,  \
				input_fname, gene_p_table, prot_interaction_table, comp_cluster_bs_table, expt_cluster_bs_table, gene_table, go_table, \
				pattern_table, output_dir, \
				comp_tf_mapping_table, expt_tf_mapping_table, tax_id, debug, report)
			instance.run()
	else:
		print __doc__
		sys.exit(2)
