#!/usr/bin/env python
"""
functions or classes common to all programs

"""

import csv, os, sys, cPickle
sys.path += [os.path.join(os.path.expanduser('~/script/annot/bin'))]
import psycopg
from sets import Set

def index_plus_one(i):
	'''
	This small function is used in map() to increase the indices by 1.
	'''
	return int(i)+1
	
def divided_by_1000(i):
	return '%1.3f'%(float(i)/float(1000))


def get_haiyan_no2gene_no(table_file):
	sys.stderr.write("Getting haiyan_no2gene_no...")
	reader = csv.reader(file(table_file), delimiter='\t')
	haiyan_no2gene_no={}
	no = 0
	for row in reader:
		gene_no = int(row[1])
		haiyan_no2gene_no[no] = gene_no
		no += 1
	del reader
	sys.stderr.write("Done\n")
	return haiyan_no2gene_no


def get_haiyan_no2gene_id(table_file):
	sys.stderr.write("Getting haiyan_no2gene_id...")
	reader = csv.reader(file(table_file), delimiter='\t')
	haiyan_no2gene_id={}
	no = 0
	for row in reader:
		gene_id = row[0]
		haiyan_no2gene_id[no] = gene_id
		no += 1
	del reader
	sys.stderr.write("Done\n")
	return haiyan_no2gene_id


def dict_map(dict, ls, type=1):
	"""
	10-13-05
		add type 2 to return item itself if mapping is not available
	"""
	new_list = []
	for item in ls:
		value = dict.get(item)
		if value:
			new_list.append(value)
		elif type==2:
			new_list.append(item)
	return new_list

def db_connect(hostname, dbname, schema=None):
	"""
	02-28-05
		establish database connection, return (conn, curs).
		copied from CrackSplat.py
	03-08-05
		parameter schema is optional
	"""
	conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
	curs = conn.cursor()
	if schema:
		curs.execute("set search_path to %s"%schema)
	return (conn, curs)

def get_gene_no2direct_go(curs, schema=None, table='association'):
	"""
	03-02-05
		return gene_no2direct_go
	"""
	sys.stderr.write("Getting gene_no2direct_go...")
	gene_no2direct_go = {}
	if schema:
		curs.execute("set search_path to %s"%schema)
	curs.execute("select ge.gene_no, a.go_id from graph.association a, gene ge\
			where ge.gene_id=a.gene_id")
	rows = curs.fetchall()
	for row in rows:
		if row[0] in gene_no2direct_go:
			gene_no2direct_go[row[0]].append(row[1])
		else:
			gene_no2direct_go[row[0]] = [row[1]]
	sys.stderr.write("Done\n")
	return gene_no2direct_go
	

def get_gene_no2gene_id(curs, schema=None, table='gene'):
	"""
	03-02-05
		return gene_no2gene_id
	"""
	sys.stderr.write("Getting gene_no2gene_id...")
	gene_no2gene_id = {}
	if schema:
		curs.execute("set search_path to %s"%schema)
	curs.execute("select gene_no, gene_id from %s"%table)
	rows = curs.fetchall()
	for row in rows:
		gene_no2gene_id[row[0]] = row[1]
	sys.stderr.write("Done\n")
	return gene_no2gene_id


def get_gene_id2gene_no(curs, schema=None, table='gene'):
	"""
	03-02-05
		return gene_id2gene_no
	"""
	sys.stderr.write("Getting gene_id2gene_no...")
	gene_id2gene_no = {}
	if schema:
		curs.execute("set search_path to %s"%schema)
	curs.execute("select gene_no, gene_id from %s"%table)
	rows = curs.fetchall()
	for row in rows:
		gene_id2gene_no[row[1]] = row[0]
	sys.stderr.write("Done\n")
	return gene_id2gene_no


def get_go_no2go_id(curs, schema=None, table='go'):
	"""
	03-02-05
		return go_no2go_id
	"""
	sys.stderr.write("Getting go_no2go_id...")
	go_no2go_id = {}
	if schema:
		curs.execute("set search_path to %s"%schema)
	curs.execute("select go_no, go_id, depth from %s"%table)
	rows = curs.fetchall()
	for row in rows:
		go_no2go_id[row[0]] = row[1]
	sys.stderr.write("Done\n")
	return go_no2go_id

def get_go_no2name(curs, schema=None, table='go'):
	"""
	03-03-05
		return go_no2name
	"""
	sys.stderr.write("Getting go_no2name...")
	go_no2name = {}
	if schema:
		curs.execute("set search_path to %s"%schema)
	curs.execute("select go_no, name from %s"%table)
	rows = curs.fetchall()
	for row in rows:
		go_no2name[row[0]] = row[1]
	sys.stderr.write("Done\n")
	return go_no2name

def get_go_no2term_id(curs, schema=None, term_table='go.term'):
	"""
	03-04-05
		get the go_no2term_id dictionary
	03-22-05
		there're GO synonyms in table term_synonym,
		which should be linked to the source term_id of table term.
	"""
	sys.stderr.write("Getting go_no2term_id...")
	if schema:
		curs.execute("set search_path to %s"%schema)
	go_no2term_id = {}
	curs.execute("select g.go_no, t.name, t.id from go g, %s t where g.go_id=t.acc"%(term_table))
	rows = curs.fetchall()
	for row in rows:
		go_no2term_id[row[0]] = row[2]
	#acc_synonym => term_id
	curs.execute("select g.go_no, t.term_id from go.term_synonym t, go g where g.go_id=t.acc_synonym")
	rows = curs.fetchall()
	for row in rows:
		go_no2term_id[row[0]] = row[1]
	sys.stderr.write("Done\n")
	return go_no2term_id

def get_gene_no2go_no(curs, schema=None, gene_table='gene'):
	"""
	03-09-05
		get the gene_no2go_no
	"""
	sys.stderr.write("Getting gene_no2go_no...")
	if schema:
		curs.execute("set search_path to %s"%schema)
	
	gene_no2go_no = {}
	
	curs.execute("select gene_no,go_functions from %s"%gene_table)
	rows = curs.fetchall()
	for row in rows:
		gene_no2go_no[row[0]] = []
		go_functions_list = row[1][1:-1].split(',')
		for go_no in go_functions_list:
			gene_no2go_no[row[0]].append(int(go_no))
	sys.stderr.write("Done\n")
	return gene_no2go_no
	
def get_gene_no2go_no_set(curs, schema=None, gene_table='gene'):
	"""
	03-09-05
		get the gene_no2go_no
	"""
	sys.stderr.write("Getting gene_no2go_no...")
	if schema:
		curs.execute("set search_path to %s"%schema)
	
	gene_no2go_no = {}
	
	curs.execute("select gene_no,go_functions from %s"%gene_table)
	rows = curs.fetchall()
	for row in rows:
		gene_no2go_no[row[0]] = Set()
		go_functions_list = row[1][1:-1].split(',')
		for go_no in go_functions_list:
			gene_no2go_no[row[0]].add(int(go_no))
	sys.stderr.write("Done\n")
	return gene_no2go_no

def get_known_genes_dict(curs, schema=None, gene_table='gene'):
	"""
	03-09-05
		get the known_genes_dict
	"""
	sys.stderr.write("Getting known_genes_dict...")
	if schema:
		curs.execute("set search_path to %s"%schema)

	known_genes_dict = {}
	curs.execute("select gene_no,go_functions from %s where known=TRUE"%gene_table)
	rows = curs.fetchall()
	for row in rows:
		go_functions_list = row[1][1:-1].split(',')
		known_genes_dict[row[0]] = Set()
		for go_no in go_functions_list:
			known_genes_dict[row[0]].add(int(go_no))
	sys.stderr.write("Done\n")
	return known_genes_dict

def get_go_no2depth(curs, schema=None, table='go'):
	"""
	03-14-05
		get the go_no2depth
	"""
	sys.stderr.write("Getting go_no2depth...")
	if schema:
		curs.execute("set search_path to %s"%schema)
	go_no2depth = {}
	curs.execute("select go_no, depth from %s"%table)
	rows = curs.fetchall()
	for row in rows:
		go_no2depth[row[0]] = row[1]
	sys.stderr.write("Done\n")
	return go_no2depth

def get_go_term_id2go_no(curs, schema=None, go_table='go', term_table='go.term'):
	"""
	03-14-05
		get go_term_id2go_no
	"""
	sys.stderr.write("Getting go_term_id2go_no...")
	if schema:
		curs.execute("set search_path to %s"%schema)
	go_term_id2go_no = {}
	curs.execute("select g.go_no, t.name, t.id from %s g, %s t where g.go_id=t.acc"%(go_table, term_table))
	rows = curs.fetchall()
	for row in rows:
		go_term_id2go_no[row[2]] = row[0]
	
	sys.stderr.write("Done\n")
	return go_term_id2go_no
	
def get_go_term_id2depth(curs, schema=None, term_table='go.term'):
	"""
	03-14-05
		get go_term_id2depth
	"""
	sys.stderr.write("Getting go_term_id2depth...")
	if schema:
		curs.execute("set search_path to %s"%schema)
	go_term_id2depth = {}
	curs.execute("select id, depth from %s where depth NOTNULL"%(term_table))
	rows = curs.fetchall()
	for row in rows:
		go_term_id2depth[row[0]] = row[1]
	
	sys.stderr.write("Done\n")
	return go_term_id2depth


def get_go_id2go_no(curs, schema=None, table='go'):
	"""
	03-14-05
		return go_id2go_no
	"""
	sys.stderr.write("Getting go_id2go_no...")
	go_id2go_no = {}
	if schema:
		curs.execute("set search_path to %s"%schema)
	curs.execute("select go_no, go_id, depth from %s"%table)
	rows = curs.fetchall()
	for row in rows:
		go_id2go_no[row[1]] = row[0]
	sys.stderr.write("Done\n")
	return go_id2go_no

"""
07-06-05
"""
def get_gene_no2id(curs, table='gene_id_to_no', organism='Mus musculus'):
	"""
	07-06-05
		return gene_no2id
	"""
	sys.stderr.write("Getting gene_no2id...")
	gene_no2id = {}
	curs.execute("select gene_no, gene_id from graph.%s where organism='%s'"%(table, organism))
	rows = curs.fetchall()
	for row in rows:
		gene_no2id[row[0]] = row[1]
	sys.stderr.write("Done\n")
	return gene_no2id
	

def get_go_id2name(curs, table='term'):
	"""
	07-06-05
		return go_id2name
	"""
	sys.stderr.write("Getting go_id2name...")
	go_id2name = {}
	curs.execute("select acc, name from go.%s where term_type='biological_process'"%table)
	rows = curs.fetchall()
	for row in rows:
		go_id2name[row[0]] = row[1]
	sys.stderr.write("Done\n")
	return go_id2name


def get_go_id2term_id(curs, table='term'):
	"""
	07-06-05
		return go_id2name
	"""
	sys.stderr.write("Getting go_id2name...")
	go_id2term_id = {}
	curs.execute("select acc, id from go.%s where term_type='biological_process'"%table)
	rows = curs.fetchall()
	for row in rows:
		go_id2term_id[row[0]] = row[1]
	sys.stderr.write("Done\n")
	return go_id2term_id

def get_prediction_pair2lca_list(curs, schema=None, p_gene_table=None):
	"""
	03-15-05
		get prediction_pair2lca_list
	"""
	sys.stderr.write("Getting prediction_pair2lca_list...")
	prediction_pair2lca_list = {}
	if schema:
		curs.execute("set search_path to %s"%schema)
	curs.execute("select distinct gene_no, go_no, lca_list from %s where lca_list notnull"%p_gene_table)
	rows = curs.fetchall()
	for row in rows:
		prediction_pair = (row[0], row[1])
		lca_list = row[2][1:-1].split(',')
		lca_list = map(int, lca_list)
		if prediction_pair not in prediction_pair2lca_list:
			prediction_pair2lca_list[prediction_pair] = lca_list
		else:
			sys.stderr.write("Error: same prediction_pair: %s appears twice with different lca_list.\n"%repr(prediction_pair))
			sys.exit(3)
	sys.stderr.write("Done\n")
	return prediction_pair2lca_list

def get_gspan_graph(gspan_file, min_weight=None):
	"""
	04-03-05
		read in a graph from a file in gspan format
	"""
	sys.stderr.write("Getting graph from a gspan file %s..."%gspan_file)
	from graphlib import Graph
	graph = Graph.Graph()
	reader = csv.reader(open(gspan_file,'r'), delimiter=' ')
	for row in reader:
		if row[0] == 'e':
			gene_no1 = int (row[1])
			gene_no2 = int(row[2])
			weight = float(row[3])
			if min_weight:
				if weight>=min_weight:
					graph.add_edge(gene_no1, gene_no2, weight)
			else:
				graph.add_edge(gene_no1, gene_no2, weight)
	del reader
	sys.stderr.write("Done\n")
	return graph

def parse_splat_table_edge_set(edge_string):
	"""
	04-05-05
		parse the edge_string fetched from splat table and return a list
	"""
	edge_set = []
	edge_list = edge_string[2:-2].split('},{')
	for edge in edge_list:
		edge = edge.split(',')
		edge = map(int, edge)
		edge_set.append(edge)
	return edge_set

def foreach_cb(model, path, iter, pathlist):
	"""
	04-17-05
		used in gui listview, pathfinding.
	"""
	pathlist.append(path)	
	
def create_columns(treeview, label_list):
	"""
	04-17-05
		create columns in the treeview in the first refresh
	04-21-05
		remove the old columns and reset the model of treeview
	"""
	import gtk
	tvcolumn_dict = {}
	cell_dict = {}
	#remove old columns
	old_column_list = treeview.get_columns()
	for column in old_column_list:
		treeview.remove_column(column)
	treeview.set_model()
		
	for i in range(len(label_list)):
		tvcolumn_dict[i] = gtk.TreeViewColumn(label_list[i])	# create the TreeViewColumn to display the data
		treeview.append_column(tvcolumn_dict[i])	# add tvcolumn to treeview
		cell_dict[i] = gtk.CellRendererText()	# create a CellRendererText to render the data
		tvcolumn_dict[i].pack_start(cell_dict[i], True)	# add the cell to the tvcolumn and allow it to expand
		# set the cell "text" attribute to column 0 - retrieve text
		# from that column in liststore
		tvcolumn_dict[i].add_attribute(cell_dict[i], 'text', i)
		tvcolumn_dict[i].set_sort_column_id(i)	# Allow sorting on the column

def fill_treeview(treeview, liststore, list_2d, reorderable=True):
	"""
	04-17-05
	05-20-05
		expand the data in list_2d if it's too short
	"""
	import gtk
	length_of_treeview = len(treeview.get_columns())
	for ls in list_2d:
		data = ls[:]	#copy the list to avoid change the content in ls, 'data=ls' changes the content of ls
		for i in range(length_of_treeview-len(data)):
			data.append('')
		liststore.append(data)
	# set the TreeView mode to be liststore
	treeview.set_model(liststore)

	if reorderable:
		for i in range(len(list_2d[0])):
			# make it searchable
			treeview.set_search_column(i)
		
		# Allow drag and drop reordering of rows
		treeview.set_reorderable(True)
	#setting the selection mode
	treeselection = treeview.get_selection()
	treeselection.set_mode(gtk.SELECTION_MULTIPLE)

def get_no_of_total_genes(curs, schema=None, gene_table='gene'):
	"""
	04-18-05
		get the no_of_total_genes
	"""
	sys.stderr.write("Getting no_of_total_genes...")
	if schema:
		curs.execute("set search_path to %s"%schema)
	
	curs.execute("select count(gene_no) from %s"%gene_table)
	rows = curs.fetchall()
	for row in rows:
		no_of_total_genes = int(row[0])
	sys.stderr.write("Done\n")
	return no_of_total_genes

def get_edge_vector_by_id(curs, edge_id_list, edge_table='edge_cor_vector'):
	"""
	04-25-05
		return cor_2d_list and sig_2d_list given a list of edge id's 
	"""
	cor_2d_list = []
	sig_2d_list = []
	for edge_id in edge_id_list:
		curs.execute("select cor_vector, sig_vector from %s where edge_id=%s"%(edge_table,edge_id))
		rows = curs.fetchall()
		if len(rows) == 0:
			sys.stderr.write('%s not found in %s\n'%(edge_id, edge_table))
			continue
		cor_vector = rows[0][0][1:-1].split(',')
		cor_vector = map(float, cor_vector)
		cor_2d_list.append(cor_vector)
		
		sig_vector = rows[0][1][1:-1].split(',')
		sig_vector = map(int, sig_vector)
		sig_2d_list.append(sig_vector)
	return (cor_2d_list, sig_2d_list)

def get_edge_vector_by_tuple(curs, edge_set, edge_table='edge_cor_vector'):
	"""
	04-25-05
		return cor_2d_list and sig_2d_list given a list of edge tuples, which
		is like [1,7]
	"""
	cor_2d_list = []
	sig_2d_list = []
	for edge in edge_set:
		edge_string = '{' + repr(edge)[1:-1] + '}'
		curs.execute("select cor_vector, sig_vector from %s where edge_name='%s'"%(edge_table,edge_string))
		rows = curs.fetchall()
		if len(rows) == 0:
			sys.stderr.write('%s not found in %s\n'%(edge_string, edge_table))
			continue
		cor_vector = rows[0][0][1:-1].split(',')
		cor_vector = map(float, cor_vector)
		cor_2d_list.append(cor_vector)
		
		sig_vector = rows[0][1][1:-1].split(',')
		sig_vector = map(int, sig_vector)
		sig_2d_list.append(sig_vector)
	return (cor_2d_list, sig_2d_list)

def get_vertex_edge_list_by_edge_id(curs, edge_id_list, edge_table='edge_cor_vector'):
	"""
	04-25-05
		return vertex_list and edge_list given a list of edge id's 
	"""
	vertex_set = Set()
	edge_list = []
	for edge_id in edge_id_list:
		curs.execute("select edge_name from %s where edge_id=%s"%(edge_table,edge_id))
		rows = curs.fetchall()
		if len(rows) == 0:
			sys.stderr.write('%s not found in %s\n'%(edge_id, edge_table))
			continue
		edge = rows[0][0][1:-1].split(',')
		#edge = map(int, edge)	no need to cast to int. do it later
		edge_list.append(edge)
		vertex_set.add(edge[0])
		vertex_set.add(edge[1])

	return (list(vertex_set), edge_list)

def system_call(commandline):
	"""
	05-16-05
		call external program recursively based on the exit_code
		exit_code non-zero means error, needs calling again
	"""
	import os
	exit_code = os.system(commandline)
	while exit_code:
		exit_code = os.system(commandline)
	return exit_code


def mpi_synchronize(communicator):
	"""
	05-19-05
		copied from MpiBiclustering.py
	"""
	import sys
	sys.stdout.flush()
	sys.stderr.flush()
	communicator.barrier()
	
def mpi_schedule_jobs(communicator, job_list, node_function, node_parameter_list, debug=0):
	"""
	05-19-05
		a universal scheduling function, the elements in job_list
		maybe string, integer, or something else, It's 'repr'ed before send.
		WARNING: NO -1 in job_list.
		So node_function should handle it as 'repr'ed.
		
		node_function()
			input: (element of the job_list, node_parameter_list).
			output: returns a value('repr'ed)
		
		
	"""
	import sys
	from sets import Set
	node_returned_value_list = []
	node_rank = communicator.rank
	if node_rank == 0:
		sys.stderr.write("\tTotally, %d jobs to be scheduled.\n"%len(job_list))
		seed_utilized = Set()
		for node in range(1, communicator.size):
			if len(job_list)==0:	#if #nodes > #jobs, tell those nodes to break their listening loop.
				stop_signal = "-1"
				communicator.send(stop_signal, node, 0)	#no more jobs, stop that node,
				if debug:
					sys.stderr.write("node %s stopped.\n"%node)
			else:
				job = job_list.pop(0)	#the first item poped first.
				communicator.send(repr(job), node, 0)	#string format
				if debug:
					sys.stderr.write("node %s schedule a job, %s to %s\n"%(node_rank, repr(job), node))
				seed_utilized.add(node)
		
		received_value, source, tag = communicator.receiveString(None, None)	#listen
		while received_value:
			node_returned_value_list.append(received_value)
			if len(job_list) == 0:	#first check if there're still files left, otherwise pop(0) raises error.
				stop_signal = "-1"
				communicator.send(stop_signal, source, 0)	#no more jobs, stop that node,
				if debug:
					sys.stderr.write("node %s stopped.\n"%source)
				seed_utilized.remove(source)
				if len(seed_utilized) == 0:	#all seed used have finished their jobs
					break
			else:
				job = job_list.pop(0)
				communicator.send(repr(job), source, 0)	#string format,
				if debug:
					sys.stderr.write("node %s get one more job, %s\n"%(source, repr(job)) )
			received_value, source, tag = communicator.receiveString(None, None)	#listen
	else:
		received_data, source, tag = communicator.receiveString(0, None)	#get data from node 0,
		while received_data:
			if received_data=="-1":	#stop signal
				if debug:
					sys.stderr.write("node %s breaked.\n"%node_rank)
				break
			else:
				sys.stderr.write("node %s working on %s...\n"%(node_rank, received_data))
				node_return_value = node_function(received_data, node_parameter_list)
				sys.stderr.write("node %s work on %s finished.\n"%(node_rank, received_data))
				communicator.send(repr(node_return_value), 0, node_rank)
				
			received_data, source, tag = communicator.receiveString(0, None)	#get data from node 0
	
	return node_returned_value_list

def graphDotOutput(output_f, graph, label_dict, gene_no2go_no, centralnode=-1, function=0, functioncolor='green', \
	plot_type='dot', weighted=1, unknown_color='red'):
	'''
	06-14-05
		output graph in graphviz dot language.(similar to subgraph_output() of subgraph_visualize.py)
	10-12-05
		default centralnode=-1(no gene_no is -1)
		default unknown_color=red
	'''
	output_f.write('graph G  {\n')
	output_f.write('\tnode [shape=ellipse, color=black];\n')
	vertex_set = graph.node_list()
	vertex_labels = []
	for vertex in vertex_set:
		label = label_dict[vertex]
		shape = 'ellipse'
		color = 'black'
		if vertex in gene_no2go_no:
			if function in gene_no2go_no[vertex]:
				color = functioncolor
			if gene_no2go_no[vertex]==[0]:
				color = unknown_color
		if vertex==centralnode:
			shape = 'box'
			color = 'yellow'
		output_f.write('\t%s [label="%s", shape=%s, color=%s];\n'%(vertex, label, shape, color))
		
	for edge_id,edge_tuple in graph.edges.iteritems():
		if weighted:
			output_f.write('\t%s -- %s [label="%s"];\n'%(edge_tuple[0], edge_tuple[1], edge_tuple[2]))
		else:
			output_f.write('\t%s -- %s ;\n'%(edge_tuple[0], edge_tuple[1]))
	output_f.write('}\n')

def org2tax_id(organism):
	"""
	09-08-05
		return tax_id, given organism
	"""
	tax_id = None
	org2tax_id = {'Arabidopsis thaliana':3702,
		'Caenorhabditis elegans':6239,
		'Drosophila melanogaster':7227,
		'Homo sapiens':9606,
		'Mus musculus':10090,
		'Saccharomyces cerevisiae':4932,
		'Rattus norvegicus':10116}
	if organism in org2tax_id:
		tax_id = org2tax_id[organism]
	return tax_id

def org_short2long(organism):
	"""
	09-08-05
		return long-form organism
	"""
	long_organism = None
	org_short2long = {'at':'Arabidopsis thaliana',
		'ce':'Caenorhabditis elegans',
		'dm':'Drosophila melanogaster',
		'hs':'Homo sapiens',
		'mm':'Mus musculus',
		'sc':'Saccharomyces cerevisiae',
		'rn':'Rattus norvegicus',
		'Rattus norvegicus':'Rattus norvegicus',
		'Arabidopsis thaliana':'Arabidopsis thaliana',
		'Caenorhabditis elegans':'Caenorhabditis elegans',
		'Drosophila melanogaster':'Drosophila melanogaster',
		'Homo sapiens':'Homo sapiens',
		'Mus musculus':'Mus musculus',
		'Gorilla gorilla Pan paniscus Homo sapiens':'Homo sapiens',
		'Saccharomyces cerevisiae':'Saccharomyces cerevisiae'}
	if organism in org_short2long:
		long_organism = org_short2long[organism]
	return long_organism

def get_mt_id2no(curs, table='transfac.matrix'):
	"""
	09-19-05
	09-30-05
		modify it to fit for tables in different schema
	"""
	sys.stderr.write("Getting mt_id2no...")
	mt_id2no = {}
	curs.execute("select mt_id, id from %s"%table)
	rows = curs.fetchall()
	for row in rows:
		mt_id2no[row[0]] = row[1]
	sys.stderr.write("Done\n")
	return mt_id2no
	
def get_gene_id2mt_no_list(tax_id, hostname='zhoudb', dbname='graphdb', schema='graph', table='gene_id2mt_no'):
	"""
	09-19-05
	"""
	sys.stderr.write("Getting gene_id2mt_no_list...")
	gene_id2mt_no_list = {}
	(conn, curs) =  db_connect(hostname, dbname, schema)
	curs.execute("select gene_id, mt_no from %s where tax_id=%s"%(table, tax_id))
	rows = curs.fetchall()
	for row in rows:
		gene_id, mt_no = row
		if gene_id not in gene_id2mt_no_list:
			gene_id2mt_no_list[gene_id] = []
		gene_id2mt_no_list[gene_id].append(mt_no)
	del conn, curs
	sys.stderr.write("Done\n")
	return gene_id2mt_no_list

def get_global_gene_id2gene_no(curs, organism, table='graph.gene_id_to_no'):
	"""
	09-19-05
	"""
	sys.stderr.write("Getting global gene_id2gene_no...")
	gene_id2gene_no = {}
	curs.execute("select gene_id, gene_no from %s where organism='%s'"%(table, organism))
	rows = curs.fetchall()
	for row in rows:
		gene_id2gene_no[row[0]] = row[1]
	sys.stderr.write("Done.\n")
	return gene_id2gene_no
	
def get_mt_no2tf_name(hostname='zhoudb', dbname='graphdb', schema='transfac', table='matrix'):
	"""
	09-26-05
	"""
	sys.stderr.write("Getting mt_no2tf_name...")
	mt_no2tf_name = {}
	(conn, curs) =  db_connect(hostname, dbname, schema)
	curs.execute("select id,tf_name from %s"%table)
	rows = curs.fetchall()
	for row in rows:
		mt_no2tf_name[row[0]] = row[1]
	del conn, curs
	sys.stderr.write("Done\n")
	return mt_no2tf_name

def get_mcl_id2tf_set(curs, table, mt_no2tf_name):
	"""
	09-26-05
	09-27-05
		only score_type=1(hypergeometric), and take the score out
	"""
	sys.stderr.write("Getting mcl_id2tf_set...")
	mcl_id2tf_set = {}
	curs.execute("select mcl_id, bs_no_list, score, local_ratio from %s \
		where bs_no_list is not null and score_type=1"%table)	#score_type==1 only hypergeometric
	rows = curs.fetchall()
	for row in rows:
		mcl_id, bs_no_list, score, local_ratio = row
		bs_no_list = bs_no_list[1:-1].split(',')
		bs_no_list = map(int, bs_no_list)
		if local_ratio>0:
			tf_name_tuple = tuple(dict_map(mt_no2tf_name, bs_no_list))
			ratio_tuple = tuple([score])
			value_tuple = tuple([tf_name_tuple, ratio_tuple])
			if mcl_id not in mcl_id2tf_set:
				mcl_id2tf_set[mcl_id] = Set([value_tuple])
			else:
				mcl_id2tf_set[mcl_id].add(value_tuple)
	sys.stderr.write("Done.\n")
	return mcl_id2tf_set

def get_gene_id2gene_symbol(curs, tax_id, table='gene.gene'):
	"""
	09-28-05
		NCBI Gene id to Gene symbol
	"""
	sys.stderr.write("Getting gene_id2gene_symbol...")
	gene_id2gene_symbol = {}
	curs.execute("select gene_id, gene_symbol from %s where tax_id=%s"%(table, tax_id))
	rows = curs.fetchall()
	for row in rows:
		gene_id = row[0]
		gene_id2gene_symbol[gene_id] = row[1]
	sys.stderr.write("Done.\n")
	return gene_id2gene_symbol

def dict_transfer(dict1, dict2):
	"""
	09-28-05
		transfer the mapping between dict1 and dict2
	"""
	dict3 = {}
	for key1, value1  in dict1.iteritems():
		if value1 in dict2:
			dict3[key1] = dict2[value1]
	return dict3

"""
09-29-05
	class to represent tables related to a schema
"""
class schema:
	def __init__(self):
		self.splat_table = None
		self.mcl_table = None
		self.p_gene_table = None
		self.lm_table = None
		self.good_p_gene_table = None
		self.gene_p_table = None
		self.good_cluster_table = None
		self.cluster_bs_table = None
		self.prediction_suffix = None
		self.lm_suffix = None

"""
09-29-05
	function to form the table names of a schema given ofname and acc_cut_off
10-10-05
	add d_matrix_table
10-12-05
	add lm_bit
10-14-05
	add  pattern_table
"""
def form_schema_tables(ofname, acc_cut_off=0.6, lm_bit='111'):
	schema_instance = schema()
	schema_instance.pattern_table = 'pattern_%s'%ofname
	schema_instance.splat_table = 'splat_%s'%ofname
	schema_instance.mcl_table = 'mcl_%s'%ofname
	schema_instance.d_matrix_table = 'd_matrix_%s'%ofname
	schema_instance.prediction_suffix = '%s_e5'%ofname
	schema_instance.p_gene_table = 'p_gene_%s_e5'%ofname
	acc_int=int(acc_cut_off*100)
	if lm_bit=='111':
		lm_bit = ''	#backward compatibility
	schema_instance.lm_suffix = '%s_e5_%sa%s'%(ofname, lm_bit, acc_int)
	schema_instance.lm_table = 'lm_' + schema_instance.lm_suffix
	schema_instance.good_p_gene_table = 'p_gene_' + schema_instance.lm_suffix
	schema_instance.gene_p_table='gene_p_'  + schema_instance.lm_suffix
	schema_instance.good_cluster_table = 'good_cl_' + schema_instance.lm_suffix
	schema_instance.cluster_bs_table = 'cluster_bs_' + schema_instance.lm_suffix
	schema_instance.good_bs_table = 'good_bs_' + schema_instance.lm_suffix #11-01-05
	return schema_instance

"""
09-30-05
	mapping tax_id to unigene_prefix,
		like 9606(human) to Hs
"""
def tax_id2unigene_prefix(tax_id):
	tax_id2unigene_prefix = {3702:'At',
		6239:'Cel',
		7227:'Dm',
		9606:'Hs',
		10090:'Mm',
		10116:'Rn'}
	return tax_id2unigene_prefix.get(tax_id)

"""
09-30-05
	contruct unigene2gene_list mapping dict from gene2unigene file
"""
def get_unigene2gene_list(inputfile, tax_id):
	sys.stderr.write("Starting to get unigene2gene_list...")
	import csv
	unigene2gene_list = {}
	unigene_prefix = tax_id2unigene_prefix(tax_id)
	reader = csv.reader(open(inputfile,'r'), delimiter='\t')
	for row in reader:
		unigene_id = row[1]
		gene_id = row[0]
		if unigene_id.find(unigene_prefix)==0:
			if unigene_id not in unigene2gene_list:
				unigene2gene_list[unigene_id] = []
			unigene2gene_list[unigene_id].append(gene_id)
	sys.stderr.write("End to get unigene2gene_list.\n")
	return unigene2gene_list

"""
10-15-05
"""
def p_gene_id_set_from_gene_p_table(curs, gene_p_table):
	sys.stderr.write("Starting to get p_gene_id_set...")
	from sets import Set
	p_gene_id_set = Set()
	curs.execute("DECLARE crs CURSOR FOR select p_gene_id from %s"%gene_p_table)
	curs.execute("fetch 10000 from crs")
	rows = curs.fetchall()
	while rows:
		for row in rows:
			p_gene_id = row[0]
			p_gene_id_set.add(p_gene_id)
		curs.execute("fetch 10000 from crs")
		rows = curs.fetchall()
	curs.execute("close crs")
	sys.stderr.write("End to get p_gene_id_set.\n")
	return p_gene_id_set

def get_go_no2edge_counter_list(curs, gene_no2go_no_set, \
	edge_type2index={(0,0):0, (0,1):1, (1,0):1, (0,2):2, (2,0):2, (1,1):3, (1,2):4, (2,1):4, (2,2):5 }, \
	edge_table='edge_cor_vector', debug=0):
	"""
	10-20-05
		construct a edge_counter list which represents the counter of six types of edges for each function
	
	"""
	sys.stderr.write("Getting go_no2edge_counter_list...\n")
	go_no2edge_counter_list = {}	#value is a six integer list, each corresponds to below edge_type
	edge_type = [0,0]
	curs.execute("DECLARE crs CURSOR FOR select edge_name from %s "%(edge_table))
	curs.execute("fetch 5000 from crs")
	rows = curs.fetchall()
	no_of_total_edges = 0
	while rows:
		if debug:
			is_continue=raw_input("Continue?Y/n")
			if is_continue=='n':
				break
		for row in rows:
			edge_name = row[0][1:-1].split(',')
			edge_name = map(int, edge_name)
			gene_no1, gene_no2 = edge_name
			go_no_set1 = gene_no2go_no_set[gene_no1]
			go_no_set2 = gene_no2go_no_set[gene_no2]
			go_no_set = go_no_set1 | go_no_set2
			if debug:
				print "gene_no1", gene_no1, "go_no_set1", go_no_set1
				print "gene_no2", gene_no2, "go_no_set2", go_no_set2
				print "go_no_set", go_no_set
				is_continue=raw_input("Continue?Y/n")
			for go_no in go_no_set:
				if go_no not in go_no2edge_counter_list:
					go_no2edge_counter_list[go_no] = [0]*6
				
				if go_no in go_no_set1:
					edge_type[0] = 0
				elif 0 in go_no_set1:
					edge_type[0] = 2
				else:
					edge_type[0] = 1
				
				if go_no in go_no_set2:
					edge_type[1] = 0
				elif 0 in go_no_set2:
					edge_type[1] = 2
				else:
					edge_type[1] = 1
				go_no2edge_counter_list[go_no][edge_type2index[tuple(edge_type)]] += 1
				if debug:
					print "go_no", go_no, "edge_type", edge_type, "its edge_counter_list",go_no2edge_counter_list[go_no]
			no_of_total_edges += 1
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
	curs.execute("close crs")
	no_of_edges_half_unknown = go_no2edge_counter_list[0][1]
	if debug:
		print "no_of_total_edges",no_of_total_edges,"no_of_edges_half_unknown",no_of_edges_half_unknown
	for go_no in go_no2edge_counter_list:
		if go_no == 0:	#0 is different
			#no of edges, both are other function(1-1)
			go_no2edge_counter_list[go_no][3] = no_of_total_edges - sum(go_no2edge_counter_list[go_no])
			#0-2,1-2 are 0-0 are same
			go_no2edge_counter_list[go_no][2] = go_no2edge_counter_list[go_no][5] = go_no2edge_counter_list[go_no][0]
		else:
			#no of total unknown edges (2-2)
			go_no2edge_counter_list[go_no][5] = go_no2edge_counter_list[0][0]
			#no of edges, one is other function, one is unknown(1-2)
			go_no2edge_counter_list[go_no][4] = no_of_edges_half_unknown - go_no2edge_counter_list[go_no][2]
			#no of edges, both are other function(1-1)
			go_no2edge_counter_list[go_no][3] = no_of_total_edges - sum(go_no2edge_counter_list[go_no])
	sys.stderr.write("End getting go_no2edge_counter_list.\n")
	return go_no2edge_counter_list

def combine_numerator_a_denominator_dict(numerator_dict, denominator_dict):
	"""
	10-20-05
		based on the key in denominator_dict.
		If one key appears in numerator_dict, not in denominator_dict, consider its
		corresponding denominator is 1.
		If one key appears in denominator_dict, not in numerator_dict, the corresponding
		numerator is 0.
	"""
	final_value_dict = {}
	for key in denominator_dict:
		if key in numerator_dict:
			final_value_dict[key] = numerator_dict[key] / float(denominator_dict[key])
		else:
			final_value_dict[key] = 0.0
	return final_value_dict

def one_dim_list2string(ls):
	"""
	10-21-05
		for database submission
	"""
	ls_string = '{' + repr(ls)[1:-1] + '}'
	return ls_string

	
def get_no_of_unknown_genes(gene_no2go):
	"""
	10-20-05
	"""
	sys.stderr.write("Getting no_of_unknown_genes from gene_no2go...\n")
	no_of_unknown_genes = 0
	for gene_no, go_no_set in gene_no2go.iteritems():
		if 0 in go_no_set:
			no_of_unknown_genes += 1
	sys.stderr.write("Done getting no_of_unknown_genes.\n")
	return no_of_unknown_genes

def get_go_no2gene_no_set(curs, schema=None, go_table='go'):
	"""
	10-21-05
		get the go_no2gene_no_set
	"""
	sys.stderr.write("Getting go_no2gene_no_set...")
	if schema:
		curs.execute("set search_path to %s"%schema)
	go_no2gene_no_set = {}
	curs.execute("select go_no,gene_array from %s"%go_table)
	rows = curs.fetchall()
	for row in rows:
		go_no, gene_array = row
		gene_array = gene_array[1:-1].split(',')
		gene_array = map(int, gene_array)
		go_no2gene_no_set[go_no] = Set(gene_array)
	sys.stderr.write("Done\n")
	return go_no2gene_no_set

"""
10-21-05
	below output_node(), fetch_cluster_block(),input_node(), computing_node() are shared
	functions for Mpi Programs
"""

def output_node(communicator, free_computing_nodes, parameter_list, handler, report=0, type=1):
	"""
	10-20-05
	10-21-05
		handle the situation when free_computing_node list is exhausted
		add type to specify the communicating data type
		based on type, different ways to get stop_signal
	10-21-05
		a common output_node() function
		two jobs:
		1. give node 0 the free_computing_node
		2. handle the data from the computing_nodes
		
		handler(communicator, parameter_list, data)
	"""
	node_rank = communicator.rank
	sys.stderr.write("Node no.%s ready to accept output...\n"%node_rank)
	if type==1:
		data, source, tag = communicator.receiveString(None, 1)
	else:
		data, source, tag, count = communicator.receive(type, None, 1)
	no_of_computing_nodes = len(free_computing_nodes)
	no_of_resting_nodes = 0	#to keep track how many computing_nodes have rested
	free_computing_node_request = 0	#10-21-05 Flag used to check whether node 0 is requesting
	request_node = 0	#default request_node is 0
	while 1:
		if source==0:	#10-19-05 the input_node is asking me for free computing_node WATCH: it's array
			if len(free_computing_nodes)==0:	#10-21-05, it's empty. Flag the request.
				free_computing_node_request = 1
				request_node = source	#record the request_node
			else:
				free_computing_node = free_computing_nodes.pop(0)
				communicator.send(str(free_computing_node), source, 2)	#WATCH tag is 2.
		elif type==1 and data=="-1":
			no_of_resting_nodes += 1
			if report:
				sys.stderr.write("node %s(%s-th) rested.\n"%(source, no_of_resting_nodes))
			if no_of_resting_nodes==no_of_computing_nodes:	#WATCH: its' size-3
				break
				if report:
					sys.stderr.write("node %s output finished.\n"%node_rank)
		elif type!=1 and data.toscalar()==-1:	#10-21-05 for non-string type
			no_of_resting_nodes += 1
			if report:
				sys.stderr.write("node %s(%s-th) rested.\n"%(source, no_of_resting_nodes))
			if no_of_resting_nodes==no_of_computing_nodes:	#WATCH: its' size-3
				break
				if report:
					sys.stderr.write("node %s output finished.\n"%node_rank)
		else:
			if free_computing_node_request==1:	#there's a request for free_computing_node.
				communicator.send(str(source), request_node, 2)	#WATCH tag is 2.
				free_computing_node_request = 0
			else:
				free_computing_nodes.append(source)	#append the free computing_node
			handler(communicator, parameter_list, data)
		if type==1:
			data, source, tag = communicator.receiveString(None, 1)
		else:
			data, source, tag, count = communicator.receive(type, None, 1)
	sys.stderr.write("Node no.%s output done.\n"%node_rank)

def fetch_cluster_block(curs, message_size, report=0):
	"""
	10-20-05
		commonly used by input_node()
	"""
	if report:
		sys.stderr.write("Fetching stuff...\n")
	curs.execute("fetch 50 from crs")
	rows = curs.fetchall()
	prediction_ls = []
	string_length = 0
	while rows:
		for row in rows:
			prediction_ls.append(list(row))
			string_length += len(repr(row))	#the length to control MPI message size
		if string_length>=message_size:
			break
		curs.execute("fetch 50 from crs")
		rows = curs.fetchall()
	if report:
		sys.stderr.write("Fetching done.\n")
	return prediction_ls
	

def input_node(communicator, curs, free_computing_nodes, message_size, report=0, input_handler=fetch_cluster_block):
	"""
	10-20-05
	10-22-05
		add input_handler and regard curs as parameter_list
	"""
	node_rank = communicator.rank
	sys.stderr.write("Input node(%s) working...\n"%node_rank)
	data = input_handler(curs, message_size, report)
	counter = 0
	while data:
		communicator.send("1", communicator.size-1, 1)	#WATCH: tag is 1, to the output_node.
		free_computing_node, source, tag = communicator.receiveString(communicator.size-1, 2)
			#WATCH: tag is 2, from the output_node
		data_pickle = cPickle.dumps(data, -1)
		communicator.send(data_pickle, int(free_computing_node), 0)	#WATCH: int()
		if report:
			sys.stderr.write("block %s sent to %s.\n"%(counter, free_computing_node))
		data = input_handler(curs, message_size, report)
		counter += 1
	#tell computing_node to exit the loop
	for node in free_computing_nodes:	#send it to the computing_node
		communicator.send("-1", node, 0)
	sys.stderr.write("Input node(%s) done\n"%(node_rank))

def computing_cleanup_handler(communicator):
	"""
	10-22-05
		default cleanup_handler
	"""
	#tell the last node to stop
	communicator.send("-1", communicator.size-1, 1)	#tag is 1

def computing_node(communicator, parameter_list, node_fire_handler, cleanup_handler=computing_cleanup_handler, report=0):
	"""
	10-21-05
		0 is the node where the data is from
		size-1 is the node where output goes
	10-22-05
		make computing_cleanup_handler to be the default cleanup_handler
	"""
	node_rank = communicator.rank
	data, source, tag = communicator.receiveString(0, 0)	#get data from node 0
	while 1:
		if data=="-1":
			if report:
				sys.stderr.write("node %s breaked.\n"%node_rank)
			break
		else:
			result = node_fire_handler(communicator, data, parameter_list)
			result_pickle = cPickle.dumps(result, -1)
			communicator.send(result_pickle, communicator.size-1, 1)
		data, source, tag = communicator.receiveString(0, 0)	#get data from node 0
	cleanup_handler(communicator)

"""
10-22-05
	a function to get edge2occurrence
"""
def get_edge2occurrence(curs, min_sup=0, max_sup=200, edge_table='edge_cor_vector'):
	sys.stderr.write("Getting edge2occurrence...\n")
	curs.execute("DECLARE crs CURSOR FOR select edge_name,sig_vector \
		from %s"%(edge_table))
	edge2occurrence = {}
	curs.execute("fetch 5000 from crs")
	rows = curs.fetchall()
	no_of_datasets = 0
	while rows:
		for row in rows:
			edge = row[0][1:-1].split(',')
			edge = map(int, edge)
			sig_vector = row[1][1:-1].split(',')
			sig_vector = map(int, sig_vector)
			if no_of_datasets==0:
				no_of_datasets = len(sig_vector)
			if sum(sig_vector)>=min_sup and sum(sig_vector)<=max_sup:
				edge2occurrence[tuple(edge)] = sum(sig_vector)
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
	curs.execute("close crs")
	sys.stderr.write("Got edge2occurrence.\n")
	return edge2occurrence, no_of_datasets

"""
10-31-05 basic functions to draw images
"""
import Image, ImageDraw
def get_char_dimension():
	im = Image.new('RGB', (50,50))
	draw = ImageDraw.Draw(im)
	char_dimension = draw.textsize('a')
	del im, draw
	return char_dimension


def get_text_region(text, dimension, rotate=1, foreground=(0,0,255), background=(255,255,255)):
	"""
	10-31-05 add background and foreground
	"""
	text_im = Image.new('RGB', dimension, background)
	text_draw = ImageDraw.Draw(text_im)
	text_draw.text((0,0), text, fill=foreground)
	box = (0,0,dimension[0], dimension[1])
	text_reg = text_im.crop(box)
	if rotate:
		text_reg = text_reg.transpose(Image.ROTATE_90)	#90 is anti-clockwise
	return text_reg

"""
10-31-05
"""
def p_gene_id_src_set_from_gene_p_table(curs, gene_p_table):
	sys.stderr.write("Starting to get p_gene_id_src_set...")
	from sets import Set
	p_gene_id_src_set = Set()
	curs.execute("DECLARE crs CURSOR FOR select p_gene_id_src from %s"%gene_p_table)
	curs.execute("fetch 10000 from crs")
	rows = curs.fetchall()
	while rows:
		for row in rows:
			p_gene_id_src = row[0]
			p_gene_id_src_set.add(p_gene_id_src)
		curs.execute("fetch 10000 from crs")
		rows = curs.fetchall()
	curs.execute("close crs")
	sys.stderr.write("End to get p_gene_id_src_set.\n")
	return p_gene_id_src_set

"""
10-31-05
"""
def get_gene_no2p_go_no_set_given_p_gene_id_set(curs, p_gene_table, given_p_gene_id_set, report=0):
	sys.stderr.write("Getting gene_no2p_go_no_set...\n")
	gene_no2p_go_no_set = {}
	curs.execute("DECLARE crs CURSOR FOR SELECT p_gene_id, gene_no, go_no \
		from %s"%p_gene_table)
	curs.execute("fetch 10000 from crs")
	rows = curs.fetchall()
	counter = 0
	real_counter = 0
	while rows:
		for row in rows:
			p_gene_id, gene_no, go_no = row
			if p_gene_id in given_p_gene_id_set:
				if gene_no not in gene_no2p_go_no_set:
					gene_no2p_go_no_set[gene_no] = Set()
				gene_no2p_go_no_set[gene_no].add(go_no)
				real_counter += 1
			counter += 1
		if report:
			sys.stderr.write("%s%s/%s"%('\x08'*20, counter, real_counter))
		curs.execute("fetch 10000 from crs")
		rows = curs.fetchall()
	curs.execute("close crs")
	return gene_no2p_go_no_set

"""
11-01-05
"""
def cluster_bs_id_set_from_good_bs_table(curs, good_bs_table):
	sys.stderr.write("Getting cluster_bs_id_set...\n")
	cluster_bs_id_set = Set()
	curs.execute("select cluster_bs_id from %s"%good_bs_table)
	rows = curs.fetchall()
	for row in rows:
		cluster_bs_id_set.add(row[0])
	sys.stderr.write("End getting cluster_bs_id_set.\n")
	return cluster_bs_id_set

"""
11-01-05
11-29-05
	modify it to read ensembl2no_of_events_table and ensembl_id2gene_id_table
"""
def get_gene_no2no_of_events(curs, ensembl2no_of_events_table='graph.ensembl2no_of_events', \
	ensembl_id2gene_id_table='graph.ensembl_id2gene_id', source_table='graph.gene_id2no_of_events', gene_table='gene', schema=None):
	sys.stderr.write("Getting gene_no2no_of_events...\n")
	gene_no2no_of_events = {}
	curs.execute("SELECT g.gene_no,avg(e2.no_of_events) as no_of_events from %s e1, %s e2, %s g \
		where e1.ensembl_id=e2.ensembl_id and e1.gene_id=g.gene_id group by g.gene_no"%(ensembl_id2gene_id_table,\
		ensembl2no_of_events_table, gene_table))
	#curs.execute("select distinct g.gene_no, e1.no_of_events from %s g, %s e1, %s e2 where g.gene_id=e.gene_id"%\
		#(gene_table, source_table))
	rows = curs.fetchall()
	for row in rows:
		gene_no, no_of_events = row
		gene_no2no_of_events[gene_no] = no_of_events
	sys.stderr.write("End getting gene_no2no_of_events.\n")
	return gene_no2no_of_events


"""
11-29-05
	if one gene belongs to >1 families, take the bigger family size
"""
def get_gene_no2family_size(curs, tax_id, ensembl_id2family_id_table='graph.ensembl_id2family_id', ensembl_id2gene_id_table='graph.ensembl_id2gene_id', \
	gene_table='gene', schema=None):
	sys.stderr.write("Getting gene_no2family_size...\n")
	#1st, get the family_id2size
	family_id2size = {}
	curs.execute("select count(ensembl_id), family_id from %s where tax_id=%s group by family_id"%\
		(ensembl_id2family_id_table, tax_id))
	rows = curs.fetchall()
	for row in rows:
		size, family_id = row
		family_id2size[family_id] = size
	#2nd, get the gene_no2family_size
	gene_no2family_size = {}
	curs.execute("select distinct g.gene_no, ef.family_id from %s ef,  %s eg, %s g where ef.ensembl_id=eg.ensembl_id\
		and eg.gene_id=g.gene_id"%(ensembl_id2family_id_table, ensembl_id2gene_id_table, gene_table))
	rows = curs.fetchall()
	for row in rows:
		gene_no, family_id = row
		if gene_no not in gene_no2family_size:
			gene_no2family_size[gene_no] = 0
		if family_id2size[family_id] > gene_no2family_size[gene_no]:	#always take the larger one
			gene_no2family_size[gene_no]  = family_id2size[family_id]
	sys.stderr.write("End getting gene_no2family_size.\n")
	return gene_no2family_size


"""
11-03-05
	Draw a pattern via graphviz commands.
	1st call graphDotOutput() to output the input file.
"""
def draw_pattern(vertex_set, edge_set, go_no, gene_no2gene_id, gene_no2go_no, cluster_info_instance, graphFname,\
	centralnode=-1, plot_type_command='neato -Goverlap=false'):
	subgraph = cluster_info_instance.graph_from_node_edge_set(vertex_set, edge_set)
	graphSrcFname = '%s.dot'%graphFname
	graphSrcF = open(graphSrcFname, 'w')
	graphDotOutput(graphSrcF, subgraph, gene_no2gene_id, gene_no2go_no, centralnode=centralnode,\
		function=go_no, weighted=0)
	graphSrcF.close()
	commandline = '%s -Tpng %s -o %s'%(plot_type_command, graphSrcFname, graphFname)
	system_call(commandline)
	os.remove(graphSrcFname)


def cal_hg_p_value(gene_no, go_no, vertex_list, no_of_total_genes, go_no2gene_no_set, r, debug=0):
	"""
	11-10-05 'from rpy import r' should be executed in upper context
	11-18-05 fix a bug in " n=no_of_total_genes - cluster_size"
	"""
	cluster_size = len(vertex_list)
	no_of_local_associated_genes = 0
	for vertex in vertex_list:
		if vertex in go_no2gene_no_set[go_no] and vertex!=gene_no:	#leave one out
			no_of_local_associated_genes += 1
	no_of_global_associated_genes = len(go_no2gene_no_set[go_no])
	if gene_no in go_no2gene_no_set[go_no]:	#leave one out
		no_of_global_associated_genes -= 1
	x = no_of_local_associated_genes
	m = no_of_global_associated_genes
	n = no_of_total_genes - m
	k = cluster_size
	"""
	if debug:
		print "no_of_local_associated_genes",no_of_local_associated_genes
		print "no_of_global_associated_genes",no_of_global_associated_genes
		print "no_of_total_genes",no_of_total_genes
		print "cluster_size", cluster_size
	"""
	p_value = r.phyper(x-1,m,n,k,lower_tail = r.FALSE)
	return p_value

"""
11-12-05 get a specified segment from chromosome sequence table
	reverse is handled but complement(strand) is not. Upper level function should take care of this.
11-15-05 improve it to be more robust, add acc_ver and report if not found in raw_sequence_table
"""
def get_sequence_segment(curs, gi, start, stop, annot_assembly_table='sequence.annot_assembly', \
	raw_sequence_table='sequence.raw_sequence', chunk_size=10000):
	need_reverse = int(start>stop)
	if need_reverse:
		start, stop = stop, start
	curs.execute("select acc_ver, start, stop, raw_sequence_start_id from %s where gi=%s"%(annot_assembly_table, gi))
	rows = curs.fetchall()
	acc_ver, orig_start, orig_stop, raw_sequence_start_id = rows[0]
	if stop>orig_stop:	#11-14-05 to avoid exceeding the boundary
		stop = orig_stop
	no_of_chunks_before = start/chunk_size	#how many chunks are before this segment
	segment_size = stop - start +1
	no_of_chunks_segment = segment_size/chunk_size + 1	#how many chunks could include this segment
	raw_sequence_start_id += no_of_chunks_before	#the first chunk which contains teh segment
	offset = no_of_chunks_segment + 1	#add one more chunk to ensure the segment is within
	#get the sequence from raw_sequence_table
	seq = ''
	for i in range(offset):
		curs.execute("select sequence from %s where acc_ver='%s' and id=%s"%\
			(raw_sequence_table, acc_ver, raw_sequence_start_id+i))
		rows = curs.fetchall()
		if rows:	#11/14/05 it's possible to exceed the whole raw_sequence table because the offset adds one more chunk
			seq += rows[0][0]
		else:
			sys.stderr.write("id %s missing in raw_sequence table.\n"%(raw_sequence_start_id+i))
			sys.stderr.write("gi: %s, start: %s, stop: %s, raw_sequence_start_id: %s\n"%(gi, start, stop, raw_sequence_start_id))
	relative_start = start - no_of_chunks_before*chunk_size
	segment = seq[relative_start-1:relative_start-1+segment_size]	#WATCH index needs -1
	if need_reverse:
		segment = list(segment)
		segment.reverse()	#only 
		segment = ''.join(segment)
	return segment

"""
11-13-05
"""
def pg_1d_array2python_ls(pg_array, type_code=int):
	pg_array = pg_array[1:-1].split(',')
	pg_array = map(type_code, pg_array)
	return pg_array

def tax_id2org(tax_id):
	"""
	11-14-05 reverse of org2tax_id()
	"""
	organism = None
	tax_id2org = {3702:'Arabidopsis thaliana',
		6239:'Caenorhabditis elegans',
		7227:'Drosophila melanogaster',
		9606:'Homo sapiens',
		10090:'Mus musculus',
		4932:'Saccharomyces cerevisiae',
		10116:'Rattus norvegicus'}
	if tax_id in tax_id2org:
		organism = tax_id2org[tax_id]
	return organism


"""
11-15-05
	return a mt_id_set from a profile_filename
"""
def get_mt_id_set_from_profile(profile_filename):
	sys.stderr.write("Getting mt_id_set from %s..."%profile_filename)
	mt_id_set = Set()
	profile_inf = open(profile_filename, 'r')
	#skip the first four lines
	for i in range(4):
		profile_inf.readline()
	for line in profile_inf:
		line_ls = line.split()
		if len(line_ls)==5:	#to avoid the last line
			mt_id = line_ls[4]
			mt_id_set.add(mt_id)
	del profile_inf
	sys.stderr.write("Done.\n")
	return mt_id_set

"""
11-15-05
	the 1st bit of sites_ls determines whether it's consensus(0) or real aligned sites(1)
"""
def get_mt_id2sites_ls(curs, mt_id_set=None, matrix_table='transfac.matrix', site_in_matrix_table='transfac.site_in_matrix'):
	sys.stderr.write("Getting mt_id2sites_ls ...\n")
	mt_id2sites_ls = {}
	curs.execute("DECLARE mt_crs CURSOR FOR SELECT mt_id, mt_acc, consensus, site_in_matrix_accs from %s"%matrix_table)
	curs.execute("fetch 1000 from mt_crs")
	rows = curs.fetchall()
	while rows:
		for row in rows:
			mt_id, mt_acc, consensus, site_in_matrix_accs = row
			if mt_id_set:	#if mt_id_set is available
				if mt_id not in mt_id_set:
					continue
			if site_in_matrix_accs:
				mt_id2sites_ls[mt_id] = [1]
				curs.execute("select sequence from %s where mt_acc='%s'"%(site_in_matrix_table, mt_acc))
				sequence_rows = curs.fetchall()
				for sequence_row in sequence_rows:
					mt_id2sites_ls[mt_id].append(sequence_row[0].upper())
			elif consensus:
				mt_id2sites_ls[mt_id] = [0, consensus.upper()]
			else:
				sys.stderr.write("mt_id %s has no site_in_matrix_accs and consensus.\n")
		curs.execute("fetch 1000 from mt_crs")
		rows = curs.fetchall()
	curs.execute("close mt_crs")
	sys.stderr.write("Done getting mt_id2sites_ls.\n")
	return mt_id2sites_ls

"""
11-28-05
	draw grid to an Image object.
"""
def draw_grid(image_object, draw_object, region_to_draw, x_gap, y_gap, color='black'):
	start_x, start_y, stop_x, stop_y = region_to_draw
	#draw horizontal 1st
	for i in range(start_y, stop_y, y_gap):
		draw_object.line((start_x, i+y_gap, stop_x, i+y_gap), fill=color)
	#draw vertical 2nd
	for i in range(start_x, stop_x, x_gap):
		draw_object.line((i+x_gap, start_y, i+x_gap, stop_y), fill=color)


"""
11-30-05
	prerequisite: from sets import Set

	given the pattern id and gene_no, return the neighbor genes within the given distance
"""
def get_neighbor_set(curs, pattern_id, pattern_table, gene_no, distance):
	curs.execute("SELECT vertex_set, d_matrix from %s where id=%s"%(pattern_table, pattern_id))
	rows = curs.fetchall()
	vertex_set, d_matrix = rows[0]
	vertex_set = vertex_set[1:-1].split(',')
	vertex_set = map(int, vertex_set)
	gene_index = vertex_set.index(gene_no)
	d_matrix = d_matrix[2:-2].split('},{')
	d_row = d_matrix[gene_index]
	d_row = d_row.split(',')
	d_row = map(int, d_row)
	neighbor_set = Set()
	for i in range(len(d_row)):
		if d_row[i] <= distance and d_row[i]>0:	#not itself
			neighbor_set.add(vertex_set[i])
	return neighbor_set


"""
11-30-05
	
	Each neighbor_set, represented by a go_no is a vertex. Two vertices are connected by an edge if
	the two neighbor_sets are similar enough.
"""
from graph.cc_from_edge_list import cc_from_edge_list
def distinct_go_no_list_based_on_neighbor_set_graph(neighbor_set_ls, go_no_ls, similarity_cutoff, debug=0):
	no_of_go_nos = len(go_no_ls)
	edge_list = []
	go_no_in_graph_flag_ls = [0]*no_of_go_nos	#corresponding to go_no_ls, a flag of whether go_no is in the graph or not
	for i in range(no_of_go_nos):
		"""
		if debug:
			sys.stderr.write("Examining go_no %s with neighbor_set %s.\n"%(go_no_ls[i], neighbor_set_ls[i]))
		"""
		for j in range(i+1, no_of_go_nos):
			intersection_set = neighbor_set_ls[i]&neighbor_set_ls[j]
			intersection_set_len = float(len(intersection_set))
			if intersection_set_len/len(neighbor_set_ls[i])>=similarity_cutoff or intersection_set_len/len(neighbor_set_ls[j])>=similarity_cutoff:
				"""
				if debug:
					sys.stderr.write("\t go_no %s with neighbor_set %s connected to it.\n"%(go_no_ls[j], neighbor_set_ls[j]))
				"""
				#either one of them is very close to the intersection_set
				edge_list.append([go_no_ls[i], go_no_ls[j]])
				#these two go_nos are in the graph
				go_no_in_graph_flag_ls[i] = 1
				go_no_in_graph_flag_ls[j] = 1
	#handle the singleton go_no first
	singleton_go_no_list = []
	for i in range(no_of_go_nos):
		if go_no_in_graph_flag_ls[i] ==0:
			singleton_go_no_list.append(go_no_ls[i])
	#check the edge_list
	if len(edge_list)>0:
		cf_instance = cc_from_edge_list()
		cf_instance.run(edge_list)
		cc_list = cf_instance.cc_list
	else:
		cc_list = []
	
	return cc_list, singleton_go_no_list

"""
11-30-05

	--get_neighbor_set()
	--distinct_go_no_list_based_on_neighbor_set_graph()
"""
def get_gene_no2no_of_topologies(curs, schema_instance, similarity_cutoff, distance, debug=0):
	sys.stderr.write("Getting gene_no2no_of_topologies...\n")
	curs.execute("DECLARE gncrs  CURSOR FOR SELECT p.p_gene_id, p.gene_no, p.go_no, p.mcl_id from %s p, %s g\
		where g.p_gene_id=p.p_gene_id"%\
		(schema_instance.p_gene_table, schema_instance.gene_p_table))
	curs.execute("fetch 5000 from gncrs")
	rows = curs.fetchall()
	gene_no2mcl_id_ls = {}
	while rows:
		for row in rows:
			p_gene_id, gene_no, go_no, mcl_id = row
			if gene_no not in gene_no2mcl_id_ls:
				gene_no2mcl_id_ls[gene_no] = []
			gene_no2mcl_id_ls[gene_no].append(mcl_id)
		curs.execute("fetch 5000 from gncrs")
		rows = curs.fetchall()
	curs.execute("close gncrs")
	
	gene_no2no_of_topologies = {}
	for gene_no, mcl_id_list in gene_no2mcl_id_ls.iteritems():
		neighbor_set_ls = []
		for mcl_id in mcl_id_list:
			neighbor_set_ls.append(get_neighbor_set(curs, mcl_id, schema_instance.pattern_table, gene_no, distance))
		cc_list, singleton_go_no_list = distinct_go_no_list_based_on_neighbor_set_graph(neighbor_set_ls, mcl_id_list, similarity_cutoff, debug)
		gene_no2no_of_topologies[gene_no] = len(cc_list)+len(singleton_go_no_list)
		"""
		if debug:
			sys.stderr.write("gene_no: %s, with %s topologies.\n"%(gene_no, len(cc_list)+len(singleton_go_no_list)))
			raw_input("continue?(Y/n)")
		"""
	sys.stderr.write("Done getting gene_no2no_of_topologies.\n")
	return gene_no2no_of_topologies

"""
12-06-05
"""
def get_gene_no2incidence_array(gim_inputfname, gene_id2no):
	sys.stderr.write("Getting gene_no2incidence_array...\n")
	reader = csv.reader(open(gim_inputfname), delimiter='\t')
	gene_no2incidence_array = {}
	for row in reader:
		no_of_occurrences, occ_array, gene_id = row[0], row[1:-1], row[-1]
		if gene_id in gene_id2no:
			gene_no = gene_id2no[gene_id]
			occ_array = map(int, occ_array)
			gene_no2incidence_array[gene_no] = occ_array
	sys.stderr.write("End getting gene_no2incidence_array.\n")
	return gene_no2incidence_array

"""
12-06-05
"""
def get_vertex_set_gim_array(gene_no2incidence_array, vertex_set):
	gim_array = [0]*len(gene_no2incidence_array[vertex_set[0]])
	total_no_of_genes = float(len(vertex_set))
	for i in range(len(gim_array)):
		no_of_incident_genes = 0
		for gene_no in vertex_set:
			if gene_no2incidence_array[gene_no][i] == 1:
				no_of_incident_genes += 1
		gim_array[i] = no_of_incident_genes/total_no_of_genes
	return gim_array

"""
12-07-05
	get the edge2encodedOccurrence from the database
"""
def fill_edge2encodedOccurrence(curs, min_sup, max_sup, edge_table='edge_cor_vector'):
	sys.stderr.write("Getting edge2encodedOccurrence...\n")
	sys.path += [os.path.expanduser('~/script/annot/bin')]
	from MpiFromDatasetSignatureToPattern import encodeOccurrenceBv
	edge2encodedOccurrence = {}
	curs.execute("DECLARE crs CURSOR FOR select edge_name,sig_vector \
		from %s"%(edge_table))
	curs.execute("fetch 5000 from crs")
	rows = curs.fetchall()
	no_of_datasets = 0
	counter = 0
	while rows:
		for row in rows:
			edge = row[0][1:-1].split(',')
			edge = map(int, edge)
			sig_vector = row[1][1:-1].split(',')
			sig_vector = map(int, sig_vector)
			if no_of_datasets==0:
				no_of_datasets = len(sig_vector)
			if sum(sig_vector)>=min_sup and sum(sig_vector)<=max_sup:
				edge2encodedOccurrence[tuple(edge)] = encodeOccurrenceBv(sig_vector)
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
	curs.execute("close crs")
	sys.stderr.write("Done.\n")
	return edge2encodedOccurrence, no_of_datasets

def get_gene_symbol2gene_id(curs, tax_id, table='gene.gene_symbol2id'):
	"""
	12-09-05
		reverse of get_gene_id2gene_symbol()
	"""
	sys.stderr.write("Getting gene_symbol2gene_id...")
	gene_symbol2gene_id = {}
	curs.execute("select gene_id, gene_symbol from %s where tax_id=%s"%(table, tax_id))
	rows = curs.fetchall()
	for row in rows:
		gene_symbol2gene_id[row[1]] = row[0]
	sys.stderr.write("Done.\n")
	return gene_symbol2gene_id


def get_ensembl_id2gene_id(curs, tax_id, table='graph.ensembl_id2gene_id'):
	"""
	12-11-05
	"""
	sys.stderr.write("Getting ensembl_id2gene_id...")
	ensembl_id2gene_id = {}
	curs.execute("select ensembl_id, gene_id from %s where tax_id=%s"%(table, tax_id))
	rows = curs.fetchall()
	for row in rows:
		ensembl_id, gene_id = row
		gene_id = int(gene_id)	#convert it to integer
		if ensembl_id not in ensembl_id2gene_id:
			ensembl_id2gene_id[ensembl_id] = []
		ensembl_id2gene_id[ensembl_id].append(gene_id)
	sys.stderr.write("Done.\n")
	return ensembl_id2gene_id

def get_entrezgene_coord(curs, gene_id, entrezgene_mapping_table='sequence.entrezgene_mapping',\
	annot_assembly_table='sequence.annot_assembly'):
	"""
	12-11-05
	"""
	curs.execute("select e.genomic_gi, a.chromosome, e.strand, e.start, e.stop from %s e, %s a \
		where e.genomic_gi = a.gi and e.gene_id=%s"%(entrezgene_mapping_table, \
		annot_assembly_table, gene_id))
	rows = curs.fetchall()
	if rows:
		genomic_gi, chromosome, strand, start, stop = rows[0]
		if strand=='1':
			strand = '+'
		else:
			strand = '-'
		return genomic_gi, chromosome, strand, start, stop
	else:
		return None

"""
12-12-05
"""
def get_gene_no2no_of_components_given_gene_id_set(curs, gene_id2gene_no, organism, table='graph.cc_association'):
	sys.stderr.write("Getting gene_no2no_of_components...")
	gene_no2no_of_components = {}
	curs.execute("SELECT go_id, gene_id from %s where organism='%s'"%(table, organism))
	rows = curs.fetchall()
	for row in rows:
		go_id, gene_id = row
		if gene_id in gene_id2gene_no:
			gene_no = gene_id2gene_no[gene_id]
			if gene_no not in gene_no2no_of_components:
				gene_no2no_of_components[gene_no] = 0
			gene_no2no_of_components[gene_no] += 1
	sys.stderr.write("Done.\n")
	return gene_no2no_of_components

"""
12-12-05
"""
def get_tg_tax_id2ca_depth(curs, src_tax_id, tax_id_relationship_table = 'homologene.tax_id_relationship'):
	sys.stderr.write("Getting tg_tax_id2ca_depth ...")
	tg_tax_id2ca_depth = {}
	curs.execute("SELECT tg_tax_id, common_ancestor_depth from %s where src_tax_id=%s"%\
		(tax_id_relationship_table, src_tax_id))
	rows = curs.fetchall()
	for row in rows:
		tg_tax_id, common_ancestor_depth = row
		tg_tax_id2ca_depth[tg_tax_id] = common_ancestor_depth
	sys.stderr.write("Done.\n")
	return tg_tax_id2ca_depth

"""
12-12-05
"""
def get_gene_id2ca_depth(curs, src_tax_id, tg_tax_id2ca_depth, homologene_table='homologene.homologene'):
	sys.stderr.write("Getting gene_id2ca_depth ...")
	gene_id2ca_depth = {}
	curs.execute("SELECT hid, tax_id, gene_id from %s order by hid"%homologene_table)
	rows = curs.fetchall()
	prev_hid = None
	prev_src_gene_id_ls = []
	prev_depth = 100
	for row in rows:
		hid, tax_id, gene_id = row
		if hid!=prev_hid:	#not the 1st row, prev_src_gene_id_ls not empty
			for gene_id in prev_src_gene_id_ls:
				gene_id2ca_depth[gene_id] = prev_depth
			prev_hid = hid
			prev_src_gene_id_ls = []
			prev_depth = 100
		if tax_id == src_tax_id:
			prev_src_gene_id_ls.append(gene_id)
		prev_depth = min(tg_tax_id2ca_depth[tax_id], prev_depth)
	#don't forget the last block
	for gene_id in prev_src_gene_id_ls:
		gene_id2ca_depth[gene_id] = prev_depth
	sys.stderr.write("Done.\n")
	return gene_id2ca_depth

"""
12-12-05
"""
def get_gene_no2no_of_tfbs(curs, schema='harbison2004'):
	sys.stderr.write("Getting gene_no2no_of_tfbs ... ")
	gene_no2no_of_tfbs = {}
	curs.execute("DECLARE motif_crs CURSOR FOR SELECT mt_id, prom_id from %s.binding_site"%\
		(schema))
	curs.execute("fetch 5000 from motif_crs")
	rows = curs.fetchall()
	while rows:
		for row in rows:
			mt_id, prom_id = row
			prom_id = int(prom_id)
			if prom_id not in gene_no2no_of_tfbs:
				gene_no2no_of_tfbs[prom_id] = 0
			gene_no2no_of_tfbs[prom_id] += 1
		curs.execute("fetch 5000 from motif_crs")
		rows = curs.fetchall()
	curs.execute("close motif_crs")
	sys.stderr.write("Done.\n")
	return gene_no2no_of_tfbs
