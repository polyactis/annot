#!/usr/bin/env python
"""
functions or classes common to all programs

"""

import csv, sys
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


def dict_map(dict, ls):
	new_list = []
	for item in ls:
		value = dict.get(item)
		if value:
			new_list.append(value)
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
	"""
	import gtk
	tvcolumn_dict = {}
	cell_dict = {}

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
	"""
	import gtk
	for data in list_2d:
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
