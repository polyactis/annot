#!/usr/bin/env python
"""
Usage:	DrawGONodeStrucGraph.py -k xx -i xx -g xx -o xx [OPTIONS}

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-p ..., p_gene table
	-n ...,	gene_p table
	-g ...,	GO informative node file, output of 7-go/go_informative_node.py
	-o ...,	output picture fname
	-y ...,	type of running, 1(go_id2gene_set from -g, default), 2(from -p and -n)
	-e ...,	sharing_threshold for an edge between two go nodes, 15(default)
	-u,	debug
	-r,	report the progress(a number)
	-h, --help      show this help
	
Examples:
	DrawGONodeStrucGraph.py -k hs_fim_65 -g /tmp/yuhuang/hs_fim_65_n2s125.go -o tmp/go_node_struc.n2s125.png
	DrawGONodeStrucGraph.py -k hs_fim_65 -p p_gene_hs_fim_65_n2s200_m5x65s4l5_ft2_e5 -n gene_p_hs_fim_65_n2s200_m5x65s4l5_ft2_e5_000001a60 -y 2 -o tmp/go_node_struc.n2s200.pred.png
	
Description:
	draw a graph , with nodes as GO node,  size proportional to how many genes it has,
	edge due to sharing of genes between GO nodes, edge width proportional to how many
	genes are shared
	

"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))

from codense.common import db_connect, get_go_id2name
import matplotlib; matplotlib.use("Agg")
import networkx as nx
import pylab
import csv, getopt
from sets import Set

class DrawGONodeStrucGraph:
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, p_gene_table=None, \
		gene_p_table=None, go_fname=None, fig_output_fname=None, running_type=1,\
		sharing_threshold=15, debug=0, report=0):
		"""
		2006-11-15
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.p_gene_table = p_gene_table
		self.gene_p_table = gene_p_table
		self.go_fname = go_fname
		self.fig_output_fname = fig_output_fname
		self.running_type = int(running_type)
		self.sharing_threshold = int(sharing_threshold)
		self.debug = int(debug)
		self.report = int(report)


	def draw_GO_node_structure_graph(self, go_id2gene_set, fig_output_fname, sharing_threshold=30,\
		label_map=None, is_to_show=1):
		"""
		2006-09-12
			draw a graph , with nodes as GO node,  size proportional to how many genes it has,
			edge due to sharing of genes between GO nodes, edge width proportional to how many
			genes are shared
			
			call get_go_id2gene_set() from upstairs
		"""
		sys.stderr.write("Drawing GO node structure graph...")
		g = nx.XGraph()
		go_id_list = go_id2gene_set.keys()
		no_of_go_ids = len(go_id_list)
		for i in range(no_of_go_ids):
			for j in range(i+1, no_of_go_ids):
				go_id1 = go_id_list[i]
				go_id2 = go_id_list[j]
				no_of_sharing_genes = len(go_id2gene_set[go_id1]&go_id2gene_set[go_id2])
				if no_of_sharing_genes>=sharing_threshold:
					g.add_edge(go_id1, go_id2, no_of_sharing_genes)
		edge_width_list = []
		for (u, v, d) in g.edges():
			edge_width_list.append(d/5)
		node_size_list = []
		for v in g:
			node_size_list.append(len(go_id2gene_set[v])*4)
		pos=nx.graphviz_layout(g)
		pylab.figure(figsize=(20,20))
		nx.draw_networkx_edges(g,pos,
			alpha=0.3,
			width=edge_width_list,
			edge_color='m')
		nx.draw_networkx_nodes(g,pos,
			node_size=node_size_list,
			node_color='r',
			alpha=0.4)
		nx.draw_networkx_edges(g,pos,
			alpha=0.4,
			node_size=0,
			width=1,
			edge_color='k')
		abr_label_map = {}
		if label_map:
			#don't give a dictionary with keys more than the nodes, (it'll get error)
			for key in g:
				if key in label_map:
					abr_label_map[key] = label_map[key]
		if abr_label_map:
			nx.draw_networkx_labels(g, pos, abr_label_map, fontsize=8)
		else:
			nx.draw_networkx_labels(g, pos, fontsize=8)
		#give a title
		pylab.title(os.path.basename(fig_output_fname))
		#write the legend
		xmin, xmax, ymin, ymax = pylab.axis()
		dx = xmax - xmin
		dy = ymax - ymin
		x = 0.02*dx + xmin
		y = 0.95*dy + ymin
		pylab.text(x, y, "edge cutoff: %s"%sharing_threshold)
		pylab.savefig(fig_output_fname)
		if is_to_show:
			pylab.show()
		sys.stderr.write("Done.\n")
	
	def get_informative_node2gene_set(self, go_fname):
		"""
		2006-11-15
			it's haifeng_annot/7-go/go_informative_node.py's output format
		"""
		sys.stderr.write("Getting go_id2gene_set...")
		reader = csv.reader(open(go_fname), delimiter='\t')
		go_id2gene_set = {}
		go_id_row = reader.next()[2:]
		go_index2id = dict(zip(range(len(go_id_row)), go_id_row))
		for row in reader:
			gene_id = int(row[0])
			if row[1]=='1':	#it's known gene
				for i in range(2, len(row)):
					if row[i]=='1':
						go_id = go_index2id[i-2]
						if go_id not in go_id2gene_set:
							go_id2gene_set[go_id] = Set()
						go_id2gene_set[go_id].add(gene_id)
		del reader
		sys.stderr.write("Done.\n")
		return go_id2gene_set
	
	def get_go_id2gene_set_from_prediction_cut_file(self, prediction_cut_fname):
		"""
		2006-09-13
			prediction_cut_fname is output of PredictionAccCurve.py
		2006-11-16
			not used, copied from misc.py with other functions
		"""
		sys.stderr.write("Getting go_id2gene_set...")
		reader = csv.reader(open(prediction_cut_fname), delimiter='\t')
		go_id2gene_set  = {}
		for row in reader:
			vertex_id, pattern_id, go_id, is_known, is_correct, m1_score, p_value = row
			if go_id not in go_id2gene_set:
				go_id2gene_set[go_id] = Set()
			go_id2gene_set[go_id].add(int(vertex_id))
		del reader
		sys.stderr.write("Done.\n")
		return go_id2gene_set
	
	def get_go_id2gene_set_from_prediction_table(self, curs, p_gene_table, gene_p_table):
		"""
		2006-11-15
			
			
		"""
		sys.stderr.write("Getting go_id2gene_set from db...")
		go_id2gene_set  = {}
		curs.execute("DECLARE crs CURSOR for select go.go_id, p.gene_no from %s p, %s g, go\
			where p.p_gene_id=g.p_gene_id and go.go_no=p.go_no"%(p_gene_table,\
			gene_p_table))
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		counter = 0
		while rows:
			for row in rows:
				go_id, gene_no = row
				if go_id not in go_id2gene_set:
					go_id2gene_set[go_id] = Set()
				go_id2gene_set[go_id].add(gene_no)
				counter +=1
			#sys.stderr.write('%s%s'%('\x08'*20, counter))
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		curs.execute("close crs")
		sys.stderr.write("Done.\n")
		return go_id2gene_set
	
	def run(self):
		"""
		2006-11-15
		"""
		conn, curs = db_connect(self.hostname, self.dbname, self.schema)
		if self.running_type==1:
			go_id2gene_set = self.get_informative_node2gene_set(self.go_fname)
		elif self.running_type==2:
			go_id2gene_set = self.get_go_id2gene_set_from_prediction_table(curs, self.p_gene_table, self.gene_p_table)
		
		go_id2name = get_go_id2name(curs)
		self.draw_GO_node_structure_graph(go_id2gene_set, self.fig_output_fname,\
			self.sharing_threshold, label_map=go_id2name, is_to_show=0)
		
		
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:p:n:g:o:y:e:ur", ["help", "hostname=", \
		"dbname=", "schema="])
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	p_gene_table = None
	gene_p_table = None
	go_fname = None
	fig_output_fname = None
	running_type = 1
	sharing_threshold = 15
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
		elif opt in ("-p",):
			p_gene_table = arg
		elif opt in ("-n",):
			gene_p_table = arg
		elif opt in ("-g",):
			go_fname = arg
		elif opt in ("-o",):
			fig_output_fname = arg
		elif opt in ("-y",):
			running_type = int(arg)
		elif opt in ("-e",):
			sharing_threshold = int(arg)
		elif opt in ("-u",):
			debug = 1
		elif opt in ("-r",):
			report = 1
	if schema and ((p_gene_table and gene_p_table) or go_fname) and fig_output_fname:
		instance = DrawGONodeStrucGraph(hostname, dbname, schema, p_gene_table, gene_p_table,\
			go_fname, fig_output_fname, running_type, sharing_threshold, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)