#!/usr/bin/env python
"""
Usage: graph_reorganize.py -k SCHEMA [OPTION] DATADIR NEWDIR

Option:
	DATADIR is the directory of the graph construction results.
	NEWDIR is the directory to store the filtered results.
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --type=...	0(default, only known), 1(known+unknown) or 2(all)
			if 2 is selected, specify the organism
	-g ..., --organism=...	two letter organism abbreviation
	-h, --help              show this help
	
Examples:
	graph_reorganize.py -k shu normal/ known/
	graph_reorganize.py -k shu -t 2 -g sc normal/ known/	take all yeast genes

Description:
	This program depends on schema.gene and graph.gene_id_to_no.
	It does the filtration for the graph construction results.
	Three types of filtration:
	0: only known genes in that schema are preserved.
	1: only known genes and unknown genes in that schema are preserved.
	2: all the genes in table gene_id_to_no are preserved.

Note:
	There're some known genes that are not in schema.
	So 1 & 2 are different.
"""


import sys, pickle, os, psycopg, getopt

class vertex_attributes:
	"a class wrapping the attributes of a vertex"
	def __init__(self, no = None):
		self.freq = 1
		self.no = no
	
	def inc_freq(self):
		self.freq += 1

class graph_attributes:
	"a class wrapping the attributes of a graph"
	def __init__(self, no=None):
		self.no = no
		self.graph_dict = {}
		self.vertex_set = {} #key stores vertex's label, value is a list which stores the adjacent vertices.
	def vertex_set_init(self):
		for i in self.graph_dict:
			if self.vertex_set.has_key(i[0]):
				self.vertex_set[i[0]].append(i[1])
			if i[0] not in self.vertex_set:
				self.vertex_set[i[0]] = []
			if self.vertex_set.has_key(i[1]):
				self.vertex_set[i[1]].append(i[0])
			if i[1] not in self.vertex_set:
				self.vertex_set[i[1]] = []
				
class graph_reorganize_old:
	'''
	Read the output from graph_construct class(above), convert it to gSpan or Splat input format.
	Store the mapping of vertex_label v.s. number and graph_label v.s. number in file no_label_map in user's home directory.
	Based on this class, optimization of all the graphs is possible.
	
	!!!Defuncted. This class requires too much memory.
	Replaced by the class graph_reorganize.
	'''
	def __init__(self):
		self.vertex_list_dict = {}
		self.graph_list_dict = {}
		self.vertex_block=''		#a block of vertices in the output file. "v M L" means the Mth vertex with L label.
		
	def init(self):
		self.vertex_list_dict = {}
		self.graph_list_dict = {}
		self.vertex_block = ''
		
	def parse(self, inf):
		self.init()

		#loop below initilizes graph_list_dict.
		line = inf.readline()
		while line:
			list = line[:-1].split('\t')
			if line[0] == 't':
				graph_label = list[2]
				no_of_graphs = len(self.graph_list_dict) + 1
				self.graph_list_dict[graph_label] = graph_attributes(no=no_of_graphs)
			if line[0] == 'e':
				vertex1 = list[1]
				vertex2 = list[2]
				self.graph_list_dict[graph_label].graph_dict[(vertex1,vertex2)] = float(list[3])
			line = inf.readline()
		
		#loop below initilizes vertex_list_dict.
		for graph_label in self.graph_list_dict:
			graph_attr = self.graph_list_dict[graph_label]
			graph_attr.vertex_set_init()
			#initlize vertex_set of each graph

			for vertex in graph_attr.vertex_set:
				if self.vertex_list_dict.has_key(vertex):
					self.vertex_list_dict[vertex].inc_freq()
				else:
					no_of_vertices = len(self.vertex_list_dict) + 1
					self.vertex_list_dict[vertex] = vertex_attributes(no=no_of_vertices)
					self.vertex_block+='v %d %s\n'%(no_of_vertices,vertex,)
		self.output()
		
	def output(self):
		import os
		no_label_map_filename = os.path.join(os.path.expanduser('~'),'no_label_map')
		map_fhandler = open(no_label_map_filename, 'w')
		map_fhandler.write('>no_graphlabel_mapping\n')
		for graph_label in self.graph_list_dict:
			outf = open('splat_'+graph_label, 'w')
			graph_attr = self.graph_list_dict[graph_label]
			outf.write('t # %s\n'%graph_label)
			
			map_fhandler.write('%d\t%s\n'%(graph_attr.no, graph_label,))
			
			outf.write(self.vertex_block)

			for edge in graph_attr.graph_dict:
				vertex1_no = self.vertex_list_dict[edge[0]].no
				vertex2_no = self.vertex_list_dict[edge[1]].no
				outf.write('e %d %d %f\n'%(vertex1_no, vertex2_no, graph_attr.graph_dict[edge],))
			outf.close()

		map_fhandler.write('>no_vertexlabel_mapping\n')
		for vertex_label in self.vertex_list_dict:
			vertex_attr = self.vertex_list_dict[vertex_label]
			map_fhandler.write('%d\t%s\n'%(vertex_attr.no, vertex_label,))
			
class graph_reorganize:
	'''
	This class plays two major roles.
	Collect labels.
	Transform graph.cc results into gspan input.
	'''
	def __init__(self, dbname, schema, type, orgn):
		self.conn = psycopg.connect('dbname=%s'%dbname)
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.type = int(type)
		self.org_short2long = {'at':'Arabidopsis thaliana',
			'ce':'Caenorhabditis elegans',
			'dm':'Drosophila melanogaster',
			'hs':'Homo sapiens',
			'mm':'Mus musculus',
			'sc':'Saccharomyces cerevisiae',
			'Arabidopsis thaliana':'Arabidopsis thaliana',
			'Caenorhabditis elegans':'Caenorhabditis elegans',
			'Drosophila melanogaster':'Drosophila melanogaster',
			'Homo sapiens':'Homo sapiens',
			'Mus musculus':'Mus musculus',
			'Gorilla gorilla Pan paniscus Homo sapiens':'Homo sapiens',
			'Saccharomyces cerevisiae':'Saccharomyces cerevisiae'}
		if self.type == 2:
			self.organism = self.org_short2long[orgn]
			
		self.global_vertex_dict = {}
		self.global_graph_dict = {}
		self.pickle_fname = os.path.join(os.path.expanduser('~'), 'pickle/yeast_global_struc')
		self.vertex_block = ''
		
		#following four variables keeping track of the change
		self.global_vertex_dict_before = {}
		self.global_vertex_dict_after = {}
		self.no_of_edges_before = 0
		self.no_of_edges_after = 0
		
	def global_mapping_construct(self, inf):
		'''
		global_mapping_construct() collects the graph label and vertex label from inf.
		mapping_batch() will iteratively call this function to collect all labels pertaining
		to one organism.
		The collection will be stored in a file through pickle.
		'''
		line = inf.readline()
		while line:
			list = line[:-1].split()
			if line[0] == 't':
				graph_label = list[2]
				if graph_label in self.global_graph_dict:
					sys.stderr.write('Error. graph with this name, "%s" appears twice\n'%graph_label)
				self.global_graph_dict[graph_label]=1
			if line[0] == 'e':
				vertex1_label = list[1]
				vertex2_label = list[2]
				if vertex1_label not in self.global_vertex_dict:
					self.global_vertex_dict[vertex1_label] = 1
				if vertex2_label not in self.global_vertex_dict:
					self.global_vertex_dict[vertex2_label] = 1
			line = inf.readline()

	def dstruc_loadin(self):
		'''
		This method loads in the data structures from database.
		'''
		if self.type == 0:
			self.curs.execute("select gene_id, gene_no from gene where known=TRUE order by gene_no")
		if self.type == 1:
			self.curs.execute("select gene_id, gene_no from gene order by gene_no")
		if self.type == 2:
			self.curs.execute("select gene_id, gene_no from graph.gene_id_to_no \
				where organism='%s' order by gene_no"%self.organism)
			
		rows = self.curs.fetchall()
		for row in rows:
			self.global_vertex_dict[row[0]] = row[1]
			self.vertex_block += 'v %d %s\n'%(row[1], row[0])
		return 0
		
	def transform(self, inf, outf):
		'''
		This method transforms a graph.cc output file(inf) into a gspan input file(outf).
		It requires self.global_vertex_list and other data structures to be filled in first.
		So self.dstruc_loadin() must be called before this method.
		transform_batch() will iteratively call this method to transform all graph.cc output files.
		'''
		line = inf.readline()
		while line:
			list = line[:-1].split()
			if line[0] == 't':
				outf.write(line.replace('\t', ' '))
				outf.write(self.vertex_block)
			if line[0] == 'e':
				vertex1_label = list[1]
				vertex2_label = list[2]
				self.no_of_edges_before += 1
				self.global_vertex_dict_before[vertex1_label] = 1
				self.global_vertex_dict_before[vertex2_label] = 1
				if vertex1_label in self.global_vertex_dict and vertex2_label in self.global_vertex_dict:
					outf.write('e %d %d %s\n'% \
								(self.global_vertex_dict[vertex1_label], \
								self.global_vertex_dict[vertex2_label], \
								list[3],)
								)
					self.no_of_edges_after += 1
					self.global_vertex_dict_after[vertex1_label] = 1
					self.global_vertex_dict_after[vertex2_label] = 1
			line = inf.readline()
	
	def stat_output(self):
		sys.stderr.write("\tedges: %d(before), %d(after)\n"%(self.no_of_edges_before, self.no_of_edges_after))
		sys.stderr.write("\tvertices: %d(before), %d(after)\n"%(len(self.global_vertex_dict_before), len(self.global_vertex_dict_after)) )
	
def transform_batch(dbname, schema, type, orgn, dir, output_dir):
	'''
	See comments of graph_reorganize.transform().
	'''
	files = os.listdir(dir)
	sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
	if not os.path.isdir(output_dir):
		os.makedirs(output_dir)
	
	instance = graph_reorganize(dbname, schema, type, orgn)
	data_not_loaded = instance.dstruc_loadin()
	if data_not_loaded:
		return
	
	for f in files:
		pathname = os.path.join(dir, f)
		sys.stderr.write("%d/%d:\t%s\n"%(files.index(f)+1,len(files),f))
		inf = open(pathname, 'r')
		outf = open(os.path.join(output_dir,f), 'w')
		instance.transform(inf, outf)
		inf.close()
		outf.close()
	instance.stat_output()
	
def mapping_batch(dbname, schema, type, orgn, dir):
	'''
	See comments of graph_reorganize.global_mapping_construct().
	'''
	files = os.listdir(dir)
	sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
	instance = graph_reorganize(dbname, schema, type, orgn)	
	for f in files:
		pathname = os.path.join(dir, f)
		sys.stderr.write("%d/%d:\t%s\n"%(files.index(f)+1,len(files),f))
		inf = open(pathname, 'r')
		instance.global_mapping_construct(inf)
		inf.close()

	global_vertex_list = instance.global_vertex_dict.keys()
	global_vertex_list.sort()
	for i in range(len(global_vertex_list)):
		instance.global_vertex_dict[global_vertex_list[i]] = i+1
	global_struc = {'vertex_dict': instance.global_vertex_dict,
			'vertex_list': global_vertex_list,
			'graph_list': instance.global_graph_dict}
	pickle.dump(global_struc, open(instance.pickle_fname, 'w') )
	

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hd:k:t:g:", ["help", "dbname=", "schema=", "type=", "organism="])
	except:
		print __doc__
		sys.exit(2)
	
	dbname = 'graphdb'
	schema = ''
	type = 0
	organism = ''
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-k", "--schema"):
			schema = arg
		elif opt in ("-t", "--type"):
			type = int(arg)
		elif opt in ("-g", "-organism"):
			organism = arg
			
	if schema and len(args) == 2:
		if type == 2 and organism == '':
			print __doc__
			sys.exit(2)
		transform_batch(dbname, schema, type, organism, args[0], args[1])
	else:
		print __doc__
		sys.exit(2)
