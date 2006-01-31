#!/usr/bin/env python
"""
Usage: PermuteClusterGeneID.py -k -i -o

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-i ...,	input fname
	-o ...,	output fname
	-b,	enable debug flag
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	~/script/annot/bin/PermuteClusterGeneID.py -k hs_fim_40 -i -o 

Description:
	Program to permute the Gene Id of output of MpiBFSCluster.py
	 
"""

import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
import sys, os, getopt, csv, random
from codense.common import db_connect

class PermuteClusterGeneID:
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None,input_fname=None, output_fname=None, \
		debug=0, report=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_fname = input_fname
		self.output_fname = output_fname
		self.debug = int(debug)
		self.report = int(report)
	
	def get_gene_id_permutation_dict(self, curs, gene_table='gene'):
		sys.stderr.write("Getting gene_id_list...")
		curs.execute("select gene_no from %s"%gene_table)
		rows = curs.fetchall()
		gene_id_list = []
		for row in rows:
			gene_id_list.append(row[0])
		
		index_list = range(len(gene_id_list))
		random.shuffle(index_list)
		old_gene_id2new = {}
		for i in range(len(gene_id_list)):
			old_gene_id2new[gene_id_list[i]] = gene_id_list[index_list[i]]
		sys.stderr.write("Done.\n")
		return old_gene_id2new
	
	def permute_cluster_file(self, input_fname, output_fname, old_gene_id2new):
		sys.stderr.write("Permuting cluster file %s..."%os.path.basename(input_fname))
		reader = csv.reader(open(input_fname, 'r'), delimiter = '\t')
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		for row in reader:
			vertex_list, edge_list, recurrence_array, d_matrix = row
			vertex_list = vertex_list[1:-1].split(',')
			for i in range(len(vertex_list)):
				vertex_list[i] = old_gene_id2new[int(vertex_list[i])]
			vertex_list.sort()
			edge_list = edge_list[2:-2].split('], [')
			for i in range(len(edge_list)):
				edge_tuple = edge_list[i].split(',')				
				for j in range(len(edge_tuple)):
					edge_tuple[j] = old_gene_id2new[int(edge_tuple[j])]
				edge_tuple.sort()
				edge_list[i]  = edge_tuple
			edge_list.sort()
			writer.writerow([vertex_list, edge_list, recurrence_array, d_matrix])
		del reader, writer
		sys.stderr.write("Done.\n")
	
	def run(self):
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		old_gene_id2new  = self.get_gene_id_permutation_dict(curs)
		self.permute_cluster_file(self.input_fname, self.output_fname, old_gene_id2new)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:o:br", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	input_fname = None
	output_fname = None
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
		elif opt in ("-o",):
			output_fname = arg
		elif opt in ("-b",):
			debug = 1
		elif opt in ("-r",):
			report = 1
	if schema and input_fname and output_fname:
		instance = PermuteClusterGeneID(hostname, dbname, schema, input_fname, output_fname, \
			debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
