#!/usr/bin/env python
"""
Usage: mcl_result2fp_input.py -k SCHEMA [OPTION]

Option:
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-s ..., --size=...	the number of mcl records in one fetch, 5000(default)
	-n ..., --connectivity_cut_off=...	0.5(default), minimum connectivity of a mcl cluster
	-h, --help              show this help
	
Examples:
	mcl_result2fp_input.py -k sc_yh60_fp
	
Description:
	Output schema.mcl_result into two files. (fpall format)
	'dataset': transaction file
	'spec': dataset specification file

"""

import sys, os, psycopg, getopt

class mcl_result2fp_input:
	def __init__(self, dbname, schema, connectivity_cut_off):
		self.vertex_dict = {}
		self.conn = psycopg.connect('dbname=%s'%dbname)
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.connectivity_cut_off = float(connectivity_cut_off)
		self.no_of_vertices =0	#redundant
		self.max_cluster_size = 0
		self.no_of_clusters = 0
		self.datasetf = open('dataset', 'w')
		self.specf = open('spec', 'w')

	def get_from_db(self, package_size):
		self.curs.execute("begin")
		self.curs.execute("DECLARE crs CURSOR FOR select mcl_id, splat_id, vertex_set from mcl_result\
			where connectivity>=%f"%self.connectivity_cut_off)
		self.curs.execute("fetch %d from crs"%(package_size))
		rows = self.curs.fetchall()
		no_of_packages = 0
		while rows:
			for row in rows:
				#write the cluster into the 'dataset' file
				self.mcl_singleton_write(row)
				self.no_of_clusters += 1
			no_of_packages += 1
			string = 'Packages: %d\t Mcl_clusters: %d'%(no_of_packages, self.no_of_clusters)
			sys.stderr.write('%s%s'%('\x08'*80,string))
			self.curs.execute("fetch %d from crs"%(package_size))
			rows = self.curs.fetchall()
		self.curs.execute("end")
		#dataset specification output
		self.spec_write()
		sys.stderr.write('\n\tTotal patterns: %d\n'%self.no_of_clusters)
		
	def mcl_singleton_write(self, row):
		mcl_id = row[0]
		splat_id = row[1]
		vertex_set = row[2][1:-1]
		vertex_list = vertex_set.split(',')
		len_of_vertex_list = len(vertex_list)
		self.no_of_vertices += len_of_vertex_list
		if len_of_vertex_list > self.max_cluster_size:
			self.max_cluster_size = len_of_vertex_list
		#mcl's clusters are already sorted.
		self.datasetf.write('%s\n'%('\t'.join(vertex_list)))
		#self.datasetf.write('%d\t%s\n'%(len_of_vertex_list, '\t'.join(vertex_list)))
	
	def spec_write(self):
		self.specf.write('dataset\n')
		self.specf.write('%d\n'%self.no_of_vertices)
		self.specf.write('%d\n'%self.max_cluster_size)
		self.specf.write('%d\n'%self.no_of_clusters)
		
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hd:k:s:n:", ["help", "dbname=", "schema=", "size=", "connectivity_cut_off"])
	except:
		print __doc__
		sys.exit(2)
	
	dbname = 'graphdb'
	schema = ''
	size = 5000
	connectivity_cut_off = 0.5
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-k", "--schema"):
			schema = arg
		elif opt in ("-s", "--size"):
			size = int(arg)
		elif opt in ("-n", "--connectivity_cut_off"):
			connectivity_cut_off = float(arg)

	if schema:
		instance = mcl_result2fp_input(dbname, schema, connectivity_cut_off)
		instance.get_from_db(size)

	else:
		print __doc__
		sys.exit(2)
