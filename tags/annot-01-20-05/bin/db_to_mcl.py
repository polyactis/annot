#!/usr/bin/env python
"""
Usage: db_to_mcl.py -k SCHEMA -s SIZE [OPTION] MCLINPUTDIR

Option:
	MCLINPUTDIR a directory to save the mcl input files.
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	the table storing the splat results, splat_result(default)
	-f ..., --format=...	output format, 0(default, mcl), 1(kMax)
	-s ..., --size=...	the number of splat records in one output file
	-h, --help              show this help
	
Examples:
	db_to_mcl.py -k shu -s 10000 gph_result/mcl_input
	
Description:
	Output schema.splat_result into several files. (mcl_input format)
	Each file is with the number of splat patterns.
	Output format:
		0,	MCL(not direct input of MCL. MCL wrapper program: mcl_batch_run.py)
		1,	kMax(also not direct input of kMax. kMax wrapper program)

"""

import sys, os, psycopg, getopt
from sets import Set

class db_to_mcl:
	def __init__(self, hostname, dbname, schema, table, format):
		self.vertex_dict = {}
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.table = table
		self.format = int(format)
		self.format_dict = {0:self.mcl_singleton_write,
			1:self.kMax_singleton_write}

	def get_from_db(self, outdir, package_size):
		if not os.path.isdir(outdir):
			os.makedirs(outdir)
		self.curs.execute("begin")
		self.curs.execute("DECLARE crs CURSOR FOR select splat_id,edge_set from %s"%self.table)
		self.curs.execute("fetch %d from crs"%(package_size))
		rows = self.curs.fetchall()
		no_of_splat_records = 0
		no_of_packages = 0
		while rows:
			outfname = 'mcl_input%d'%no_of_packages
			path = os.path.join(outdir, outfname)
			outf = open(path, 'w')
			for row in rows:
				self.format_dict[self.format](row, outf)
				no_of_splat_records += 1
			outf.close()
			no_of_packages += 1
			string = 'Packages: %d\t Splat_patterns: %d'%(no_of_packages, no_of_splat_records)
			sys.stderr.write('%s%s'%('\x08'*80,string))
			self.curs.execute("fetch %d from crs"%(package_size))
			rows = self.curs.fetchall()
		self.curs.execute("end")
		sys.stderr.write('\n\tTotal patterns: %d\n'%no_of_splat_records)
		
	def mcl_singleton_write(self, row, outf):
		self.vertex_dict = {}	#intialize the vertex_dict
		splat_id = row[0]
		edge_set = row[1]
		edge_list = edge_set[2:-2].split('},{')
		for edge in edge_list:
			vertex_list = edge.split(',')
			vertex1 = vertex_list[0]
			vertex2 = vertex_list[1]
			if vertex1 in self.vertex_dict:
				self.vertex_dict[vertex1].append(vertex2)
			else:
				self.vertex_dict[vertex1] = [vertex2]
			if vertex2 in self.vertex_dict:
				self.vertex_dict[vertex2].append(vertex1)
			else:
				self.vertex_dict[vertex2] = [vertex1]
		dim = len(self.vertex_dict)
		out_block = '(splat_id %s )\n'%splat_id		# here it is '=' not '+='
		out_block += '(mclheader\n'
		out_block += 'mcltype matrix\n'
		out_block += 'dimensions %dx%d\n)\n'%(dim,dim)
		out_block += '(mcldoms\n'
		vertex_list = self.vertex_dict.keys()
		out_block += '%s $\n)\n'%' '.join(vertex_list)
		out_block += '(mclmatrix\nbegin\n'
		for vertex in self.vertex_dict:
			out_block += '%s '%vertex
			out_block += '%s $\n'%' '.join(self.vertex_dict[vertex])
		out_block += ')\n\n'
		outf.write(out_block)
	
	def kMax_singleton_write(self, row, outf):
		vertex_set = Set()	#intialize the vertex_set
		splat_id = row[0]
		edge_set = row[1]
		edge_list = edge_set[2:-2].split('},{')
		for edge in edge_list:
			vertex_list = edge.split(',')
			vertex_set.add(int(vertex_list[0]))
			vertex_set.add(int(vertex_list[1]))
		outf.write('t %d\n'%splat_id)
		#cast to list type
		vertex_list = list(vertex_set)
		#sort it before output
		vertex_list.sort()
		for vertex in vertex_list:
			outf.write('v %d %d\n'%(vertex, vertex))
		for edge in edge_list:
			#replace the comma with a blank.
			edge = edge.replace(',', ' ')
			outf.write('e %s 0\n'%(edge))
		outf.write('\n')

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:f:s:", ["help", "hostname=", "dbname=", "schema=", "table=", "format=", "size="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = 'splat_result'
	format = 0
	size = 0
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
		elif opt in ("-f", "--format"):
			format = int(arg)
		elif opt in ("-s", "--size"):
			size = int(arg)

	if schema and size and len(args)==1:
		instance = db_to_mcl(hostname, dbname, schema, table, format)
		instance.get_from_db(args[0], size)

	else:
		print __doc__
		sys.exit(2)
