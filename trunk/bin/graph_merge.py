#!/usr/bin/env python
"""
Usage: graph_merge.py [OPTION] INPUTDIR OUTPUTFILE

Option:
	INPUTDIR is the directory containing all the graph files in gspan format.
	OUTPUTFILE is the file to store the merged graphs also in gspan format.
	-s ..., --support=...	minimum support for the edge to be kept. 5(default)
	-h, --help              show this help
	
Examples:
	graph_merge.py -s 6 gph_result/sc/ gph_result/sc_mcl/mcl_gph_dataset1

Description:
	This program merges all the graphs in gspan format, generated by
	graph_reorganize.py. The ouput is also in gspan format.
	After this, either run reverse+kMax or gspan2mcl_input.py+mcl to
	get the dense clusters.

"""


import sys, os, re, getopt
from codense.common import db_connect

class graph_merge:
	'''
	02-26-05
		use the kjDict to reduce the memory usage, but no effect.
	03-08-05
		use a database table to handle 
	'''
	def __init__(self, support, dir, ofname):
		self.support = int(support)
		self.dir = dir
		self.of = open(ofname, 'w')
		
		self.hostname = 'zhoudb'
		self.dbname = 'template1'
		self.schema = 'public'
		self.counter_table = 'counter_table'
	
	def dstruc_loadin(self, dir):
		#output block put before edges
		first_block = ''
		#data structure to store the merged graph. key is the edge.
		#value is the recurrence
		from kjbuckets import kjDict
		graph_dict = kjDict()
		
		files = os.listdir(dir)
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
		for f in files:
			pathname = os.path.join(dir, f)
			sys.stderr.write("%d/%d:\t%s\n"%(files.index(f)+1,len(files),f))
			file_no = files.index(f)+1
			inf = open(pathname, 'r')
			for line in inf:
				if line[0] == 'e':
					#edge here, like 'e 3807 3859 0.804645'
					line_list = line[:-1].split()
					vertex1 = int(line_list[1])
					vertex2 = int(line_list[2])
					if vertex1 <= vertex2:
						edge = (vertex1, vertex2)
					else:
						edge = (vertex2, vertex1)
					if graph_dict.has_key(edge):
						graph_dict[edge] += 1
					else:
						graph_dict[edge] = 1
				elif file_no == 1:
					first_block += line
			inf.close()
		
		return (first_block, graph_dict)

	def output(self, first_block, graph_dict, support):
		#output the preceding block first
		self.of.write(first_block)
		for edge in graph_dict.keys():
			recurrence = graph_dict[edge]
			if recurrence >= support:
				self.of.write("e %d %d %d\n"%(edge[0], edge[1], recurrence))

	def loadin_edges(self, dir, curs, counter_table):
		#output block put before edges
		first_block = ''

		files = os.listdir(dir)
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
		for f in files:
			pathname = os.path.join(dir, f)
			sys.stderr.write("%d/%d:\t%s\n"%(files.index(f)+1,len(files),f))
			file_no = files.index(f)+1
			inf = open(pathname, 'r')
			for line in inf:
				if line[0] == 'e':
					#edge here, like 'e 3807 3859 0.804645'
					line_list = line[:-1].split()
					vertex1 = int(line_list[1])
					vertex2 = int(line_list[2])
					edge = [vertex1, vertex2]
					self.add_one_edge(curs, counter_table, edge)
				elif file_no == 1:
					first_block += line
			inf.close()
		return first_block

	def create_counter_table(self, curs, counter_table):
		"""
		03-08-05
			create a counter table
		"""
		curs.execute("create table %s(\
			edge integer[] primary key,\
			frequency	integer)"%counter_table)
	
	def add_one_edge(self, curs, counter_table, edge):
		"""
		03-08-05
			add an edge to the table, if edge already exists, increase its frequency
		"""
		#keep in ascending order
		if edge[0]>edge[1]:
			edge.reverse()
		#format the edge
		edge = '{%s,%s}'%(edge[0], edge[1])
		curs.execute("select frequency from %s where edge='%s'"%(counter_table, edge))
		rows = curs.fetchall()
		if len(rows)==1:
			frequency = rows[0][0]
			#increase it by 1
			curs.execute("update %s set frequency=%d where edge='%s'"%(counter_table, frequency+1, edge))
		elif len(rows)==0:
			#insert it
			curs.execute("insert into %s(edge, frequency) values('%s', %d)"%(counter_table, edge, 1))
		else:
			sys.stderr.write("The edge, %s gets more than one frequency in the table.\n"%edge)
			sys.exit(1)
	
	def output_from_db(self, curs, counter_table, first_block, support):
		"""
		03-08-05
			output the merged graph by getting edges from database
		"""
		sys.stderr.write("Outputting merged graph...")
		#output the preceding block first
		self.of.write(first_block)
		
		curs.execute("DECLARE crs CURSOR FOR select edge, frequency \
			from %s"%(counter_table))
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		while rows:
			for row in rows:
				edge = row[0][1:-1].split(',')
				frequency = row[1]
				if frequency>=support:
					self.of.write("e %s %s %d\n"%(edge[0], edge[1], frequency))
			
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		sys.stderr.write("Done\n")
		
	
	def old_run(self):
		(first_block, graph_dict) = self.dstruc_loadin(self.dir)
		self.output(first_block, graph_dict, self.support)
	
	def run(self):
		"""
		03-08-05
			new run(), database approach to store the frequency of edges.
		
		--db_connect()
		--create_counter_table()
		--loadin_edges()
			--add_one_edge()
		--output_from_db()
		"""
		(conn, curs) = db_connect(self.hostname, self.dbname, self.schema)
		self.create_counter_table(curs, self.counter_table)
		first_block = self.loadin_edges(self.dir, curs, self.counter_table)
		self.output_from_db(curs, self.counter_table, first_block, self.support)
		
	
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "s:h", ["support=", "help"])
	except:
		print __doc__
		sys.exit(2)
	
	support = 5
	for opt, arg in opts:
		if opt in ("-s", "--support"):
			support = int(arg)
		elif opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)

			
	if len(args) == 2:
		instance = graph_merge(support, args[0], args[1])
		instance.run()
	else:
		print __doc__
		sys.exit(2)
