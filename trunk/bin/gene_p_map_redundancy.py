#!/usr/bin/env python
"""
Usage: gene_p_map_redundancy.py -k SCHEMA -t P_GENE_TABLE -n GENE_P_TABLE [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --p_gene_table=...	the p_gene table
	-n ..., --gene_p_table=...	update the gene_p table, needed if needcommit
	-c, --commit	commit this database transaction
	-r, --report	report flag
	-u, --debug debug flag
	-h, --help              show this help

Examples:
	gene_p_map_redundancy.py -k sc_54 -t p_gene_repos_2_e5
		-n gene_p_repos_2_e5
	
Description:
	This module merges the p_gene_ids whose predicted functions are
	parent-child, for each gene in the gene_p table
"""

import sys, os, getopt
from codense.common import db_connect, get_go_no2term_id

class gene_p_map_redundancy:
	"""
	03-01-05
		This module merges the p_gene_ids whose predicted functions are
		parent-child, for each gene in the gene_p table
	03-03-05
		fix an important bug
		

	run()
		--db_connect()
		--get_go_no2term_id()		#for get_distance(), needs self.go_no2term_id
		--data_fetch()
			--gene_no2p_gene_setup()
		--p_gene_id_map()
			--_p_gene_map()
				--get_distance()		#touches self.go_no2distance
		--submit()		
	
	"""
	def __init__(self, hostname=None, dbname=None, schema=None, p_gene_table=None, \
		gene_p_table=None, needcommit=0, report=0, debug=0):
		
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.p_gene_table = p_gene_table
		self.gene_p_table = gene_p_table
		self.needcommit = int(needcommit)
		self.report = int(report)		
		self.debug = int(debug)

		self.gene_no2p_gene = {}	#key is gene_no, value is p_gene_id2go_no
		self.p_gene_id_map = {}
		#mapping between a pair of go_no's and its associated distances
		self.go_no2distance = {}
		
		self.distance_table = 'go.node_dist'	#for get_distance()
		self.term_table = 'go.term'	#for get_go_no2term_id()
	
	def data_fetch(self, curs, p_gene_table, gene_p_table):
		"""
		03-01-05
			borrowed from p_gene_lm.py
			
			--gene_no2p_gene_setup()
		"""
		curs.execute("DECLARE crs CURSOR FOR select g.p_gene_id, p.gene_no, p.go_no\
			from %s g, %s p where g.p_gene_id=p.p_gene_id"%(gene_p_table, p_gene_table))
		no_of_records = 0
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		while rows:
			for row in rows:
				self.gene_no2p_gene_setup(row)
				no_of_records += 1
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, no_of_records))
			
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		if self.report:
			sys.stderr.write("\n")
	
	def gene_no2p_gene_setup(self, row):
		"""
		03-01-05
			initial
		"""
		p_gene_id = row[0]
		gene_no = row[1]
		go_no = row[2]
		if gene_no not in self.gene_no2p_gene:
			self.gene_no2p_gene[gene_no] = {}
		
		#pass the value to ease programming
		p_gene_id2go_no = self.gene_no2p_gene[gene_no]
		
		if p_gene_id not in p_gene_id2go_no:
			p_gene_id2go_no[p_gene_id] = go_no
		else:
			sys.stderr.write("Error, duplicate p_gene_id:%s found in gene_p_table\n"%p_gene_id)
			sys.exit(127)
			
	def p_gene_map(self, gene_no2p_gene, p_gene_id_map, curs, distance_table, go_no2distance, go_no2term_id):
		"""
		03-01-05
			initial
		"""
		sys.stderr.write("Computing the p_gene_id_map...")
		for (gene_no, p_gene_id2go_no) in gene_no2p_gene.iteritems():
			self._p_gene_map(p_gene_id2go_no, p_gene_id_map, curs, distance_table, \
				go_no2distance, go_no2term_id)
		sys.stderr.write("Done\n")
	
	def _p_gene_map(self, p_gene_id2go_no, p_gene_id_map, curs, distance_table, go_no2distance, go_no2term_id):
		"""
		03-01-05
			initial, modeled after return_distinct_functions() of gene_stat_plot.py
		"""
		p_gene_id_list = p_gene_id2go_no.keys()
		for i in range(len(p_gene_id_list)):
			p_gene_id1 = p_gene_id_list[i]
			if p_gene_id1 not in p_gene_id_map:
				#not mapped, first encounter, it's the source, mapp to itself
				p_gene_id_map[p_gene_id1] = p_gene_id1
				if self.debug:
					print "%s not in p_gene_id_map yet,mapped to itself"%p_gene_id1
				go_no1 = p_gene_id2go_no[p_gene_id1]
				for j in range(i+1, len(p_gene_id_list)):
					p_gene_id2 = p_gene_id_list[j]
					if p_gene_id2 not in p_gene_id_map:
						#the p_gene_id hasn't found its source
						go_no2 = p_gene_id2go_no[p_gene_id2]
						if go_no1 < go_no2:
							key= (go_no1, go_no2)
						else:
							key = (go_no2, go_no1)
						if key in go_no2distance:
							jasmine_distance = go_no2distance[key][2]	#03-03-05 fix an important bug here. [2] was missing
						else:
							jasmine_distance = self.get_distance(curs, go_no1,\
								go_no2, distance_table, go_no2distance, go_no2term_id)
						if self.debug:
							print "jasmine_distance of %s and %s is %s"%(go_no1, go_no2, jasmine_distance)
						if jasmine_distance == 0:
							#jasmine_distance=0 means they are parent-child
							p_gene_id_map[p_gene_id2] = p_gene_id1
							if self.debug:
								print "%s not in p_gene_id_map, mapped to %s"%(p_gene_id2, p_gene_id1)

					
			
			
	def get_distance(self, curs, go_no1, go_no2, distance_table, go_no2distance, go_no2term_id):
		"""
		03-01-05
			borrowed from gene_stat_plot.py
			this function is only envoked when the distance of (go_no1, go_no2) is not available.
			
			fill in self.go_no2distance
			return jasmine_distance
		"""
		if go_no1==go_no2:
			#this is possible because several p_gene_ids point to the same go_no.
			return 0
			
		go_id1 = go_no2term_id.get(go_no1)
		go_id2 = go_no2term_id.get(go_no2)
		
		if go_id1 > go_id2:
			#swap it
			go_id1, go_id2 = go_id2, go_id1
			
		
		if go_id1==None or go_id2==None:
			sys.stderr.write("go_no: %s or %s don't have term_id\n"%(go_no1, go_no2))
			sys.exit(13)
		
		curs.execute("select raw_distance, lee_distance, jasmine_distance, \
			common_ancestor_list from %s where go_id1=%d and go_id2=%d"%\
			(distance_table, go_id1, go_id2))
		rows = curs.fetchall()
		if len(rows) == 0:
			sys.stderr.write("go_id1: %s and go_id2: %s, distance not present\n"%(go_id1, go_id2))
			sys.exit(14)
		for row in rows:
			#key tuple in ascending order
			if go_no1<go_no2:
				go_no2distance[(go_no1, go_no2)] = (row[0], row[1], row[2])
			else:
				go_no2distance[(go_no2, go_no1)] = (row[0], row[1], row[2])
		return row[2]
	
	def submit(self, curs, gene_p_table, p_gene_id_map):
		"""
		03-01-05
			update the gene_p_table to reflect the mapping
		"""
		sys.stderr.write("Updating table %s..."%gene_p_table)
		for (p_gene_id, p_gene_id_src) in p_gene_id_map.iteritems():
			curs.execute("update %s set p_gene_id_src=%d where p_gene_id=%d"%\
				(gene_p_table, p_gene_id_src, p_gene_id))
		sys.stderr.write("Done\n")
		
	def run(self):
		"""
		03-01-05
			initial
		"""	
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		curs.execute("begin")	#because of cursor usage
		self.go_no2term_id = get_go_no2term_id(curs, self.schema, self.term_table)
		self.data_fetch(curs, self.p_gene_table, self.gene_p_table)
		self.p_gene_map(self.gene_no2p_gene, self.p_gene_id_map, curs,\
			self.distance_table, self.go_no2distance, self.go_no2term_id)
		if self.needcommit:
			self.submit(curs, self.gene_p_table, self.p_gene_id_map)
			curs.execute("end")


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)	
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:n:cru", ["help", "hostname=", \
			"dbname=", "schema=", "p_gene_table=", "gene_p_table=", \
			"commit", "report", "debug"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	p_gene_table = None
	gene_p_table = None
	commit = 0
	report = 0
	debug = 0
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
		elif opt in ("-t", "--p_gene_table"):
			p_gene_table = arg
		elif opt in ("-n", "--gene_p_table"):
			gene_p_table = arg
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-u", "--debug"):
			debug = 1
	if schema and p_gene_table and gene_p_table:
		instance = gene_p_map_redundancy(hostname, dbname, schema, p_gene_table, gene_p_table,\
			commit, report, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
