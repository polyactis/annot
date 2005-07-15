#!/usr/bin/env python
"""
Usage: HomoEdges.py -i INPUTFILE -o OUTPUTFILE

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database(IGNORE)
	-i ..., --inputfile=...	the inputfile, output of MpiGraphModeling.py
	-o ..., --outputfile=...	store the information of homo pairs
	-r, --report	report flag
	-u, --debug debug flag
	
Examples:
	HomoEdges.py -i gph_hs_GDS760 -o stat
	
Description:
	For each graph file(graph.bin output format), calculates the ratio of homogeneous pairs and outputs go_id2no_of_homo_pairs in outputfile.
"""

import sys, os, getopt, csv
from codense.common import db_connect
from sets import Set

def get_gene_id2go_id(curs, table='association', organism='Homo sapiens'):
	""""
	07-14-05
		get from graph.association, no parent nodes, disregard GO:0000004(unknown)
	"""
	sys.stderr.write("Getting gene_id2go_id...")
	gene_id2go_id = {}
	curs.execute("select gene_id, go_id from graph.%s where organism='%s' and go_id!='GO:0000004'"%(table, organism))
	rows = curs.fetchall()
	for row in rows:
		if row[0] not in gene_id2go_id:
			gene_id2go_id[row[0]] = Set()
		gene_id2go_id[row[0]].add( row[1])
	sys.stderr.write("Done\n")
	return gene_id2go_id

class HomoEdges:
	def __init__(self, hostname=None, dbname=None, schema=None, infname=None, \
		outfname=None, report=0, debug=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.infname = infname
		self.outfname = outfname
		self.report = int(report)		
		self.debug = int(debug)
	
	def find_homo_pairs(self, infname, gene_id2go_id):
		"""
		07-14-05
			return no_of_homo_pairs and go_id2no_of_homo_pairs
		"""
		sys.stderr.write("Finding homo pairs...")
		reader = csv.reader(open(infname,'r'), delimiter='\t')
		no_of_total_pairs = 0
		no_of_homo_pairs = 0
		go_id2no_of_homo_pairs = {}
		for row in reader:
			if row[0] == 'e':
				no_of_total_pairs += 1
				edge = row[1:3]
				common_go_id_set = self.return_common_go_id_set(edge, gene_id2go_id)
				if common_go_id_set:	#non-empty Set
					no_of_homo_pairs += 1
					for go_id in common_go_id_set:
						if go_id not in go_id2no_of_homo_pairs:
							go_id2no_of_homo_pairs[go_id] = 1
						else:
							go_id2no_of_homo_pairs[go_id] += 1
				if self.debug:
					print "no_of_homo_pairs: ",no_of_homo_pairs
					print "go_id2no_of_homo_pairs: ",go_id2no_of_homo_pairs
		del reader
		sys.stderr.write("Done.\n")
		return (no_of_total_pairs, no_of_homo_pairs, go_id2no_of_homo_pairs)

	def return_common_go_id_set(self, edge, gene_id2go_id):
		"""
		07-14-05
			find the common_go_id_set from two genes' associated go_id_set
		"""
		common_go_id_set = Set()
		go_id_set1 = gene_id2go_id.get(edge[0])
		go_id_set2 = gene_id2go_id.get(edge[1])
		if self.debug:
			print edge
			print edge[0], go_id_set1
			print edge[1], go_id_set2
			is_continue = raw_input("Continue?")
		if go_id_set1 and go_id_set2:
			common_go_id_set = go_id_set1 & go_id_set2
		return common_go_id_set
	
	def output_go_id2no_of_homo_pairs(self, outfname, no_of_homo_pairs, no_of_total_pairs, go_id2no_of_homo_pairs):
		"""
		07-14-05
			output go_id2no_of_homo_pairs

		"""
		sys.stderr.write("Outputting go_id2no_of_homo_pairs...")
		writer = csv.writer(open(outfname, 'w'), delimiter='\t')
		#writer.writerow(["No of homo pairs: %s/%s(%s)"%(no_of_homo_pairs, no_of_total_pairs, float(no_of_homo_pairs)/no_of_total_pairs)])
			#Don't uncomment the line above. Otherwise the output is hard to be processed further by other programs.
		for go_id, no_of_homo_pairs in go_id2no_of_homo_pairs.iteritems():
			writer.writerow([go_id, no_of_homo_pairs])
		sys.stderr.write("Done.\n")
		del writer
	
	def run(self):
		"""
		07-14-05
			--db_connect()
			--get_gene_id2go_id()
			--find_homo_pairs()
				--return_common_go_id_set()
			--output_go_id2no_of_homo_pairs()
		"""
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		gene_id2go_id = get_gene_id2go_id(curs)
		no_of_total_pairs, no_of_homo_pairs, go_id2no_of_homo_pairs = self.find_homo_pairs(self.infname, gene_id2go_id)
		print "No of homo pairs: %s/%s(%s)"%(no_of_homo_pairs, no_of_total_pairs, float(no_of_homo_pairs)/no_of_total_pairs)
		
		self.output_go_id2no_of_homo_pairs(self.outfname, no_of_homo_pairs, no_of_total_pairs, go_id2no_of_homo_pairs)


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)	
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:o:ru", ["help", "hostname=", \
			"dbname=", "schema=", "inputfile=", "outputfile=", "report", "debug"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = None
	infname = None
	outfname = None
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
		elif opt in ("-i", "--inputfile"):
			infname = arg
		elif opt in ("-o", "--outputfile"):
			outfname = arg
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-u", "--debug"):
			debug = 1
	if infname and outfname:
		instance = HomoEdges(hostname, dbname, schema, infname, outfname, report, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
