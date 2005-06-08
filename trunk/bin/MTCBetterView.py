#!/usr/bin/env python
"""
Usage: MTCBetterView.py -k SCHEMA -i INPUTFILE -o OUTPUTFILE [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-i ..., --infname=...	the input file(outputed by MpiTightClust.py)
	-o ..., --outfname=...	the output file.
	-b, --debug	enable debugging, no debug by default
	-h, --help              show this help

Examples:
	MTCBetterView.py -k mm_73 -i tightClust -o /tmp/datasetPartition.txt

Description:
	This program beautifies the output results of MpiTightClust.py.

"""

import sys, os, psycopg, getopt, csv, re
from sets import Set
from codense.common import db_connect


class MTCBetterView:
	"""
	06-08-05
		a class to beautify the output results of MpiTightClust.py.
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None,\
		infname=None, outfname=None, debug=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.infname = infname
		self.outfname = outfname
		self.debug = int(debug)
		
		self.p_go_no = re.compile(r'edge_data_(\d+)\.\d')
	
	def headerOutput(self, curs, outf):
		"""
		06-08-05
			output the header information
			return the number of datasets
		"""
		sys.stderr.write("Outputing header...")
		dataset_no_list = []
		dataset_id_list = []
		dataset_desc_list = []
		curs.execute("select d.*,dd.description from dataset_no2id d, graph.dataset_desc dd\
			where d.dataset_id=dd.dataset_id  order by dataset_no")
		rows = curs.fetchall()
		for row in rows:
			dataset_no_list.append(row[0])
			dataset_id_list.append(row[1])
			dataset_desc_list.append(row[2])
		dataset_no_list.append('dataset_no')
		dataset_id_list.append('dataset_id')
		dataset_desc_list.append('dataset_desc')
		
		outf.writerow(dataset_desc_list)
		outf.writerow(dataset_id_list)
		outf.writerow(dataset_no_list)
		sys.stderr.write("Done.\n")
		return len(dataset_no_list)-1
		
	def datasetClustOutput(self, curs, outf, row, no_of_datasets):
		"""
		06-08-05
			output the selected datasets, the funtion category from which it's from
		
			global variable used
				self.p_go_no
			
			--return_go_name()
		"""
		dataset_cluster = row[:-1]	#the last one is blank
		if len(dataset_cluster)==1:	#no dataset_no, only the function category
			return
		has_go_no = self.p_go_no.search(dataset_cluster[0])
		if has_go_no:
			go_no = has_go_no.groups()[0]
		else:
			sys.stderr.write("%s doesn't have go_no embeded. Abort.\n"%dataset_cluster[0])
			sys.exit(2)
		go_name = self.return_go_name(curs, go_no)
		
		dataset_no_set = Set(map(int, dataset_cluster[1:]))
			#integer form, to judge whether one dataset is included in the cluster.
		output_list = []	#this'll be outputed into the outfname
		for i in range(no_of_datasets):
			if i in dataset_no_set:
				output_list.append(i+1)	#Note: i+1(starting from 1)
			else:
				output_list.append('')
		output_list.append('%s(%s)'%(go_name, go_no))

		if self.debug:
			print output_list
			continue_here = raw_input("Continue?(Y/n)")
			if continue_here=='n':
				sys.exit(3)
		outf.writerow(output_list)
		
	def return_go_name(self, curs, go_no, go_table='go'):
		"""
		06-08-05
			return the name of the GO function category
		"""
		curs.execute("select name from go where go_no=%s"%go_no)
		rows = curs.fetchall()
		go_name=  rows[0][0]
		return go_name
	
	def run(self):
		"""
		06-08-05
			
			--db_connect()
			--headerOutput()
			--datasetClustOutput()
				--return_go_name()
		"""
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		outf = csv.writer(open(self.outfname, 'w'), delimiter='\t')
		no_of_datasets = self.headerOutput(curs, outf)
		
		reader = csv.reader(open(self.infname, 'r'), delimiter='\t')
		for row in reader:
			if self.debug:
				print row
			self.datasetClustOutput(curs, outf, row, no_of_datasets)
		
		del reader
		del outf

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema=", "infname=", "outfname=", \
		"debug"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:o:b", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	infname = None
	outfname = None
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
		elif opt in ("-i", "--infname"):
			infname = arg
		elif opt in ("-o", "--outfname"):
			outfname = arg
		elif opt in ("-b", "--debug"):
			debug = 1

	if schema and infname and outfname:
		instance = MTCBetterView(hostname, dbname, schema, infname, outfname, \
			debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
