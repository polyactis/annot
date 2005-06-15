#!/usr/bin/env python
"""
Usage: MTCBetterView.py -k SCHEMA -i INPUTFILE -o OUTPUTFILE [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-i ..., --infname=...	the input file(outputed by MpiTightClust.py)
	-o ..., --outfname=...	the output file.
	-y ..., --type=...	1(default, output reformat), 2(group dataset clusters)
	-s ..., --similar_score=...	percentage of sharing between two
		dataset_clusters, 0.8(default)
	-b, --debug	enable debugging, no debug by default
	-h, --help              show this help

Examples:
	MTCBetterView.py -k mm_73 -i tightClust -o /tmp/datasetPartition.txt

Description:
	This program beautifies the output results of MpiTightClust.py.

"""

import sys, os, psycopg, getopt, csv, re
from sets import Set
from codense.common import db_connect, index_plus_one
from graph.cc_from_edge_list import cc_from_edge_list
from CcFromBiclusteringOutput import CcFromBiclusteringOutput

class MTCBetterView:
	"""
	06-08-05
		a class to beautify the output results of MpiTightClust.py.
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None,\
		infname=None, outfname=None, type=1, similar_score=0.8, debug=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.infname = infname
		self.outfname = outfname
		self.type = int(type)
		self.similar_score = float(similar_score)
		self.debug = int(debug)
		
		self.p_go_no = re.compile(r'edge_data_(\d+)\.\d')
	
	def return_dataset_no2id_desc(self, curs):
		"""
		06-12-05
		"""
		dataset_no2id_desc = {}
		curs.execute("select d.*,dd.description from dataset_no2id d, graph.dataset_desc dd\
			where d.dataset_id=dd.dataset_id  order by dataset_no")
		rows = curs.fetchall()
		for row in rows:
			dataset_no2id_desc[row[0]] = list(row[1:])
		return dataset_no2id_desc
	
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
	
	def datasetClustOutput2(self, curs, outf, row, dataset_no2id_desc):
		"""
		06-12-05
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
		outf.writerow(['%s(function category)'%go_name, '%s(function number)'%go_no])
		
		dataset_no_list = list(map(index_plus_one, dataset_cluster[1:]))
		self._datasetClustOutput(outf, dataset_no_list, dataset_no2id_desc)
	
	def _datasetClustOutput(self, outf, dataset_no_list, dataset_no2id_desc):
		"""
		06-12-05
		"""
		dataset_no_list.sort()
		for i in range(len(dataset_no_list)):
			dataset_no = dataset_no_list[i]
			outf.writerow([dataset_no]+dataset_no2id_desc[dataset_no])
		outf.writerow([])	#one blank line
	
	def return_go_name(self, curs, go_no, go_table='go'):
		"""
		06-08-05
			return the name of the GO function category
		"""
		curs.execute("select name from go where go_no=%s"%go_no)
		rows = curs.fetchall()
		go_name=  rows[0][0]
		return go_name
	
	def id2dataset_cluster_setConstruct(self, infname):
		"""
		06-09-05
			contruct id2dataset_cluster_set
		"""
		sys.stderr.write("Constructing id2dataset_cluster_set...")
		id2dataset_cluster_set = {}
		reader = csv.reader(open(infname, 'r'), delimiter='\t')
		id = 0
		for row in reader:
			dataset_cluster = row[1:-1]
			if len(dataset_cluster)==0:	#empty
				continue
			dataset_cluster = map(index_plus_one, dataset_cluster)	#plus one here
			id2dataset_cluster_set[id] = Set(dataset_cluster)
			id+=1
		del reader
		sys.stderr.write("Done.\n")
		return id2dataset_cluster_set
	
	def dataset_clusterGraphConstruct(self, id2dataset_cluster_set, similar_score):
		"""
		06-09-05
		"""
		sys.stderr.write("Constructing dataset_cluster graph...")
		edge_list = []
		no_of_dataset_clusters = len(id2dataset_cluster_set)
		for i in range(no_of_dataset_clusters):
			dataset_cluster_set1 = id2dataset_cluster_set[i]
			edge_list.append([i,i])	#append the self-loop to avoid the singleton to be left out.
			for j in range(i+1, no_of_dataset_clusters):
				dataset_cluster_set2 = id2dataset_cluster_set[j]
				join_set = dataset_cluster_set1 & dataset_cluster_set2
				ratio1 = len(join_set)/float(len(dataset_cluster_set1))
				ratio2 = len(join_set)/float(len(dataset_cluster_set2))
				if ratio1>=similar_score or ratio2>=similar_score:
					if self.debug:
						l1 = list(dataset_cluster_set1)
						l2 = list(dataset_cluster_set2)
						l1.sort()
						l2.sort()
						print "edge valid with similar_score",ratio1,ratio2,i,j
						print l1
						print l2
						#continue_here = raw_input("Continue:(Y/n)?")
						#if continue_here=='n':
						#	sys.exit(3)
					edge_list.append([i,j])
		sys.stderr.write("Done.\n")
		return edge_list
	
	def returnBigDatasetClust(self, id2dataset_cluster_set, id_set):
		"""
		06-09-05
		"""
		big_dataset_cluster_set = Set()
		for id in id_set:
			big_dataset_cluster_set |= id2dataset_cluster_set[id]
		return big_dataset_cluster_set
	
	def run(self):
		"""
		06-08-05
		
		06-09-05
			add type 2: group dataset clusters
			
			--db_connect()
			--headerOutput()
			if self.type==1:
				--datasetClustOutput()
					--return_go_name()
			elif self.type==2:
				--id2dataset_cluster_setConstruct()
				--dataset_clusterGraphConstruct()
				--<cc_edge_list>
				--<CcFromBiclusteringOutput>
					--returnBigDatasetClust()
		"""
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		outf = csv.writer(open(self.outfname, 'w'), delimiter='\t')
		#no_of_datasets = self.headerOutput(curs, outf)
		dataset_no2id_desc = self.return_dataset_no2id_desc(curs)
		if self.type==1:
			reader = csv.reader(open(self.infname, 'r'), delimiter='\t')
			for row in reader:
				if self.debug:
					print row
				self.datasetClustOutput2(curs, outf, row, dataset_no2id_desc)
			del reader
		elif self.type==2:
			id2dataset_cluster_set = self.id2dataset_cluster_setConstruct(infname)	#06-09-05	mapping between an id and a dataset cluster set
			if self.debug:
				print "id2dataset_cluster_set is:"
				print id2dataset_cluster_set
			edge_list = self.dataset_clusterGraphConstruct(id2dataset_cluster_set, self.similar_score)
			if self.debug:
				print "The constructed graph has %s edges"%len(edge_list)
			cfe_instance= cc_from_edge_list()
			cfe_instance.run(edge_list)
			cfbo_instance = CcFromBiclusteringOutput()
			for cc_edge_list in cfe_instance.cc_list:
				id_set = cfbo_instance.vertex_set_from_cc_edge_list(cc_edge_list)
				if self.debug:
					print cc_edge_list
					print id_set
				big_dataset_cluster_set = self.returnBigDatasetClust(id2dataset_cluster_set, id_set)
				big_dataset_cluster = list(big_dataset_cluster_set)
				big_dataset_cluster.sort()
				self._datasetClustOutput(outf, big_dataset_cluster, dataset_no2id_desc)
				
		del outf

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema=", "infname=", "outfname=", \
		"type=", "similar_score=", "debug"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:o:y:s:b", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	infname = None
	outfname = None
	type = 1
	similar_score = 0.8
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
		elif opt in ("-y", "--type"):
			type = int(arg)
		elif opt in ("-s", "--similar_score"):
			similar_score = float(arg)
		elif opt in ("-b", "--debug"):
			debug = 1

	if schema and infname and outfname:
		instance = MTCBetterView(hostname, dbname, schema, infname, outfname, \
			type, similar_score, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
