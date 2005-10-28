#!/usr/bin/env mpipython
"""
Usage: MpiRecurrenceFilter.py -k SCHEMA -i INPUTFILE -o OUTPUTFILE [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-i ...,	--inputfile=...	the input file, dataset signature output by fim
	-o ..., --outputfile=...	the output file
	-m ..., --min_sup=...	minimum support of an edge in the database, 0(default)
	-x ..., --max_sup=...	maximum support of an edge, 200(default)
	-s ..., --min_size=...	minimum size of a cluster(5 default)
	-a ...,	alpha value(0.05, default)
	-v ...,	message size(10,000,000, default)
	-b,	debug version.
	-r,	enable report flag
	-h, --help	Display the usage infomation.
	
Examples:
	mpirun -np 20 -machinefile ~/hostfile /usr/bin/mpipython ~/script/annot/bin/MpiRecurrenceFilter.py
	-k rn_fim_33 -i ~/bin/hhu_clustering/data/output/netmine/rn_fim_33m3x100
	-o ~/bin/hhu_clustering/data/output/netmine/rn_fim_33m3x100.filter -m 3
	
Description:
	Program to filter clusters based on recurrence and no_of_edges.
	10-28-05 p-value log(ln)ged.
"""

import sys, os, getopt, csv, math, cPickle
sys.path += [os.path.expanduser('~/script/annot/bin')]
from Scientific import MPI
from codense.common import mpi_synchronize, db_connect, output_node, get_edge2occurrence,\
	computing_node, input_node
from sets import Set
from rpy import r

class MpiRecurrenceFilter:
	def __init__(self,hostname='zhoudb', dbname='graphdb', schema=None, inputfile=None,\
		outputfile=None, min_sup=0, max_sup=200,  min_size=5, alpha=0.05, \
		message_size=10000000, debug=0, report=0):
		"""
		10-22-05
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.inputfile = inputfile
		self.outputfile = outputfile
		self.min_sup = int(min_sup)
		self.max_sup = int(max_sup)
		self.min_size = int(min_size)
		self.alpha = float(alpha)
		self.message_size = int(message_size)
		self.debug = int(debug)
		self.report = int(report)
	
	def input_handler(self, parameter_list, message_size, report=0):
		"""
		10-22-05
		"""
		if report:
			sys.stderr.write("Fetching stuff...\n")
		inf = parameter_list[0]
		block = []
		string_length = 0
		for row in inf:
			block.append(row)
			string_length += len(repr(row))	#the length to control MPI message size
			if string_length>=message_size:
				break
		if report:
			sys.stderr.write("Fetching done.\n")
		return block
	
	def node_fire(self, communicator, data, parameter_list):
		"""
		10-22-05
		"""
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		min_size, alpha, edge2occurrrence, no_of_datasets = parameter_list
		result = []
		for row in data:
			vertex_set = row[0][1:-1].split(',')
			if len(vertex_set)<min_size:	#too small
				continue
			recurrence_array = row[2][1:-1].split(',')
			recurrence_array = map(float, recurrence_array)
			recurrence = sum([int(occurrence==1) for occurrence in recurrence_array])
			prob = 1.0
			edge_set = row[1][2:-2].split('], [')
			for edge in edge_set:
				edge = edge.split(',')
				edge = map(int, edge)
				edge.sort()
				prob *= edge2occurrrence[tuple(edge)]/float(no_of_datasets)
			p_value = r.pbinom(recurrence-1, no_of_datasets, prob, lower_tail=r.FALSE,log_p=r.TRUE)
			if p_value<=alpha:
				result.append(row)
		sys.stderr.write("Node no.%s done with %s/%s clusters left.\n"%(node_rank, len(result), len(data)))
		return result
	
	def output_handler(self, communicator, parameter_list, data):
		"""
		10-22-05
		"""
		writer = parameter_list[0]
		data = cPickle.loads(data)
		for row in data:
			writer.writerow(row)
		
	
	def run(self):
		"""
		10-22-05
			
		"""
		communicator = MPI.world.duplicate()
		node_rank = communicator.rank
		free_computing_nodes = range(1,communicator.size-1)
		print "this is node",node_rank
		if node_rank == 0:
			(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
			edge2occurrrence, no_of_datasets = get_edge2occurrence(curs, self.min_sup, self.max_sup)
			edge2occurrrence_pickle = cPickle.dumps((edge2occurrrence, no_of_datasets), -1)
			for node in free_computing_nodes:	#send it to the computing_node
				communicator.send(edge2occurrrence_pickle, node, 0)
			del conn, curs
		elif node_rank in free_computing_nodes:	#exclude the last node
			data, source, tag = communicator.receiveString(0, 0)
			edge2occurrrence, no_of_datasets = cPickle.loads(data)
		
		mpi_synchronize(communicator)
		if node_rank == 0:
			inf = csv.reader(open(self.inputfile,'r'), delimiter='\t')
			parameter_list = [inf]
			input_node(communicator, parameter_list, free_computing_nodes, self.message_size, self.report, input_handler=self.input_handler)
			del inf
		elif node_rank in free_computing_nodes:
			parameter_list = [self.min_size, self.alpha, edge2occurrrence, no_of_datasets]
			computing_node(communicator, parameter_list, self.node_fire, report=self.report)
		elif node_rank == communicator.size-1:
			writer = csv.writer(open(self.outputfile, 'w'), delimiter='\t')
			parameter_list = [writer]
			output_node(communicator, free_computing_nodes, parameter_list, self.output_handler, self.report)
			del writer


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:o:m:x:s:a:v:br", ["help", "hostname=", \
			"dbname=", "schema=", "inputfile=", "outputfile=", "min_sup", "max_sup=", "min_size=", \
			"debug", "report"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	inputfile = None
	outputfile = None
	min_sup = 0
	max_sup = 200
	min_size = 5
	alpha = 0.05
	message_size = 10000000
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
		elif opt in ("-i", "--inputfile"):
			inputfile = arg
		elif opt in ("-o", "--outputfile"):
			outputfile = arg
		elif opt in ("-m", "--min_sup"):
			min_sup = int(arg)
		elif opt in ("-x", "--max_sup"):
			max_sup = int(arg)
		elif opt in ("-s", "--min_size"):
			min_size = int(arg)
		elif opt in ("-a"):
			alpha = float(arg)
		elif opt in ("-v"):
			message_size = int(arg)
		elif opt in ("-b"):
			debug = 1
		elif opt in ("-r"):
			report = 1
	if schema and inputfile and outputfile:
		instance = MpiRecurrenceFilter(hostname, dbname, schema, inputfile, outputfile, \
			min_sup, max_sup, min_size, alpha, message_size, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
