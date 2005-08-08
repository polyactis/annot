#!/usr/bin/env mpipython
"""
Usage: MpiFromDatasetSignatureToPattern.py -k SCHEMA -i INPUTFILE -o OUTPUTFILE [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-i ...,	--inputfile=...	the input file, dataset signature output by fim
	-o ..., --outputfile=...	the output file
	-m ..., --min_sup=...	minimum support of an edge in the database, 0(default)
	-x ..., --max_sup=...	maximum support of an edge, 200(default)
	-n, --no_cc	no connected components extraction
	-b, --debug	debug version.
	-r, --report	enable report flag
	-h, --help	Display the usage infomation.
	
Examples:
	mpirun.mpich -np 30 -machinefile ~/hostfile /usr/bin/mpipython ~/script/annot/bin/MpiFromDatasetSignatureToPattern.py
	-k mm_fim_97 -i /scratch/00/yuhuang/tmp/mm_fim_97.s4.l5.output -o /scratch/00/yuhuang/tmp/mm_fim_97.o
	
Description:
	Program to handle the fim_closed output over the edge-transaction mining
	Restore the patterns(vertex_set, edge_set) based on the dataset signature.
	
"""

import sys, os, getopt, csv, math, Numeric
sys.path += [os.path.expanduser('~/script/annot/bin')]
from Scientific import MPI
from graph.cc_from_edge_list import cc_from_edge_list
from codense.common import system_call, mpi_schedule_jobs, mpi_synchronize, db_connect
from netmine_wrapper import netmine_wrapper
from codense.codense2db import codense2db
from sets import Set

def patternFormation(offset, parameter_list):
	"""
	08-07-05
		datasetSignatureFname is outputed by fim_closed
		intermediateFile is outputed by outputEdgeData()
	"""
	offset = int(offset)
	datasetSignatureFname, intermediateFile, outputfile, offset_list, no_cc, debug = parameter_list
	node_outputfile = '%s.%s'%(outputfile, offset)
	#initialize the signatures in signature2pattern
	signature2pattern = {}
	dSF = csv.reader(open(datasetSignatureFname, 'r'), delimiter=' ')
	lower_bound = offset_list[offset]+1
	upper_bound = offset_list[offset+1]
	counter = 0
	for row in dSF:
		counter += 1
		if counter >= lower_bound and counter<= upper_bound:
			sig_row = map(int, row[:-1])
			frequency = int(row[-1][1:-1])
			occurrenceBinaryForm = encodeOccurrence(sig_row)
			if debug:
				sys.stderr.write("dataset signature is %s with frequency %s\n"%(repr(sig_row), frequency))
				sys.stderr.write("its binary form is %s\n"%occurrenceBinaryForm)
			"""raw_input will stop the paralle program
			is_continue = raw_input("Continue?(Y/n)")
			if is_continue=='n':
				sys.exit(2)
			"""
			if occurrenceBinaryForm not in signature2pattern:
				signature2pattern[occurrenceBinaryForm] = [frequency]
			else:
				sys.stderr.write("Something wrong, %s already appears in signature2pattern\n"%occurrenceBinaryForm)
				sys.exit(2)
	del dSF
	#fill patterns in signature2pattern
	iF = csv.reader(open(intermediateFile,'r'), delimiter='\t')
	of = open(node_outputfile, 'w')
	codense2db_instance = codense2db()
	counter = 0
	for row in iF:
		counter += 1
		edge = map(int, row[:2])
		occurrence_vector = map(int, row[2:])
		occurrenceBinaryForm = encodeOccurrence(occurrence_vector)
		signatureToBeDeleted = []
		for signature in signature2pattern:
			if (occurrenceBinaryForm&signature)==signature:
				signature2pattern[signature].append(edge)
				if debug:
					sys.stderr.write("the occurrence_vector of edge %s is %s\n"%(repr(edge), repr(occurrence_vector)))
					sys.stderr.write("occurrence_vector's binary form is %s, signature is %s\n"%(occurrenceBinaryForm, signature))
				if len(signature2pattern[signature]) == signature2pattern[signature][0]+1:
					signatureToBeDeleted.append(signature)
					if debug:
						sys.stderr.write("signature %s to be deleted, its pattern is %s\n"%(signature, repr(signature2pattern[signature])))
				"""raw_input will stop the paralle program
				is_continue = raw_input("Continue?(Y/n)")
				if is_continue=='n':
					sys.exit(2)
				"""
		for signature in signatureToBeDeleted:
			outputCcFromEdgeList(of, signature, signature2pattern[signature][1:], codense2db_instance, no_cc)
			del signature2pattern[signature]
	if len(signature2pattern)>1:
		sys.stderr.write('Weird %s signatures are still available\n'%len(signature2pattern))
		if debug:
			sys.stderr.write('%s\n'%repr(signature2pattern))
	del iF
	of.close()
	return node_outputfile

def outputCcFromEdgeList(of, signature, edge_list, codense2db_instance, no_cc):
	"""
	08-07-05
	"""
	if no_cc:
		vertex_set = codense2db_instance.vertex_set_from_cc_edge_list(edge_list)
		vertex_set.sort()
		of.write('%s\t%s\n'%(repr(vertex_set), repr(edge_list) ) )
	else:
		cf_instance = cc_from_edge_list()
		cf_instance.run(edge_list)
		cc_list = cf_instance.cc_list
		for cc_edge_list in cc_list:
			vertex_set = codense2db_instance.vertex_set_from_cc_edge_list(cc_edge_list)
			vertex_set.sort()
			cc_edge_list = map(list, cc_edge_list)	#change the tuple type to list
			for i in range(len(cc_edge_list)):
				cc_edge_list[i].sort()	#sort it
			of.write('%s\t%s\n'%(repr(vertex_set), repr(cc_edge_list) ) )
	
def encodeOccurrence(ls):
	"""
	08-06-05
		encode an occurrence vector to a binary number
	"""
	binary_number = 0
	for digit in ls:
		binary_number += int(math.pow(2, digit-1))	#int() because math.pow() returns float 
			#IMPORTANT: later makes the binary_number approximate cause it's too large
	return binary_number

def decodeOccurrence(signature):
	"""
	08-07-05
		decode a decimal number into a binary_list, and get the occurrence_vector
	"""
	occurrence_vector = []
	binary_list = []
	if signature <=0:
		return [0]
	while signature!=1:
		signature, m = divmod(signature, 2)
		binary_list.append(m)
		if m==1:
			occurrence_vector.append(len(binary_list))
	binary_list.append(1)	#the last 1
	occurrence_vector.append(len(binary_list))
	return occurrence_vector

class MpiFromDatasetSignatureToPattern:
	def __init__(self,hostname='zhoudb', dbname='graphdb', schema=None, inputfile=None,\
		outputfile=None, min_sup=0, max_sup=200, no_cc=0, debug=0, report=0):
		"""
		08-07-05
			Program to handle the fim_closed output over the edge-transaction mining
			Restore the patterns(vertex_set, edge_set) based on the dataset signature.
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.inputfile = inputfile
		self.outputfile = outputfile
		self.min_sup = int(min_sup)
		self.max_sup = int(max_sup)
		self.no_cc = int(no_cc)
		self.debug = int(debug)
		self.report = int(report)
	
	def outputEdgeData(self, curs, outputfile, min_sup, max_sup, edge_table='edge_cor_vector'):
		"""
		08-06-05
			Output edge data into an intermediate file to let other nodes to read.
			(edge pair + occurrence vector)
			
		"""
		writer = csv.writer(open(outputfile, 'w'), delimiter='\t')
		sys.stderr.write("Getting edge matrix for all functions...\n")
		curs.execute("DECLARE crs CURSOR FOR select edge_name,sig_vector \
			from %s"%(edge_table))
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		counter = 0
		while rows:
			for row in rows:
				edge = row[0][1:-1].split(',')
				edge = map(int, edge)
				sig_vector = row[1][1:-1].split(',')
				sig_vector = map(int, sig_vector)
				
				if sum(sig_vector)>=min_sup and sum(sig_vector)<=max_sup:
					new_row = edge
					for i in range(len(sig_vector)):
						if sig_vector[i]==1:
							new_row.append(i+1)
					writer.writerow(new_row)
				counter +=1
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, counter))
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		sys.stderr.write("Done\n")
		del writer

	def createOffsetList(self, inputfile, no_of_nodes):
		"""
		08-06-05
			create an offset list to separate the inputfile
				the offset_list has no_of_nodes+1 entries, last entry is no_of_lines
				first entry is 0.
		"""
		inf = open(inputfile, 'r')
		no_of_lines = 0
		for line in inf:
			no_of_lines += 1
		del inf
		offset_list = []
		step = no_of_lines/(no_of_nodes-1)	#this is integer division
		for i in range(no_of_nodes):
			offset_list.append(step*i)
		offset_list.append(no_of_lines)
		return offset_list
	
	def uniqueSort(self, inputfile, outputfile):
		"""
		08-07-05
			sort the file (unique it simutaneouly)
		"""
		sys.stderr.write("Starting to sort and unique %s..."%inputfile)
		commandline = 'sort -u %s > %s'%(inputfile, outputfile)
		exit_code = system_call(commandline)
		#os.remove(inputfile)
		sys.stderr.write("Done.\n")
	
	def run(self):
		"""
		08-06-05
		"""
		communicator = MPI.world.duplicate()
		
		intermediateFile = '%s.int'%self.outputfile
		offset_list = Numeric.zeros((communicator.size), Numeric.Int)	#not communicator.size-1, 
		if communicator.rank == 0:
			sys.stderr.write("this is node %s\n"%communicator.rank)
			(conn, curs) = db_connect(self.hostname, self.dbname, self.schema)
			self.outputEdgeData(curs, intermediateFile, self.min_sup, self.max_sup)
			offset_list[:] = self.createOffsetList(self.inputfile, communicator.size-1)
			if self.debug:
				sys.stderr.write("offset_list: %s"%repr(offset_list))
		
		communicator.broadcast(offset_list, 0)	#share the offset_list
		
		mpi_synchronize(communicator)
		job_list = range(communicator.size-1)	#corresponding to the indices in the offset_list
		parameter_list =[self.inputfile, intermediateFile, self.outputfile, offset_list, self.no_cc, self.debug]
		if self.debug:
			sys.stderr.write("The common parameter_list is %s.\n"%repr(parameter_list))
		of_name_list = mpi_schedule_jobs(communicator, job_list, patternFormation, parameter_list, self.debug)
		
		mpi_synchronize(communicator)
		
		#collecting
		if communicator.rank==0:
			intermediateFile = '%s.unsorted'%self.outputfile	#intermediateFile to store concatenated results
			netmine_wrapper_instance = netmine_wrapper()
			netmine_wrapper_instance.collect_and_merge_output(of_name_list, intermediateFile)
			self.uniqueSort(intermediateFile, self.outputfile)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:o:y:m:x:nbr", ["help", "hostname=", \
			"dbname=", "schema=", "inputfile=", "outputfile=", "min_sup", "max_sup=", \
			"no_cc", "debug", "report"])
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
	no_cc = 0
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
		elif opt in ("-n", "--no_cc"):
			no_cc = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	if schema and inputfile and outputfile:
		instance = MpiFromDatasetSignatureToPattern(hostname, dbname, schema, \
			inputfile, outputfile, min_sup, max_sup, no_cc, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
