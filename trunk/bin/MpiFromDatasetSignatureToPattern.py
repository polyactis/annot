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

edge2occurrence_vector = {}	#08-09-05	global structure to store the occurrence vectors of edges appearing in the patterns

def patternFormation(offset, parameter_list):
	"""
	08-07-05
		datasetSignatureFname is outputed by fim_closed
		intermediateFile is outputed by outputEdgeData()
	08-08-05
		store the occurrence_vector in edge2occurrence_vector
	08-24-05
		edge_sig_matrix replaces the intermediateFile
		
		(loop)
			--encodeOccurrence()
		(loop)
			--decodeOccurrenceBv()
			--outputCcFromEdgeList()
				--get_combined_vector()
				--codense2db_instance.parse_recurrence()
	"""
	offset = int(offset)
	datasetSignatureFname, intermediateFile, outputfile, offset_list, no_cc, \
		edge_sig_matrix, no_of_datasets, debug = parameter_list
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
			"""
			raw_input will stop the paralle program
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
	#iF = csv.reader(open(intermediateFile,'r'), delimiter='\t')	#08-24-05 intermediateFile defunct
	of = open(node_outputfile, 'w')
	codense2db_instance = codense2db()
	counter = 0
	for row in edge_sig_matrix:
		counter += 1
		edge = map(int, row[:2])
		occurrenceBinaryForm = row[2]	#08-24-05	already encoded when edge_sig_matrix is filled in
		occurrence_vector = decodeOccurrenceToBv(occurrenceBinaryForm, no_of_datasets)
		signatureToBeDeleted = []
		for signature in signature2pattern:
			if (occurrenceBinaryForm&signature)==signature:
				signature2pattern[signature].append(edge)
				if debug:
					sys.stderr.write("the occurrence_vector of edge %s is %s\n"%(repr(edge), repr(occurrence_vector)))
					sys.stderr.write("occurrence_vector's binary form is %s, signature is %s\n"%(occurrenceBinaryForm, signature))
				if len(signature2pattern[signature]) == signature2pattern[signature][0]+1:	#the 1st entry is not edge
					signatureToBeDeleted.append(signature)
					if debug:
						sys.stderr.write("signature %s to be deleted, its pattern is %s\n"%(signature, repr(signature2pattern[signature])))
				edge_tuple = tuple(edge)
				if edge_tuple not in edge2occurrence_vector:
					edge2occurrence_vector[edge_tuple] = [1]
					edge2occurrence_vector[edge_tuple].append(occurrence_vector)
				else:
					edge2occurrence_vector[edge_tuple][0] += 1
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
	#del iF
	of.close()
	return node_outputfile

def outputCcFromEdgeList(of, signature, edge_list, codense2db_instance, no_cc):
	"""
	08-07-05
	08-09-05
		calculate recurrence array for codense2db.py
	"""
	if no_cc:
		vertex_set = codense2db_instance.vertex_set_from_cc_edge_list(edge_list)
		vertex_set.sort()
		combined_vector = get_combined_vector(edge_list)
		recurrence_array = codense2db_instance.parse_recurrence(combined_vector)
		of.write('%s\t%s\t%s\n'%(repr(vertex_set), repr(edge_list), repr(recurrence_array) ) )
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
			combined_vector = get_combined_vector(cc_edge_list)
			recurrence_array = codense2db_instance.parse_recurrence(combined_vector)
			of.write('%s\t%s\t%s\n'%(repr(vertex_set), repr(cc_edge_list), repr(recurrence_array) ) )

def get_combined_vector(edge_list):
	"""
	08-09-05
		get combined_vector from global structure: edge2occurrence_vector
		Reduce the counter of the edge, if the counter == 0, delete it.
	"""
	combined_vector = []
	for edge in edge_list:
		edge_tuple = tuple(edge)
		combined_vector.append(edge2occurrence_vector[edge_tuple][1])
		edge2occurrence_vector[edge_tuple][0] -= 1
		if edge2occurrence_vector[edge_tuple][0] == 0:
			del edge2occurrence_vector[edge_tuple]
	return combined_vector

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

def encodeOccurrenceBv(ls):
	"""
	08-06-05
		encode an occurrence vector to a binary number
	"""
	binary_number = 0
	for i in range(len(ls)):
		if ls[i] == 1:
			binary_number += int(math.pow(2, i))
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

def decodeOccurrenceToBv(signature, no_of_datasets):
	"""
	08-24-05
		similar to decodeOccurrence(), but it's binary vector of length no_of_datasets
	"""
	occurrence_vector = [0]*no_of_datasets
	binary_list = []
	if signature <=0:
		return [0]
	while signature!=1:
		signature, m = divmod(signature, 2)
		binary_list.append(m)
		if m==1:
			occurrence_vector[len(binary_list)-1] = 1
	binary_list.append(1)	#the last 1
	occurrence_vector[len(binary_list)-1] = 1
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
	
	def get_no_of_edges_and_no_of_datasets(self, curs, edge_table='edge_cor_vector'):
		"""
		08-24-05
			get the no_of_edges and no_of_datasets
				no_of_edges got from edge_cor_vector_edge_id_seq
				no_of_datasets got from the first sig_vector
		"""
		sys.stderr.write("Getting no_of_edges and no_of_datasets...")
		array_to_return = Numeric.zeros(2, Numeric.Int)
		curs.execute("select last_value from %s_edge_id_seq"%edge_table)
		rows = curs.fetchall()
		array_to_return[0] = rows[0][0]
		if array_to_return[0]<=1:
			sys.stderr.write("edge_cor_vector_edge_id_seq says there's <=1 edge. something wrong, exit.\n")
			sys.exit(2)
		curs.execute("select sig_vector from %s limit 1"%edge_table)
		rows = curs.fetchall()
		for row in rows:
			sig_vector = row[0][1:-1].split(',')
			no_of_datasets = len(sig_vector)
		array_to_return[1] = no_of_datasets
		sys.stderr.write("Done.\n")
		return array_to_return
	
	def outputEdgeData(self, curs, outputfile, min_sup, max_sup, edge_table='edge_cor_vector'):
		"""
		08-06-05
			Output edge data into an intermediate file to let other nodes to read.
			(edge pair + occurrence vector)
		08-08-05
			output sig_vector directly, not occurrence_vector
		"""
		sys.stderr.write("Getting edge matrix for all functions...\n")
		writer = csv.writer(open(outputfile, 'w'), delimiter='\t')
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
					"""
					new_row = edge
					for i in range(len(sig_vector)):
						if sig_vector[i]==1:
							new_row.append(i+1)
					"""
					writer.writerow(edge+sig_vector)
				counter +=1
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, counter))
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		sys.stderr.write("Done\n")
		del writer
	
	def fillEdgeSigMatrix(self, curs, edge_sig_matrix, min_sup, max_sup, edge_table='edge_cor_vector'):
		"""
		08-24-05
			similar to outputEdgeData() but read it into edge_sig_matrix, instead of output
		"""
		sys.stderr.write("Getting edge matrix for all functions...\n")
		curs.execute("DECLARE crs CURSOR FOR select edge_name,sig_vector \
			from %s"%(edge_table))
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		counter = 0
		pointer = 0
		while rows:
			for row in rows:
				edge = row[0][1:-1].split(',')
				edge = map(int, edge)
				sig_vector = row[1][1:-1].split(',')
				sig_vector = map(int, sig_vector)
				
				if sum(sig_vector)>=min_sup and sum(sig_vector)<=max_sup:
					"""
					new_row = edge
					for i in range(len(sig_vector)):
						if sig_vector[i]==1:
							new_row.append(i+1)
					"""
					edge_sig_matrix[pointer] = edge[0], edge[1], encodeOccurrenceBv(sig_vector)
					pointer += 1
				else:
					if self.debug:
						sys.stderr.write("edge %s with support %s, not included.\n"%(repr(edge), sum(sig_vector)))
				counter +=1
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, counter))
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		sys.stderr.write("Done\n")
		return pointer	#return the number of edges satisfying the min_sup and max_sup
		
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
		08-24-05
			read all edge data into matrix
		"""
		communicator = MPI.world.duplicate()
		
		intermediateFile = '%s.int'%self.outputfile
		offset_list = Numeric.zeros((communicator.size), Numeric.Int)	#not communicator.size-1, 
		two_number_array = Numeric.zeros(2, Numeric.Int)	#08-24-05, no_of_edges and no_of_datasets
		if communicator.rank == 0:
			sys.stderr.write("this is node %s\n"%communicator.rank)
			(conn, curs) = db_connect(self.hostname, self.dbname, self.schema)
			offset_list[:] = self.createOffsetList(self.inputfile, communicator.size-1)
			if self.debug:
				sys.stderr.write("offset_list: %s\n"%repr(offset_list))
			two_number_array = self.get_no_of_edges_and_no_of_datasets(curs)
			if self.debug:
				sys.stderr.write("two_number_array is %s\n"%repr(two_number_array))
			del conn, curs
		
		communicator.broadcast(offset_list, 0)	#share the offset_list
		mpi_synchronize(communicator)
		
		communicator.broadcast(two_number_array, 0)	#broadcast no_of_edges and no_of_datasets
		mpi_synchronize(communicator)
		no_of_edges, no_of_datasets = two_number_array
		#read all edges and its encoded occurrence into edge_sig_matrix and broadcast it
		edge_sig_matrix = Numeric.zeros((no_of_edges, 3), Numeric.Int)
		real_no_of_edges = Numeric.zeros(1,Numeric.Int)
		if communicator.rank == 0:
			sys.stderr.write("this is node %s\n"%communicator.rank)
			(conn, curs) = db_connect(self.hostname, self.dbname, self.schema)
			real_no_of_edges[0] = self.fillEdgeSigMatrix(curs, edge_sig_matrix, self.min_sup, self.max_sup)
			if self.debug:
				sys.stderr.write("real_no_of_edges is %s\n"%repr(real_no_of_edges))
			del conn, curs
		communicator.broadcast(edge_sig_matrix, 0)
		mpi_synchronize(communicator)
		
		#broadcast real_no_of_edges
		communicator.broadcast(real_no_of_edges, 0)
		mpi_synchronize(communicator)
		#resize the edge_sig_matrix, throw away additional zeros
		edge_sig_matrix = Numeric.resize(edge_sig_matrix, (real_no_of_edges[0], 3))
		
		job_list = range(communicator.size-1)	#corresponding to the indices in the offset_list
		parameter_list =[self.inputfile, intermediateFile, self.outputfile, offset_list, self.no_cc, edge_sig_matrix,  no_of_datasets, self.debug]
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
