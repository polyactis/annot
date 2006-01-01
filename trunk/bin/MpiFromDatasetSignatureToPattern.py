#!/usr/bin/env mpipython
"""
Usage: MpiFromDatasetSignatureToPattern.py -i INPUTFILE -s SIG_VECTOR_FILE -o OUTPUTFILE [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default, USELESS)
	-d ..., --dbname=...	the database name, graphdb(default, USELESS)
	-k ..., --schema=...	which schema in the database(USELESS)
	-i ...,	--inputfile=...	the input file, dataset signature output by fim
	-s ...,	sig_vector file
	-o ..., --outputfile=...	the output file
	-m ..., --min_sup=...	minimum support of an edge in the database, 0(default)
	-x ..., --max_sup=...	maximum support of an edge, 200(default)
	-q ...,	edge_sig_vector_queue size for each computing node, 3,000,000(default)
	-n, --no_cc	no connected components extraction
	-b, --debug	debug version.
	-r, --report	enable report flag
	-h, --help	Display the usage infomation.
	
Examples:
	mpirun.mpich -np 30 -machinefile ~/hostfile /usr/bin/mpipython ~/script/annot/bin/MpiFromDatasetSignatureToPattern.py
	-s mm_fim_97_4.sig_vector -i mm_fim_97.s4.l5.output -o mm_fim_97.o
	
Description:
	Program to handle the fim_closed output over the edge-transaction mining
	Restore the patterns(vertex_set, edge_set) based on the dataset signature.
	
"""

import sys, os, getopt, csv, math, Numeric, cPickle
sys.path += [os.path.expanduser('~/script/annot/bin')]
from Scientific import MPI
from graph.cc_from_edge_list import cc_from_edge_list
from codense.common import system_call, mpi_synchronize, db_connect
from netmine_wrapper import netmine_wrapper
from codense.codense2db import codense2db
from sets import Set
if sys.version_info[:2] < (2, 3):       #python2.2 or lower needs some extra
        from python2_3 import *

from threading import *
from Queue import Queue
class PatternFormThread(Thread):
	"""
	12-31-05
		a thread on each computing node to do patternFormation()
	"""
	def __init__(self, rank, parameter_list):
		Thread.__init__(self)
		self.rank = rank
		self.parameter_list = parameter_list
		self.node_outputfile = None
	
	def run(self):
		sys.stderr.write("Thread from node %s starts ...\n"%(self.rank))
		self.node_outputfile = patternFormation(self.rank-1, self.parameter_list)
		sys.stderr.write("Thread from node %s ends.\n"%(self.rank))
	
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
	datasetSignatureFname, outputfile, offset_list, no_cc, \
		edge_sig_vector_queue, no_of_datasets, debug = parameter_list
	#initialize the signatures in signature2pattern, by reading from datasetSignatureFname
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
			signature2pattern[occurrenceBinaryForm] = [frequency]
		elif counter>upper_bound:	#12-30-05	skip the other rows
			break
	del dSF
	
	node_outputfile = '%s.%s'%(outputfile, offset)
	of = open(node_outputfile, 'w')
	codense2db_instance = codense2db()
	counter = 0
	edge_occurrenceBinaryForm_row = edge_sig_vector_queue.get()
	while edge_occurrenceBinaryForm_row!= -1:
		counter += 1
		edge = edge_occurrenceBinaryForm_row[:2]
		occurrenceBinaryForm = edge_occurrenceBinaryForm_row[2]	#08-24-05	already encoded when edge_sig_matrix is filled in
		#occurrence_vector = decodeOccurrenceToBv(occurrenceBinaryForm, no_of_datasets)
		signatureToBeDeleted = []
		for signature in signature2pattern:
			frequency = signature2pattern[signature][0]
			if (occurrenceBinaryForm&signature)==signature:
				signature2pattern[signature].append(edge)
				if debug:
					sys.stderr.write("the occurrence_vector of edge %s is %s\n"%(repr(edge), \
						repr(decodeOccurrenceToBv(occurrenceBinaryForm, no_of_datasets))))
					sys.stderr.write("occurrence_vector's binary form is %s, signature is %s\n"%(occurrenceBinaryForm, signature))
				if len(signature2pattern[signature]) == frequency+1:	#the 1st entry is frequency
					signatureToBeDeleted.append(signature)
					if debug:
						sys.stderr.write("signature %s to be deleted, its pattern is %s\n"%(signature, repr(signature2pattern[signature])))
				"""
				edge_tuple = tuple(edge)
				if edge_tuple not in edge2occurrence_vector:
					edge2occurrence_vector[edge_tuple] = [1]
					edge2occurrence_vector[edge_tuple].append(occurrence_vector)
				else:
					edge2occurrence_vector[edge_tuple][0] += 1
				"""
		for signature in signatureToBeDeleted:
			edge_list = signature2pattern[signature][1:]
			outputCcFromEdgeList(of, signature, edge_list, codense2db_instance, no_cc)
			del signature2pattern[signature]
		edge_occurrenceBinaryForm_row = edge_sig_vector_queue.get()
	if len(signature2pattern)>1:
		sys.stderr.write('Weird %s signatures are still available\n'%len(signature2pattern))
		if debug:
			sys.stderr.write('%s\n'%repr(signature2pattern))
	of.close()
	return node_outputfile

def outputCcFromEdgeList(of, signature, edge_list, codense2db_instance, no_cc):
	"""
	08-07-05
	08-09-05
		calculate recurrence array for codense2db.py
	12-31-05
		remove several time-consuming steps, but vertex_set and cc_edge_list are not sorted anymore
		no recurrence_array
		cc_edge_list is tuple-list
	"""
	if no_cc:
		vertex_set = codense2db_instance.vertex_set_from_cc_edge_list(edge_list)
		vertex_set.sort()
		#combined_vector = get_combined_vector(edge_list)
		#recurrence_array = codense2db_instance.parse_recurrence(combined_vector)
		of.write('%s\t%s\n'%(repr(vertex_set), repr(edge_list)) )
	else:
		cf_instance = cc_from_edge_list()
		cf_instance.run(edge_list)
		cc_list = cf_instance.cc_list
		for cc_edge_list in cc_list:
			vertex_set = codense2db_instance.vertex_set_from_cc_edge_list(cc_edge_list)
			vertex_set.sort()
			"""
			#12-31-05	each edge in cc_edge_list is already sorted
			cc_edge_list = map(list, cc_edge_list)  #change the tuple type to list
			for i in range(len(cc_edge_list)):
				cc_edge_list[i].sort()  #sort it
			"""
			cc_edge_list.sort()
			#combined_vector = get_combined_vector(cc_edge_list)
			#recurrence_array = codense2db_instance.parse_recurrence(combined_vector)
			of.write('%s\t%s\n'%(repr(vertex_set), repr(cc_edge_list) ) )

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
	12-30-05
		python2.2, int is can't handle long integer. replace it with long.
	"""
	binary_number = 0
	for digit in ls:
		binary_number += long(math.pow(2, digit-1))	#int() because math.pow() returns float 
			#IMPORTANT: later makes the binary_number approximate cause it's too large
	return binary_number

def encodeOccurrenceBv(ls):
	"""
	08-06-05
		encode an occurrence vector to a binary number
	12-30-05
		python2.2, int is can't handle long integer. replace it with long.
	"""
	binary_number = 0
	for i in range(len(ls)):
		if ls[i] == 1:
			binary_number += long(math.pow(2, i))
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
		sig_vector_fname=None, outputfile=None, min_sup=0, max_sup=200, queue_size=3000000,\
		no_cc=0, debug=0, report=0):
		"""
		08-07-05
			Program to handle the fim_closed output over the edge-transaction mining
			Restore the patterns(vertex_set, edge_set) based on the dataset signature.
		12-31-05
			add queue_size
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.inputfile = inputfile
		self.sig_vector_fname = sig_vector_fname
		self.outputfile = outputfile
		self.min_sup = int(min_sup)
		self.max_sup = int(max_sup)
		self.queue_size = int(queue_size)
		self.no_cc = int(no_cc)
		self.debug = int(debug)
		self.report = int(report)
	
	def get_no_of_datasets(self, sig_vector_fname):
		"""
		12-31-05
			
		"""
		sys.stderr.write("Getting no_of_datasets...")
		reader = csv.reader(open(sig_vector_fname, 'r'), delimiter='\t')
		edge_sig_vector = reader.next()
		no_of_datasets = len(edge_sig_vector)-2
		del reader
		sys.stderr.write("Got no_of_datasets.\n")
		return no_of_datasets
	
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
	
	def fillEdgeSigMatrix(self, communicator, free_computing_nodes, sig_vector_fname, edge_sig_vector_queue, \
		no_of_datasets, min_sup, max_sup):
		"""
		08-24-05
			similar to outputEdgeData() but read it into edge_sig_matrix, instead of output
		08-31-05
			Because Numeric.Int is only 32 bits. So it can't encode occurrence_vector and send it.
			If sending occurrence_vector directly, it will blow the message size. So send block by block
		12-30-05
			get edge data from a sig_vector_fname(in gene_no form), instead of edge_table
		12-31-05
			edge_sig_vector_queue replaces edge_sig_matrix
		01-01-06
			use block_size instead of 5000
		"""
		sys.stderr.write("Node %s getting edge recurrence_array...\n"%communicator.rank)
		node_rank = communicator.rank
		block_size = 10000
		edge_sig_block = Numeric.zeros((block_size, 2+no_of_datasets), Numeric.Int)
		if node_rank == 0:		
			reader = csv.reader(open(sig_vector_fname, 'r'), delimiter='\t')
			counter = 0
			for row in reader:
				edge_sig_vector = map(int, row)
				frequency = sum(edge_sig_vector[2:])
				if frequency>=min_sup and frequency<=max_sup:
					edge_sig_block[counter] = edge_sig_vector
					counter += 1
					if counter == block_size:
						for node in free_computing_nodes:	#send it to everyone
							communicator.send(edge_sig_block, node, 0)
							if self.debug:
								sys.stderr.write("a block sent to %s.\n"%node)
						counter = 0	#reset the counter
				"""
				elif self.debug:
						sys.stderr.write("edge_sig_vector %s, not included.\n"%(repr(edge_sig_vector)))
				"""
			del reader
			if counter>0 and counter<block_size:	#the last block is not empty and not sent.
				edge_sig_block[counter][0] = -1	#-1 is used to mark the end of this incomplete edge_sig_block
				for node in free_computing_nodes:	#send it to everyone
					communicator.send(edge_sig_block, node, 0)
					"""
					if self.debug:
						sys.stderr.write("a block sent to %s.\n"%node)
					"""
			#tell each node to exit the loop
			edge_sig_block[0,0] = -1
			for node in range(1, communicator.size):	#send it to everyone
				communicator.send(edge_sig_block, node, 0)
		else:
			edge_sig_block, source, tag, count = communicator.receive(edge_sig_block, 0, 0)	#get data from node 0
			while edge_sig_block:
				if edge_sig_block[0,0]==-1:
					"""
					if debug:
						sys.stderr.write("node %s breaked with %s edges.\n"%(node_rank, edge_sig_vector_queue.qsize()))
					"""
					break
				else:
					for edge_sig_vector in edge_sig_block:
						if edge_sig_vector[0]==-1:	#reach the end of this edge_sig_block
							break
						if edge_sig_vector[0]<edge_sig_vector[1]:	#in ascending order
							edge_sig_vector_queue.put([edge_sig_vector[0], edge_sig_vector[1], encodeOccurrenceBv(edge_sig_vector[2:])])
						else:
							edge_sig_vector_queue.put([edge_sig_vector[1], edge_sig_vector[0], encodeOccurrenceBv(edge_sig_vector[2:])])
				edge_sig_block, source, tag, count = communicator.receive(edge_sig_block, 0, 0)	#get data from node 0
			edge_sig_vector_queue.put(-1)	#12-31-05 the exit signal for the other thread
		sys.stderr.write("Node %s edge recurrence_array fetching Done\n"%node_rank)
	
	def createOffsetList(self, communicator, inputfile, no_of_nodes):
		"""
		08-06-05
			create an offset list to separate the inputfile
				the offset_list has no_of_nodes+1 entries, last entry is no_of_lines
				first entry is 0.
		12-31-05
			add some stderr output
		"""
		sys.stderr.write("node %s creating offset_list..."%communicator.rank)
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
		sys.stderr.write("done.\n")
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
		08-31-05
			the integer returned by encodeOccurrenceBv() could be 138-bit(human no_of_datasets)
			And Numeric.Int is only 32 bit. So Change edge_sig_matrix format.
		12-31-05
			no database connection any more
			2 threads on computing node
			
			(rank==0)
				--createOffsetList()
				--get_no_of_datasets()
			else:
				--start PatternFormThread
					--patternFormation()
					--outputCcFromEdgeList()
			
			--mpi_synchronize()
			--fillEdgeSigMatrix()
				--encodeOccurrenceBv()
			
			(rank==0)
				--receive node_outputfile
				--netmine_wrapper()
				--collect_and_merge_output()
				--uniqueSort()
			else:
				--join PatternFormThread
			
		"""
		communicator = MPI.world.duplicate()
		
		free_computing_nodes = range(1,communicator.size)	#exclude the 1st node
		edge_sig_vector_queue = Queue(self.queue_size)	#12-31-05	for thread
		if communicator.rank == 0:
			no_of_computing_nodes = len(free_computing_nodes)
			offset_list = self.createOffsetList(communicator, self.inputfile, no_of_computing_nodes)
			if self.debug:
				sys.stderr.write("offset_list: %s\n"%repr(offset_list))
			offset_list_pickle = cPickle.dumps(offset_list, -1)
			no_of_datasets = self.get_no_of_datasets(self.sig_vector_fname)
				#no_of_datasets is used in fillEdgeSigMatrix() and patternFormation()
			for node in free_computing_nodes:
				communicator.send(offset_list_pickle, node, 0)
				communicator.send(str(no_of_datasets), node, 0)
		else:
			data, source, tag = communicator.receiveString(0, 0)
			offset_list = cPickle.loads(data)			
			data, source, tag = communicator.receiveString(0, 0)
			no_of_datasets = int(data)	#take the data
			parameter_list =[self.inputfile, self.outputfile, offset_list, self.no_cc, edge_sig_vector_queue, no_of_datasets, self.debug]	
			if self.debug:
				sys.stderr.write("The common parameter_list is %s.\n"%repr(parameter_list))
			PatternFormThread_instance = PatternFormThread(communicator.rank, parameter_list)
			PatternFormThread_instance.start()
		
		mpi_synchronize(communicator)
		
		self.fillEdgeSigMatrix(communicator, free_computing_nodes, self.sig_vector_fname, \
			edge_sig_vector_queue, no_of_datasets, self.min_sup, self.max_sup)
		
		if communicator.rank == 0:
			#12-31-05 wait until of_name_list is full
			of_name_list = []
			while len(of_name_list)<no_of_computing_nodes:
				data, source, tag = communicator.receiveString(None, 1)
				of_name_list.append(data)
			#collecting
			intermediateFile = '%s.unsorted'%self.outputfile	#intermediateFile to store concatenated results
			netmine_wrapper_instance = netmine_wrapper()
			netmine_wrapper_instance.collect_and_merge_output(of_name_list, intermediateFile)
			self.uniqueSort(intermediateFile, self.outputfile)
		else:
			#12-31-05 wait for the thread to finish and return the node_outputfile
			PatternFormThread_instance.join()
			communicator.send(PatternFormThread_instance.node_outputfile, 0, 1)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:s:o:y:m:x:q:nbr", ["help", "hostname=", \
			"dbname=", "schema=", "inputfile=", "outputfile=", "min_sup", "max_sup=", \
			"no_cc", "debug", "report"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	inputfile = None
	sig_vector_fname = None
	outputfile = None
	min_sup = 0
	max_sup = 200
	queue_size = 3000000
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
		elif opt in ("-s", ):
			sig_vector_fname = arg
		elif opt in ("-o", "--outputfile"):
			outputfile = arg
		elif opt in ("-m", "--min_sup"):
			min_sup = int(arg)
		elif opt in ("-x", "--max_sup"):
			max_sup = int(arg)
		elif opt in ("-x", "--max_sup"):
			max_sup = int(arg)
		elif opt in ("-q",):
			queue_size = int(arg)
		elif opt in ("-n", "--no_cc"):
			no_cc = 1
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", "--report"):
			report = 1
	if inputfile and outputfile and sig_vector_fname:
		instance = MpiFromDatasetSignatureToPattern(hostname, dbname, schema, \
			inputfile, sig_vector_fname, outputfile, min_sup, max_sup, queue_size, no_cc, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
