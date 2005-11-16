#!/usr/bin/env mpipython
"""
Usage: MpiConfirmMatchResult.py -i -o [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database(transfac)
	-i ...,	inputdir, TRANSFAC match output files
	-o ...,	output file
	-p ...,	profile file which is used by match
	-e ...,	max_mis_match_perc, 0.1(default)
	-n ...,	min_no_of_mismatches, 1(default)
	-x ...,	max_esc_length, 6(default), sites less than it require full match
	-v ...,	message size(10,000,000, default)
	-b,	debug version.
	-r,	enable report flag
	-h, --help	Display the usage infomation.
	
Examples:
	mpirun -np 20 -machinefile ~/hostfile /usr/bin/mpipython ~/script/annot/bin/MpiConfirmMatchResult.py
	-i -o -p -e -n -x
	
Description:
	This program tries to confirm the TRANSFAC match output by directly counting
	how many mismatches exist between the found binding site and the original
	sequences(or consensus) for the PWM.
	
"""

import sys, os, getopt, csv, cPickle, fileinput, cStringIO
sys.path += [os.path.expanduser('~/script/annot/bin')]
sys.path += [os.path.expanduser('~/script/transfac/src')]
from Scientific import MPI
from codense.common import mpi_synchronize, db_connect, output_node, \
	computing_node, input_node, get_mt_id_set_from_profile, get_mt_id2sites_ls
from transfacdb import match_block_iterator
from Bio.Seq import Seq
from Bio.Data import IUPACData

class MpiConfirmMatchResult:
	def __init__(self,hostname='zhoudb', dbname='graphdb', schema=None, inputdir=None,\
		output_file=None, profile_filename=None, max_mis_match_perc=0.1, min_no_of_mismatches=1,\
		max_esc_length=6, message_size=10000000, debug=0, report=0):
		"""
		11-15-05
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.inputdir = inputdir
		self.output_file = output_file
		self.profile_filename = profile_filename
		self.max_mis_match_perc = float(max_mis_match_perc)
		self.min_no_of_mismatches = int(min_no_of_mismatches)
		self.max_esc_length = int(max_esc_length)
		self.message_size= int(message_size)
		self.debug = int(debug)
		self.report = int(report)
	
	def get_match_output_header(self, input_filename):
		sys.stderr.write("Getting match_output_header...\n")
		match_output_header = ''
		inf = open(input_filename, 'r')
		for i in range(4):
			match_output_header += inf.readline()
		del inf
		sys.stderr.write("Done getting match_output_header.\n")
		return match_output_header
	
	def input_handler(self, parameter_list, message_size, report=0):
		if report:
			sys.stderr.write("Fetching stuff...\n")
		iter = parameter_list[0]
		block = []
		string_length = 0
		for match_block in iter:
			block.append(match_block)
			string_length += len(match_block)
			if string_length>=message_size:
				break
		if report:
			sys.stderr.write("Fetching done.\n")
		return block
	
	def get_no_of_mismatches_allowed(self, sequence, max_mis_match_perc, min_no_of_mismatches, max_esc_length):
		#derive the no_of_mismatches_allowed
		sequence_length = len(sequence)
		if sequence_length<=max_esc_length:
			no_of_mismatches_allowed = 0
		else:
			no_of_mismatches_allowed = int(sequence_length*max_mis_match_perc)
			if no_of_mismatches_allowed<min_no_of_mismatches:
				no_of_mismatches_allowed = min_no_of_mismatches
		"""
		if self.debug:
			sys.stderr.write("max_mis_match_perc: %s; min_no_of_mismatches: %s; max_esc_length: %s\n"%\
				(max_mis_match_perc, min_no_of_mismatches, max_esc_length))
			sys.stderr.write("no_of_mismatches_allowed: %s.\n"%no_of_mismatches_allowed)
		"""
		return no_of_mismatches_allowed
	
	def is_good_consensus(self, consensus, max_esc_length):
		"""
		11-16-05
			the consensus should at least have >= max_esc_length unambiguous letters. 
			otherwise, it's not significant at all.
		"""
		no_of_unambiguous_letters = 0
		for i in range(len(consensus)):
			if consensus[i] in 'ACGT':
				no_of_unambiguous_letters += 1
		return no_of_unambiguous_letters>=max_esc_length
	
	def get_no_of_mismatches_for_consensus(self, sequence, consensus, no_of_mismatches_allowed, max_esc_length):
			"""
			if self.debug:
				sys.stderr.write("sequence %s vs consensus: %s \n"%(sequence, consensus))
			"""
			if not self.is_good_consensus(consensus, max_esc_length):
				return False
			no_of_mismatches = 0
			#more strict for consensus
			no_of_mismatches_allowed = 0
			for i in range(len(sequence)):
				if sequence[i] not in IUPACData.ambiguous_dna_values[consensus[i]]:
					no_of_mismatches += 1
					if no_of_mismatches>no_of_mismatches_allowed:
						break
			"""
			if self.debug:
				sys.stderr.write("consensus %s with no_of_mismatches %s.\n"%(consensus, no_of_mismatches))
			"""
			return no_of_mismatches<=no_of_mismatches_allowed
	
	def get_no_of_mismatches_for_site(self, sequence, sites_ls, no_of_mismatches_allowed, max_esc_length):
			for site in sites_ls[1:]:
				no_of_mismatches = 0	#zero before each comparison
				"""
				if self.debug:
					sys.stderr.write("sequence %s vs site %s..\n"%(sequence, site))
				"""
				for i in range(len(site)):	#site could be shorter than sequence, trailing spaces stripped
					if site[i]!=' ' and sequence[i] != site[i]:	#' '(empty) means that position is not available in site, regarded as 'N'
						no_of_mismatches += 1
						if no_of_mismatches>no_of_mismatches_allowed:
							break
				"""
				if self.debug:
					sys.stderr.write("site %s with no_of_mismatches %s.\n"%(site, no_of_mismatches))
				"""
				#if the number of mismatches between site[i] and sequence is under threshold, we got it.
				if len(site)<=max_esc_length and no_of_mismatches==0:
					#AGAIN, site could be shorter, so check its length first
					break
				elif no_of_mismatches<=no_of_mismatches_allowed:
					break
			return no_of_mismatches<=no_of_mismatches_allowed
	
	def is_site_confirmed(self, mt_id2sites_ls, line, max_mis_match_perc, min_no_of_mismatches, max_esc_length):
		### 1st parse (copied from transfacdb.py
		ls = line[:-1].split('|')
		mt_id = ls[0].strip()	#remove spaces
		bs_disp_start_strand = ls[1].strip()
		#bs_disp_start = int(bs_disp_start_strand[:-3])
		strand = bs_disp_start_strand[-2]
		#core_similarity_score = float(ls[2])
		#matrix_similarity_score = float(ls[3])
		sequence = ls[4].strip()
		
		if strand=='-':	#take the reverse_compliment()
			seq = Seq(sequence)
			sequence = seq.reverse_complement().tostring()
			"""
			if self.debug:
				sys.stderr.write("Strand is -, need reverse_compliment() from %s to %s.\n"%(seq.data, sequence))
			"""
		#transform it into upper case
		sequence = sequence.upper()
		
		no_of_mismatches_allowed = self.get_no_of_mismatches_allowed(sequence, max_mis_match_perc, \
			min_no_of_mismatches, max_esc_length)
		
		#check the no_of_mismatches
		sites_ls = mt_id2sites_ls[mt_id]
		if sites_ls[0] == 0:	#it's consensus
			return self.get_no_of_mismatches_for_consensus(sequence, sites_ls[1], no_of_mismatches_allowed,\
				max_esc_length)
		elif sites_ls[0] == 1:	#it's the sequence where the consensus is derived
			return self.get_no_of_mismatches_for_site(sequence, sites_ls, no_of_mismatches_allowed,\
				max_esc_length)
		else:
			sys.stderr.write("Wrong type of sites_ls of mt_id2sites_ls: %s.\n"%sites_ls[0])
			return None
		
	
	def computing_handler(self, communicator, data, parameter_list):
		node_rank = communicator.rank
		sys.stderr.write("Node no.%s working...\n"%node_rank)
		data = cPickle.loads(data)
		mt_id2sites_ls, max_mis_match_perc, min_no_of_mismatches, max_esc_length = parameter_list
		match_confirmed_results = ''
		for match_block in data:
			match_block = cStringIO.StringIO(match_block)
			seq_id_line = match_block.readline()
			match_confirmed_result = ''
			for line in match_block:	#now start the binding site info
				if self.is_site_confirmed(mt_id2sites_ls, line, max_mis_match_perc, min_no_of_mismatches, max_esc_length):
					match_confirmed_result += line
			if match_confirmed_result:	#not nothing
				match_confirmed_results += seq_id_line+match_confirmed_result
		sys.stderr.write("Node no.%s done.\n"%node_rank)
		return match_confirmed_results
	
	def output_handler(self, communicator, parameter_list, data):
		outf = parameter_list[0]
		data = cPickle.loads(data)
		outf.write(data)
	
	def run(self):
		"""
		11-16-05
			
			--computing_handler()
				--is_site_confirmed()
					--get_no_of_mismatches_allowed()
					--get_no_of_mismatches_for_consensus()
						--is_good_consensus()
					--get_no_of_mismatches_for_site()
		"""
		communicator = MPI.world.duplicate()
		node_rank = communicator.rank	
		free_computing_nodes = range(1,communicator.size-1)	#exclude the last node
		if node_rank == 0:
			(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
			if self.profile_filename:
				mt_id_set = get_mt_id_set_from_profile(self.profile_filename)
			else:
				mt_id_set = None
			mt_id2sites_ls = get_mt_id2sites_ls(curs, mt_id_set)
			mt_id2sites_ls_pickle = cPickle.dumps(mt_id2sites_ls, -1)
			for node in free_computing_nodes:	#send it to the computing_node
				communicator.send(mt_id2sites_ls_pickle, node, 0)
			
			input_files = os.listdir(self.inputdir)
			for i in range(len(input_files)):	#attach the directory path to the files
				input_files[i] = os.path.join(self.inputdir, input_files[i])
			#the following infomation is just header info inserted into the top of the output_file
			match_output_header = self.get_match_output_header(input_files[0])
			communicator.send(match_output_header, communicator.size-1, 0)
		elif node_rank in free_computing_nodes:
			data, source, tag = communicator.receiveString(0, 0)
			mt_id2sites_ls = cPickle.loads(data)	#take the data
		elif node_rank==communicator.size-1:
			outf = open(self.output_file, 'w')
			match_output_header, source, tag = communicator.receiveString(0, 0)
			outf.write(match_output_header)
			
		mpi_synchronize(communicator)
		if node_rank == 0:
			aggregated_inf = fileinput.input(input_files)
			iter = match_block_iterator(aggregated_inf)
			parameter_list = [iter]
			input_node(communicator, parameter_list, free_computing_nodes, self.message_size, self.report, \
				input_handler=self.input_handler)
			del iter, aggregated_inf
		elif node_rank in free_computing_nodes:
			parameter_list = [mt_id2sites_ls, max_mis_match_perc, min_no_of_mismatches, max_esc_length]
			computing_node(communicator, parameter_list, self.computing_handler, report=self.report)
		elif node_rank==communicator.size-1:
			parameter_list = [outf]
			output_node(communicator, free_computing_nodes, parameter_list, self.output_handler, self.report)
			del outf


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:o:p:e:n:x:v:br", ["help", \
			"hostname=", "dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = 'transfac'
	inputdir = None
	output_file = None
	profile_filename = None
	max_mis_match_perc = 0.1
	min_no_of_mismatches = 1
	max_esc_length = 6
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
		elif opt in ("-i"):
			inputdir = arg
		elif opt in ("-o"):
			output_file = arg
		elif opt in ("-p"):
			profile_filename = arg
		elif opt in ("-e"):
			max_mis_match_perc = float(arg)
		elif opt in ("-n"):
			min_no_of_mismatches = int(arg)
		elif opt in ("-x"):
			max_esc_length = int(arg)
		elif opt in ("-v"):
			message_size = int(arg)
		elif opt in ("-b"):
			debug = 1
		elif opt in ("-r"):
			report = 1
	if inputdir and output_file:
		instance = MpiConfirmMatchResult(hostname, dbname, schema, inputdir, output_file, profile_filename, \
			max_mis_match_perc,  min_no_of_mismatches, max_esc_length, message_size, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
