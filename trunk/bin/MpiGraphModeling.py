#!/usr/bin/env mpipython
"""
Usage: MpiGraphModeling.py -i INPUT_DIR -o OUTPUT_DIR [OPTION]

Option:
	-i ...,	--input_dir=...	the input directory
	-o ..., --output_dir=...	the output directory, outfname = gph_infname
	-p .., --p_value_cut_off=...	p_value significance cutoff,0.01(default)
	-c ..., --cor_cut_off=...	correlation cutoff, 0.6(default)
		if p_value_cut_off=0, this cut_off is used instead.
	-d ..., --max_degree=...	maximum degree of freedom(#columns-2), 10000,(default).
	-t ..., --top_percentage=...	0.01(default).
		if p_value_cut_off=0 and cor_cut_off=0, top_percentage is used to select edges.
	-b, --debug	debug version.
	-l, --leave_one_out	leave_one_out.
	-h, --help	Display the usage infomation.
	
Examples:
	mpirun.lam N MpiGraphModeling.py -i datasets/hs -o gph_result/hs
	
	mpirun.lam N MpiGraphModeling.py -l -i datasets/hs -o gph_result/hs
		leave_one_out version
	
Description:
	The MPI version of graph_modeling, which replaces batch_graph_cc.py
	
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import sys, os, getopt
from sets import Set
from Numeric import array, Float, take, zeros, Int
from Scientific import MPI
from graph import graph_modeling

class MpiGraphModeling:
	"""
	05-13-05
		the mpi version of graph_modeling.cc
	"""
	def __init__(self, input_dir=None, output_dir=None, p_value_cut_off='0.01', cor_cut_off='0.6', \
		max_degree='10000', top_percentage="0.01", leave_one_out=0, debug=0):
		"""
		05-13-05
			parameters to be passed to graph_modeling must be in string form.
		"""
		self.input_dir = input_dir
		self.output_dir = output_dir
		self.p_value_cut_off = p_value_cut_off
		self.cor_cut_off = cor_cut_off
		self.max_degree = max_degree
		self.top_percentage = top_percentage
		self.leave_one_out = int(leave_one_out)
		self.debug = int(debug)
		
	
	def schedule_jobs(self, communicator, bin_path, input_dir, output_dir, parameters):
		"""
		05-13-05
			the schedule mechanism is copied from MpiBiclustering.py
			
			if #nodes > #jobs, tell those nodes to break the listening loop.
		05-16-05
			record the exit code to conquer the weird problem(os.system() fails on the first time on some machines)
			#while exit_code not zero, means os.system() fails, rerun it.
		"""
		node_rank = communicator.rank
		if node_rank == 0:
			#if the output_dir doesn't exist
			if not os.path.isdir(output_dir):
				os.makedirs(output_dir)
			files = os.listdir(input_dir)
			seed_utilized = Set()
			for node in range(1, communicator.size):
				if len(files)==0:	#if #nodes > #jobs, tell those nodes to break their listening loop.
					stop_signal = "-1"
					communicator.send(stop_signal, node, 0)	#no more jobs, stop that node,
					if self.debug:
						sys.stderr.write("node %s stopped.\n"%node)
				else:
					input_file = files.pop(0)	#the first item poped first.
					communicator.send(input_file, node, 0)
					if self.debug:
						sys.stderr.write("Node %s schedule a job to %s\n"%(node_rank, node))
					seed_utilized.add(node)
			
			received_value, source, tag = communicator.receiveString(None, None)	#listen
			while received_value:		#??check what the received_value is
				if len(files) == 0:	#first check if there're still files left, otherwise pop(0) raises error.
					stop_signal = "-1"
					communicator.send(stop_signal, source, 0)	#no more jobs, stop that node,
					if self.debug:
						sys.stderr.write("node %s stopped.\n"%source)
					seed_utilized.remove(source)
					if len(seed_utilized) == 0:	#all seed used have finished their jobs
						break
				else:
					input_file = files.pop(0)
					if input_file:
						communicator.send(input_file, source, 0)	#more jobs
						if self.debug:
							sys.stderr.write("Node %s get one more job\n"%source)
				received_value, source, tag = communicator.receiveString(None, None)	#listen
		else:
			received_data, source, tag = communicator.receiveString(0, None)	#get data from node 0,
				#04-24-05 the array is one-dimension no matter what dimension the original array is
			while received_data:
				if received_data=="-1":	#it's array([0.0]), stopping signal, don't use 'received_data==array([0.0])' to judge.
					if self.debug:
						sys.stderr.write("node %s breaked.\n"%node_rank)
					break
				else:
					input_file = received_data
					output_file = 'gph_%s'%received_data
					input_file = os.path.join(input_dir, input_file)	#absolute path
					output_file = os.path.join(output_dir, output_file)	#absolute path
					jobrow = '%s -o %s %s %s'%(bin_path, output_file, parameters, input_file)
					sys.stderr.write("node %s working on %s...\n"%(node_rank, received_data))
					exit_code = os.system(jobrow)
						#05-16-05 record the exit code to conquer the weird problem(os.system() fails on the first time on some machines)
						#while exit_code not zero, means os.system() fails, rerun it.
					while exit_code:
						exit_code = os.system(jobrow)
					sys.stderr.write("node %s work on %s done, exit_code: %s.\n"%(node_rank, received_data, exit_code))
					communicator.send("finished", 0, node_rank)
					
				received_data, source, tag = communicator.receiveString(0, None)	#get data from node 0
	
	def run(self):
		"""
		05-13-05
			
			--schedule_jobs()
		"""
		communicator = MPI.world.duplicate()
		
		homedir = os.path.expanduser('~')
		bin_path = os.path.join(homedir,'script/annot/bin/graph/graph_modeling')
		parameters = '-p %s -c %s -d %s -t %s'%(self.p_value_cut_off, self.cor_cut_off, self.max_degree, self.top_percentage)
		if self.leave_one_out:
			parameters += ' -l'
		
		self.schedule_jobs(communicator, bin_path, self.input_dir, self.output_dir, parameters)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "input_dir=", "output_dir=", "p_value_cut_off=", "cor_cut_off=", "max_degree=",\
		"top_percentage=", "leave_one_out", "debug"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hi:o:p:c:d:t:lb", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	input_dir = None
	output_dir = None
	p_value_cut_off = "0.01"
	cor_cut_off = "0.6"
	max_degree = "10000"
	top_percentage = "0.01"
	leave_one_out = 0
	debug = 0
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-i", "--input_dir"):
			input_dir = arg
		elif opt in ("-o", "--output_dir"):
			output_dir = arg
		elif opt in ("-p", "--p_value_cut_off"):
			p_value_cut_off = arg
		elif opt in ("-c", "--cor_cut_off"):
			cor_cut_off = arg
		elif opt in ("-d", "--max_degree"):
			max_degree = arg
		elif opt in ("-t", "--top_percentage"):
			top_percentage = arg
		elif opt in ("-l", "--leave_one_out"):
			leave_one_out = 1
		elif opt in ("-b", "--debug"):
			debug = 1

	if input_dir and output_dir:
		instance = MpiGraphModeling(input_dir, output_dir, p_value_cut_off, cor_cut_off, max_degree, \
			top_percentage, leave_one_out, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
