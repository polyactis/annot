#!/usr/bin/env mpipython
"""
Usage: graph_merge_lam.py [OPTION] INPUTDIR OUTPUTFILE

Option:
	INPUTDIR is the directory containing all the graph files in gspan format.
	OUTPUTFILE is the file to store the merged graphs also in gspan format.
	-s ..., --support=...	minimum support for the edge to be kept. 5(default)
	-t ..., --threshold=...	the size of graph_dict, 10000,000(default) 26.5% of 6G memory
	-r ..., --max_rank=...	the maximum rank of nodes to initiate the program, size-1(default)
	-h, --help              show this help
	
Examples:
	graph_merge_lam.py -s 6 gph_result/sc/ gph_result/sc_mcl/mcl_gph_dataset1

Description:
	This program merges all the graphs in gspan format, generated by
	graph_reorganize.py. The ouput is also in gspan format.
	After this, either run reverse+kMax or gspan2mcl_input.py+mcl to
	get the dense clusters.
	
	03-13-05
		lam/mpi version of graph_merge.py

"""


from Scientific import MPI
import sys, os, re, getopt
import Numeric

class graph_merge:
	'''
	03-13-05
		the lam/mpi version of graph_merge
	03-16-05
		add max_rank to specify which nodes to start off the program
		or how many nodes used as memory bucket.
	'''
	def __init__(self, support, threshold, max_rank, dir, ofname):
		self.support = int(support)
		self.threshold = int(threshold)
		self.max_rank = max_rank
		self.dir = dir
		self.ofname = ofname
		#data structure to store the merged graph. key is the edge.
		#value is the recurrence
		self.graph_dict = {}
	
	def check_edge(self, edge, communicator, half_full_machine_no, dict_size_of_half_full_machine):
		"""
		03-13-05
			input: an edge
			output: either '-1'(incremented) or a threshold_indicator_list
		12-26-05 speed up
		01-29-06
			fix a bug, "no break after '1' and '2' signal			
			
			check which node has this edge, 
			return '1': that node has it, increment its counter in its graph_dict
			return '2': not have it but memory not full and add it
			return '0': not have it and memory full and can't add it
		"""
		for dest in range(1, half_full_machine_no+1):
			communicator.send(edge, dest, dest)
			data, source, tag, count = communicator.receive(Numeric.Int, dest, None)
			if data[0]==2:
				if dest!=half_full_machine_no:
					sys.stderr.write("Error: machine, %s, adding this edge, %s is not half_full_machine_no, %s.\n"\
						%(dest, edge, half_full_machine_no))
				dict_size_of_half_full_machine += 1
				break
			elif data[0] == 1:
				break
		if data[0]==0:
			sys.stderr.write("Error: edge %s neither inserted nor incremented.\n"%edge)
		return dict_size_of_half_full_machine

	def edge_loadin(self, dir, communicator, max_rank):
		"""
		12-26-05
			design a way to speed up program
		"""
		#output block put before edges
		first_block = ''
		
		files = os.listdir(dir)
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
		
		half_full_machine_no = 1	#12-26-05	the initial machine without memory being full
		dict_size_of_half_full_machine = 0	#12-26-05
		
		for f in files:
			pathname = os.path.join(dir, f)
			sys.stderr.write("%d/%d:\t%s\n"%(files.index(f)+1,len(files),f))
			file_no = files.index(f)+1
			inf = open(pathname, 'r')
			for line in inf:
				if line[0] == 'e':
					#edge here, like 'e 3807 3859 0.804645'
					line_list = line[:-1].split()
					vertex1 = int(line_list[1])
					vertex2 = int(line_list[2])
					if vertex1 <= vertex2:
						edge = Numeric.array([vertex1,vertex2])
					else:
						edge = Numeric.array([vertex2, vertex1])
					dict_size_of_half_full_machine = self.check_edge(edge, communicator, half_full_machine_no, dict_size_of_half_full_machine)
					if dict_size_of_half_full_machine>=threshold:	#12-26-05	increase the machine no
						half_full_machine_no += 1
						dict_size_of_half_full_machine = 0	#12-26-05 reset the dict size for half_full_machine_no
						if half_full_machine_no>max_rank:
							sys.stderr.write("Error: memory used up on all machines(current half_full_machine_no: %s.\n"%half_full_machine_no)
							break
				elif file_no == 1:
					first_block += line
			inf.close()
		
		return first_block

	def output(self, ofname, first_block, support, communicator, max_rank):
		"""
		03-13-05
			input: first_block, support, communicator
			no output:
		12-27-05
			change the type of communication message
		
			first output the first_block,
			then tell each node to output its graph_dict given a support
		"""
		#output the preceding block first
		of = open(ofname, 'w')
		of.write(first_block)
		of.close()
		output_signal = Numeric.array([-1,-1])	#12-27-05 No gene is named -1
		for dest in range(1, max_rank+1):
			communicator.send(output_signal, dest, dest)
			return_value, source, tag, count = communicator.receive(Numeric.Int, dest, None)	#12-27-05
			if return_value[0]==1:
				print "%s finished its output"%(source)
			else:
				print "%s encounted error: %s in its output"%(source, return_value)
	
	def node_output(self, ofname, support, communicator):
		"""
		03-13-05
			input: ofname, support, communicator
			no output:
			
			This function is for all secondary nodes.
			It outputs its graph_dict.
		"""
		sys.stderr.write("%s is outputing its graph_dict..."%communicator.rank)
		#NOTE: It's 'a' not 'w'. 'w' overwrites output() of the primary node.
		of = open(ofname, 'a')
		for edge in self.graph_dict.keys():
			recurrence = self.graph_dict[edge]
			if recurrence >= support:
				of.write("e %d %d %d\n"%(edge[0], edge[1], recurrence))
		#don't forget to close it
		of.close()
		sys.stderr.write("Done\n")
		
	def node_loop(self, ofname, support, threshold, communicator):
		"""
		03-13-05
			This is a message capturing loop for each secondary node.
			They function according to the message.
		12-25-05
			use break instead of sys.exit(0)
			sys.exit(0) causes the main node to be dead
		12-26-05
			merge the functinality of 'not full' and 'add edge'
		12-27-05
			change the communication type
		12-27-05
			define signal ahead, don't make it on the fly when sending message
		"""
		good_signal = Numeric.array([1])
		add_signal = Numeric.array([2])
		bad_signal = Numeric.array([0])
		while 1:
			data, source, tag, count = communicator.receive(Numeric.Int, 0, None)	#12-27-05
			if data[0]==-1:	#12-27-05
				self.node_output(ofname, support, communicator)
				communicator.send(good_signal,0,1)	#12-27-05
				break
			else:
				#create an edge tupple
				edge = (data[0], data[1])
				if edge in self.graph_dict:
					#present, increment its counter
					self.graph_dict[edge] += 1
					communicator.send(good_signal,0, 3)
				else:
					#not present
					if len(self.graph_dict)<threshold:
						#not full and add this edge
						self.graph_dict[edge] = 1
						communicator.send(add_signal,0,4)
					else:
						#full
						communicator.send(bad_signal,0,5)

	def run(self):
		"""
		03-13-05
			the central coordinator
			
			rank==0
				--edge_loadin()
					--check_edge()
					--add_edge()
				--ouput()
			rank==1
				--node_loop()
					--node_output()
		"""
		communicator = MPI.world.duplicate()
		#set the default max_rank
		if self.max_rank == None or self.max_rank >=communicator.size:
			self.max_rank = communicator.size-1
		if communicator.rank == 0:
			first_block = self.edge_loadin(self.dir, communicator, self.max_rank)
			self.output(self.ofname, first_block, self.support, communicator, self.max_rank)
		elif communicator.rank<=self.max_rank:
			self.node_loop(self.ofname, self.support, self.threshold, communicator)
			
	
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "s:t:r:h", ["support=", "threshold=", "max_rank=", "help"])
	except:
		print __doc__
		sys.exit(2)
	
	support = 5
	threshold = 10000000
	max_rank = None
	for opt, arg in opts:
		if opt in ("-s", "--support"):
			support = int(arg)
		elif opt in ("-t", "--threshold"):
			threshold = int(arg)
		elif opt in ("-r", "--max_rank"):
			max_rank = int(arg)
		elif opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)

			
	if len(args) == 2:
		instance = graph_merge(support, threshold, max_rank, args[0], args[1])
		instance.run()
	else:
		print __doc__
		sys.exit(2)
