#!/usr/bin/env mpipython
"""
Usage: lam_sub.py [OPTIONS] ProgramAndArguments

Option:

	-r ...	--rr=...	rank range, 0(default) but range is also ok. i.e. 0-3 means 0,1,2,3
		This is used to specify the nodes to use.
	-h, --help              show this help
	
Examples:
	mpirun C lam_sub.py -r 3 ls
	
Description:
	This program is like qsub, but need to specify node and the only way to
	kill a job is to go to that node and kill it.

	It depends on batch_haiyan_lam.py(one function borrowed)
	
	It's better to use lamexec. Simpler and faster.
	i.e. lamexec n3 ls
	
"""


from Scientific import MPI
import sys, os, getopt
from batch_haiyan_lam import batch_haiyan_lam

class lam_sub:
	"""
	03-17-05
		program submit jobs on lam/mpi system, but unfortunately,
		lamexec can do the same job better.
	"""
	def __init__(self, rank_range, args):
		"""
		03-17-05
			
		"""
		self.rank_range = rank_range
		self.args = args
		self.batch_haiyan_lam_instance = batch_haiyan_lam()
	
	def run(self):
		"""
		03-17-05
			call _node_fire in batch_haiyan_lam
		"""
		#a communicator
		comm = MPI.world.duplicate()
		if len(self.rank_range) == 1:
			node_rank_list = self.rank_range
		elif self.rank_range[0] > comm.size or self.rank_range[1] > comm.size:
			sys.stderr.write("Error: invalid rank range: %s"%repr(self.rank_range))
			sys.exit(2)
		else:
			node_rank_list = range(self.rank_range[0], self.rank_range[1]+1)
		#the first one is the program_path
		program_path = self.args[0]
		parameter_list = [self.args]
		if comm.rank in node_rank_list:
			self.batch_haiyan_lam_instance._node_fire(comm, program_path, parameter_list, node_rank_list)


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hr:", ["help","rr="])
	except:
		print __doc__
		sys.exit(2)
	rank_range = [0]
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-r", "--rr"):
			rank_range = arg.split('-')
			rank_range = map(int, rank_range)
		

	instance = lam_sub(rank_range, args)
	instance.run()
