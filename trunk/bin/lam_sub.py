#!/usr/bin/env mpipython
"""
Usage: lam_sub.py [OPTIONS] ProgramAndArguments

Option:

	-r ...	--rr=...	rank range, (IGNORE)
		NOTE use 'n#' by mpirun or mpiexec to control the node.
	-h, --help              show this help
	
Examples:
	#run the ls on n3
	mpirun n3 lam_sub.py ls
	
Description:
	This program is like qsub, but need to specify node and the only way to
	kill a job is to go to that node and kill it.
	
	It's better to use lamexec. Simpler and faster.
	i.e. lamexec n3 ls
	
	But the good thing of this program is that it keeps track of the running by mpitask
	because it invokes MPI functions(the communicator).
"""


from Scientific import MPI
import sys, os, getopt

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
	
	def run(self):
		"""
		03-17-05
			
		"""
		#a communicator
		comm = MPI.world.duplicate()
		#the first one is the program_path
		program_path = self.args[0]
		print "node %s running %s..."%(comm.rank, program_path)
		os.execvp(program_path, self.args)


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hr:", ["help","rr="])
	except:
		print __doc__
		sys.exit(2)
	rank_range = None
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-r", "--rr"):
			rank_range = arg.split('-')
			rank_range = map(int, rank_range)
		

	instance = lam_sub(rank_range, args)
	instance.run()
