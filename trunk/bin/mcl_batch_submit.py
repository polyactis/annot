#!/usr/bin/env python
import sys,os,math

class mcl_batch:
	'''
	Partition a big mcl job into 
	'''
	def __init__(self, dir, no_of_jobs, parameter=''):
		self.dir = os.path.abspath(dir)
		if not os.path.isdir(self.dir):
			sys.stderr.write('%s not present\n'%self.dir)
			sys.exit(1)
		self.jobdir = os.path.join(os.path.expanduser('~'),'qjob')
		if not os.path.isdir(self.jobdir):
			os.makedirs(self.jobdir)
		self.f_list = os.listdir(self.dir)
		self.no_of_jobs = int(no_of_jobs)
		self.parameter = parameter
	
	def submit(self):
		no_of_files = len(self.f_list)
		files_per_job = int( math.ceil(no_of_files/float(self.no_of_jobs)) )

		for i in range(self.no_of_jobs):
			job_fname = os.path.join(self.jobdir, 'mcl'+repr(i+1)+'.sh')
			job_f = open(job_fname, 'w')
			job_f.write('#!/bin/sh\n')
			job_f.write('cd %s\n'%self.dir)
			for j in range(files_per_job):
				index = i*files_per_job + j
				if index == no_of_files:
					break
				job_f.write('mcl %s %s\n'%(self.f_list[index], self.parameter,))
			job_f.close()
			wl = ['qsub', job_fname]
			os.spawnvp(os.P_WAIT, 'qsub', wl)

if __name__ == '__main__':
	def helper():
		sys.stderr.write('\
		argv[1] is directory containing mcl input files.\n\
		argv[2] is the no_of_jobs you want to submit.\n\
		argv[3] is the parameter passed to mcl. You can ignore it.\n')
		
	if len(sys.argv) == 4:
		instance = mcl_batch(sys.argv[1], sys.argv[2], sys.argv[3])
	elif len(sys.argv) == 3:
		instance = mcl_batch(sys.argv[1], sys.argv[2])
	else:
		helper()
		sys.exit(1)
	instance.submit()
