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
			out_block = '#!/bin/sh\n'		#this is '=' not '+='
			out_block += 'cd %s\n'%self.dir
			job_outfname = 'out%d'%i
			out_block += 'touch %s\n'%job_outfname
			for j in range(files_per_job):
				index = i*files_per_job + j
				if index == no_of_files:
					break
				mcl_outfname = 'mcl_%s'%self.f_list[index]
				parameter = ' -o %s '%mcl_outfname		#this is for directing mcl's output.
				parameter += self.parameter				#adding up the user_passed parameter
				out_block += 'echo "(splat_id %s )"  >>%s\n'%(self.f_list[index], job_outfname)
				out_block += 'echo "(parameter %s )"  >>%s\n'%(self.parameter, job_outfname)
				out_block += 'mcl %s %s\n'%(self.f_list[index], parameter)
				out_block += 'cat %s >>%s\n'%(mcl_outfname, job_outfname)
				out_block += 'rm %s\n'%mcl_outfname
			job_f.write(out_block)
			job_f.close()
			wl = ['qsub', job_fname]
			os.spawnvp(os.P_WAIT, 'qsub', wl)

if __name__ == '__main__':
	def helper():
		sys.stderr.write('\
		argv[1] is directory containing mcl input files.\n\
		argv[2] is the no_of_jobs you want to submit.\n\
		argv[3] is the parameter passed to mcl. Default is none. Parameter should be double quoted("").\n')
		
	if len(sys.argv) == 4:
		instance = mcl_batch(sys.argv[1], sys.argv[2], sys.argv[3])
	elif len(sys.argv) == 3:
		instance = mcl_batch(sys.argv[1], sys.argv[2])
	else:
		helper()
		sys.exit(1)
	instance.submit()
