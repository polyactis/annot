#!/usr/bin/env python
import sys,os,math

class mcl_batch_submit:
	'''
	Partition a big mcl job into 
	'''
	def __init__(self, dir, parameter=''):
		self.dir = os.path.abspath(dir)
		if not os.path.isdir(self.dir):
			sys.stderr.write('%s not present\n'%self.dir)
			sys.exit(1)
		self.jobdir = os.path.join(os.path.expanduser('~'),'qjob')
		if not os.path.isdir(self.jobdir):
			os.makedirs(self.jobdir)
		self.f_list = os.listdir(self.dir)
		self.parameter = parameter
	
	def submit(self):
		no_of_jobs = 0
		for file in self.f_list:
			if file.find('mcl_input') == 0:
				job_fname = os.path.join(self.jobdir, 'mcl'+repr(no_of_jobs)+'.sh')
				job_f = open(job_fname, 'w')
				out_block = '#!/bin/sh\n'		#this is '=' not '+='
				out_block += 'cd %s\n'%self.dir
				job_outfname = 'out%d'%no_of_jobs
				out_block += '~/script/annot/bin/mcl_batch_run.py %s %s %s\n'%\
					(file, job_outfname, self.parameter)
				job_f.write(out_block)
				job_f.close()
				wl = ['qsub', job_fname]
				os.spawnvp(os.P_WAIT, 'qsub', wl)
				no_of_jobs += 1

if __name__ == '__main__':
	def helper():
		sys.stderr.write('\
		argv[1] is directory containing mcl input files.\n\
		argv[2] is the parameter passed to mcl. Ignore to use default.\n\
			Parameter should be double quoted("").\n')
		
	if len(sys.argv) == 3:
		instance = mcl_batch_submit(sys.argv[1], sys.argv[2])
	elif len(sys.argv) == 2:
		instance = mcl_batch_submit(sys.argv[1])
	else:
		helper()
		sys.exit(1)
	instance.submit()
