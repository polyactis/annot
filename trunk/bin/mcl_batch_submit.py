#!/usr/bin/env python
"""
Usage: mcl_batch_submit.py [OPTION] MCLINPUTDIR

Option:
	MCLINPUTDIR a directory containing the mcl input files.
	-n ..., --jobnum=...	The starting number for job. default is 0.
	-p ..., --parameter=...	parameter passed to mcl, double quoted.
	-h, --help              show this help
	
Examples:
	mcl_batch_submit.py -n 11 gph_result/mcl_input
	mcl_batch_submit.py -n 11 -p "-I 5" gph_result/mcl_input
	
Description:
	Batch submit lots of mcl jobs. It needs a directory which stores
	all the mcl input files(prefix is mcl_input). Real work is done via
	mcl_batch_run.py. Latter iteratively calls mcl.

"""

import sys,os,math,getopt

class mcl_batch_submit:
	'''
	Partition a big mcl job into 
	'''
	def __init__(self, dir, no_of_job=0, parameter=''):
		self.dir = os.path.abspath(dir)
		if not os.path.isdir(self.dir):
			sys.stderr.write('%s not present\n'%self.dir)
			sys.exit(1)
		self.jobdir = os.path.join(os.path.expanduser('~'),'qjob')
		if not os.path.isdir(self.jobdir):
			os.makedirs(self.jobdir)
		self.f_list = os.listdir(self.dir)
		self.no_of_job = no_of_job
		self.parameter = parameter
	
	def submit(self):
		for file in self.f_list:
			if file.find('mcl_input') == 0:
				job_fname = os.path.join(self.jobdir, 'mcl'+repr(self.no_of_job)+'.sh')
				job_f = open(job_fname, 'w')
				out_block = '#!/bin/sh\n'		#this is '=' not '+='
				out_block += 'cd %s\n'%self.dir
				job_outfname = 'out%d'%self.no_of_job
				out_block += '~/script/annot/bin/mcl_batch_run.py %s %s %s\n'%\
					(file, job_outfname, self.parameter)
				job_f.write(out_block)
				job_f.close()
				wl = ['qsub', job_fname]
				os.spawnvp(os.P_WAIT, 'qsub', wl)
				self.no_of_job += 1

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hn:p:", ["help", "jobnum=", "parameter="])
	except:
		print __doc__
		sys.exit(2)
	
	jobnum = ''
	parameter = ''
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-n", "--jobnum"):
			jobnum = int(arg)
		elif opt in ("-p", "--parameter"):
			parameter = arg

	if jobnum != '' and len(args)==1:
		instance = mcl_batch_submit(args[0], jobnum, parameter)
		instance.submit()

	else:
		print __doc__
		sys.exit(2)
