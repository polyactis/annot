#!/usr/bin/env python
"""
Usage: graph_construct_batch_submit.py [OPTIONS]

Option:
	-d ..., --data_dir=...,	DATA_DIR is the directory where all the datasets are stored
	-g ..., --gph_dir=...,	GPH_DIR is the directory to store the graph results
	-n ..., --no_of_jobs=...,	specify how many jobs you want to submit
	-r, --restore	cleanup the directory
	-h, --help              show this help
	
Examples:
	Two kinds of usage:
	graph_construct_batch_submit.py -n 10 -d datasets/hs -g gph_result/hs
		submit 10 graph_contruct jobs.
	graph_construct_batch_submit.py -r -d datasets/hs
		restore the datasets into directory 'datasets/hs'
Description:
	This program partitions a big batch_graph_cc job into several subjobs,
	submits them via qsub.
	It moves all the files into several sub directories and calls batch_graph_cc.py.

Warning:
	If you want submit the same job again. Please run the second example first
	to restore the datasets in the data_dir.
	
"""


import sys,os,math,getopt

class batch_submit:
	'''
	Partition a big batch_graph_cc job into several small jobs
	'''
	def __init__(self, dir, dstdir, no_of_jobs, restore):
		self.dir = os.path.abspath(dir)
		if dstdir:
			self.dstdir = os.path.abspath(dstdir)
		if not os.path.isdir(self.dir):
			sys.stderr.write('%s not present\n'%self.dir)
			sys.exit(1)
		self.jobdir = os.path.join(os.path.expanduser('~'),'qjob')
		if not os.path.isdir(self.jobdir):
			os.makedirs(self.jobdir)
		self.f_list = os.listdir(self.dir)
		if no_of_jobs:
			self.no_of_jobs = int(no_of_jobs)
		self.restore = int(restore)
	
	def run(self):
		if self.restore==1:
			for f in self.f_list:
				path = os.path.join(self.dir, f)
				if os.path.isdir(path):
					sub_f_list = os.listdir(path)
					for fname in sub_f_list:
					#move the files to the parent directory
						src = os.path.join(path, fname)
						dst = os.path.join(self.dir, fname)
						os.rename(src, dst)
					os.rmdir(path)
					#remove the directory
					print '%s removed'%path
		else:
			self.submit()
	
	def submit(self):
		no_of_files = len(self.f_list)
		files_per_job = int( math.ceil(no_of_files/float(self.no_of_jobs)) )

		for i in range(self.no_of_jobs):
			job_fname = os.path.join(self.jobdir, 'graph_construct_'+repr(i+1)+'.sh')
			job_f = open(job_fname, 'w')
			out_block = '#!/bin/sh\n'		#this is '=' not '+='
			tmpdir = os.path.join(self.dir, 'tmp%d'%i)
			os.mkdir(tmpdir)
			for j in range(files_per_job):
				#move files to tmpdir
				index = i*files_per_job + j
				if index == no_of_files:
					break
				src = os.path.join(self.dir, self.f_list[index])
				dst = os.path.join(tmpdir, self.f_list[index])
				os.rename(src, dst)
			
			out_block += '~/script/annot/bin/batch_graph_cc.py %s %s\n'%(tmpdir, self.dstdir)
			job_f.write(out_block)
			job_f.close()
			wl = ['qsub', job_fname]
			os.spawnvp(os.P_WAIT, 'qsub', wl)


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hrn:d:g:", ["help", "restore", "no_of_jobs=", "data_dir=", "gph_dir="])
	except:
		print __doc__
		sys.exit(2)
	
	no_of_jobs = None
	restore = 0
	data_dir = None
	gph_dir = None
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-n", "--no_of_jobs"):
			no_of_jobs = int(arg)
		elif opt in ("-d", "--data_dir"):
			data_dir = arg
		elif opt in ("-g", "--gph_dir"):
			gph_dir = arg
		elif opt in ("-r", "--restore"):
			restore = 1


	if restore==1 and data_dir:
		instance = batch_submit(data_dir, gph_dir, no_of_jobs, restore)
		instance.run()
	elif no_of_jobs and data_dir and gph_dir:
		instance = batch_submit(data_dir, gph_dir, no_of_jobs, restore)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
