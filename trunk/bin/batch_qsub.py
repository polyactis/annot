#!/usr/bin/env python
"""
Usage: batch_qsub.py -c PROGRAM_TO_BE_INVOKED -i INIT_NUM -n NO_OF_JOBS -o OUTPUT_PREFIX [OPTIONS] DATADIR

Option:
	DATADIR is the directory containing all the files to be qsubbed.
	-c ..., --call=... Call the program. FULL PATH please.
	-i ..., --initnum=...	The initial number for all the jobs.
	-n ..., --no_of_jobs=...	specify how many jobs you want to submit
	-p ..., --parameter=...	parameters to be passed to the program invoked, DOUBLE QUOTE IT
	-o ..., --oprefix=...	output prefix, FULL PATH please.
	-f, --file_form		the input for PROGRAM_TO_BE_INVOKED is a file, not a directory(default)
	-h, --help              show this help
	
Examples:
	
	batch_qsub.py -c ~/script/annot/bin/triplet_construct.py -i 0 -n 10 -p "-q -z" -o sc_quadruplets gph_result/sc_yh60_type0
		:submit 10 quadruplet construction jobs
	batch_qsub.py -c ~/script/annot/bin/graph.bin -i 0 -n 10 -o gph_sc_dataset -f datasets/yeast_data
		:submit 10 graph construction jobs, in file form. This is a bad example. Don't use it.
	
Description:
	This program partitions a big job into several subjobs and submits them via qsub.
	It can either move all the files into several sub directories and call the program
	to run over each directory, or run the program over the files.

Warning:
	If you want submit the same job again. Please remove the tmp directories first.
	
"""


import sys, os, math, getopt

class batch_qsub:
	'''
	Partition a big job into several small jobs
	'''
	def __init__(self, dir, call, initnum, no_of_jobs, parameter, oprefix, file_form):
		self.dir = os.path.abspath(dir)
		self.program = os.path.abspath(call)
		self.initnum = int(initnum)
		self.no_of_jobs = int(no_of_jobs)
		self.parameter = parameter
		self.oprefix = oprefix
		self.file_form = int(file_form)
		
		if not os.path.isdir(self.dir):
			sys.stderr.write('%s not present\n'%self.dir)
			sys.exit(1)
		
		self.jobdir = os.path.join(os.path.expanduser('~'),'qjob')
		if not os.path.isdir(self.jobdir):
			os.makedirs(self.jobdir)
		self.f_list = os.listdir(self.dir)
	
	
	def submit(self):
		no_of_files = len(self.f_list)
		files_per_job = int( math.ceil(no_of_files/float(self.no_of_jobs)) )

		for i in range(self.no_of_jobs):
			job_fname = os.path.join(self.jobdir, 'batch_qsub_'+repr(self.initnum+i)+'.sh')
			out_block = '#!/bin/sh\n'		#this is '=' not '+='
			out_block += 'date\n'
			if not self.file_form:
				tmpdir = os.path.join(self.dir, 'tmp%d'%i)
				if os.path.isdir(tmpdir):
					sys.stderr.write("Please remove tmp directories first.\n")
					sys.exit(2)
				os.mkdir(tmpdir)
			for j in range(files_per_job):
				index = i*files_per_job + j
				if index >= no_of_files:
					break
				src = os.path.join(self.dir, self.f_list[index])
				if self.file_form:
				#run the program on each file
					command_block = "%s %s %s %s_%s"%(self.program, self.parameter, src, self.oprefix, self.f_list[index])
					out_block += 'echo "%s"\n'%command_block
					out_block += command_block
					out_block += '\n'
				else:
					#move files to tmpdir
					dst = os.path.join(tmpdir, self.f_list[index])
					os.symlink(src, dst)
			if not self.file_form:
			#run the program over the directory
				command_block = '%s %s %s %s%d'%(self.program, self.parameter, tmpdir, self.oprefix, i)
				out_block += 'echo "%s"\n'%command_block
				out_block += command_block
				out_block += '\n'
			out_block += 'date\n'
			job_f = open(job_fname, 'w')
			job_f.write(out_block)
			job_f.close()
			wl = ['qsub', job_fname]
			os.spawnvp(os.P_WAIT, 'qsub', wl)
			if index == no_of_files:
			#there're no files left to allocate, so break.
			#otherwise, you will submit 'empty' jobs
				break


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hc:i:n:p:o:f", ["help", "call=", "initnum=", "no_of_jobs=", "parameter=", "oprefix=", "file_form"])
	except:
		print __doc__
		sys.exit(2)
	
	call = None
	initnum = None
	no_of_jobs = None
	parameter = ''
	oprefix = ''
	file_form = 0
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-c", "--call"):
			call = arg
		elif opt in ("-i", "--initnum"):
			initnum = arg
			#will be integered in the class
		elif opt in ("-n", "--no_of_jobs"):
			no_of_jobs = int(arg)
		elif opt in ("-p", "--parameter"):
			parameter = arg
		elif opt in ("-o", "--oprefix"):
			oprefix = arg
		elif opt in ("-f", "--file_form"):
			file_form = 1


	if call and initnum and no_of_jobs and oprefix and len(args) == 1:
		instance = batch_qsub(args[0], call, initnum, no_of_jobs, parameter, oprefix, file_form)
		instance.submit()
	else:
		print __doc__
		sys.exit(2)
