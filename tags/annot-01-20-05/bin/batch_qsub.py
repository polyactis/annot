#!/usr/bin/env python
"""
Usage: batch_qsub.py -c PROGRAM_TO_BE_INVOKED -i INIT_NUM -n NO_OF_JOBS -o OUTPUT_PREFIX [OPTIONS] DATADIR

Option:
	DATADIR is the directory containing all the files to be 'qsub'ed or 'sh'ed.
	-c ..., --call=... Call the program. FULL PATH please.
	-j ..., --job_prefix=...	the prefix of all the jobs, the name of the program called(default)
	-i ..., --initnum=...	The initial number for all the jobs, 0 (default)
	-n ..., --no_of_jobs=...	specify how many jobs you want to submit
	-p ..., --parameter=...	parameters to be passed to the program invoked, DOUBLE QUOTE IT
	-o ..., --oprefix=...	output prefix, FULL PATH please, 'out'(default)
	-a ..., --argument_order=...	3 letter combination of 'p', 'i', 'o', '_'; 'pio'(default)
	-s, --shell	Execute the job file via call 'sh', not 'qsub'(default)
	-f, --file_form		the input for PROGRAM_TO_BE_INVOKED is a file, not a directory(default)
	-h, --help              show this help
	
Examples:
	
	batch_qsub.py -c ~/script/annot/bin/triplet_construct.py -n 10 -p "-q -z" gph_result/sc_yh60_type0
		:submit 10 quadruplet construction jobs
	batch_qsub.py -c ~/script/annot/bin/batch_graphcc.py -j graph_construct -n 10 -o gph_sc_dataset  datasets/yeast_data
		:submit 10 graph construction jobs,
	batch_qsub.py -c ~/script/annot/bin/gspan2splat_output.py -n 38 -p "-n 38"
		-f -s -o gph_result/sc_yh60_fp_splat/splat gph_result/sc_yh60_fp_type1
		:batch run gspan2splat_output over 38 files.
	batch_qsub.py -c ~/script/annot/bin/fim_closed -n 10 -p "6" -f -o fim_input/out -a 'ipo' fim_input/
		:batch submit 10 fim_closed jobs

Description:
	This program partitions a big job into several subjobs and submits them via qsub
	or runs them via 'sh'. It can either move all the files into several sub directories
	and call the program to run over each directory, or run the program over the files.
	About the argument_order:
		'p' is parameter.
		'i' stands for input.
		'o' is ouput.
		'_' is null.
	It must be three-letter long. If can't, use '_' instead.

Warning:
	If you want submit the same job again. Please remove the tmp directories first.
	
"""


import sys, os, math, getopt

class batch_qsub:
	'''
	Partition a big job into several small jobs
	'''
	def __init__(self, dir, call, job_prefix, initnum, no_of_jobs, parameter, oprefix, argument_order, shell, file_form):
		self.dir = os.path.abspath(dir)
		self.program = os.path.abspath(call)
		if job_prefix == '':
			job_prefix = os.path.basename(call)
		self.job_prefix = job_prefix
		self.initnum = int(initnum)
		self.no_of_jobs = int(no_of_jobs)
		self.parameter = parameter
		self.oprefix = oprefix
		self.argument_order = argument_order
		#'_' must be initialized first. others (p, i, o) later in the program
		self.argument_dict = {'_':''}
		self.submit_command = 'qsub'
		if shell == 1:
			self.submit_command = 'sh'
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
			job_fname = os.path.join(self.jobdir, '%s%d.sh'%(self.job_prefix, self.initnum+i))
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
					self.argument_dict['p'] = self.parameter
					self.argument_dict['i'] = src
					self.argument_dict['o'] = '%s_%s'%(self.oprefix, self.f_list[index])
					command_block = "%s %s %s %s"%(self.program, self.argument_dict[self.argument_order[0]],\
						self.argument_dict[self.argument_order[1]], self.argument_dict[self.argument_order[2]])
					out_block += 'echo "%s"\n'%command_block
					out_block += command_block
					out_block += '\n'
				else:
					#move files to tmpdir
					dst = os.path.join(tmpdir, self.f_list[index])
					os.symlink(src, dst)
			if not self.file_form:
			#run the program over the directory
				self.argument_dict['p'] = self.parameter
				self.argument_dict['i'] = tmpdir
				self.argument_dict['o'] = '%s%d'%(self.oprefix, i)
				command_block = "%s %s %s %s"%(self.program, self.argument_dict[self.argument_order[0]],\
					self.argument_dict[self.argument_order[1]], self.argument_dict[self.argument_order[2]])
				out_block += 'echo "%s"\n'%command_block
				out_block += command_block
				out_block += '\n'
			out_block += 'date\n'
			job_f = open(job_fname, 'w')
			job_f.write(out_block)
			job_f.close()
			wl = [self.submit_command, job_fname]
			sys.stderr.write('\tJob %d: %s\n'%(i, job_fname))
			os.spawnvp(os.P_WAIT, self.submit_command, wl)
			if index == no_of_files:
			#there're no files left to allocate, so break.
			#otherwise, you will submit 'empty' jobs
				break


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hc:j:i:n:p:o:a:sf", \
			["help", "call=", "job_prefix=", "initnum=", "no_of_jobs=", "parameter=", "oprefix=", "argument_order=", "shell", "file_form"])
	except:
		print __doc__
		sys.exit(2)
	
	call = None
	job_prefix = ''
	initnum = 0
	no_of_jobs = None
	parameter = ''
	oprefix = 'out'
	argument_order = 'pio'
	shell = 0
	file_form = 0
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-c", "--call"):
			call = arg
		elif opt in ("-j", "--job_prefix"):
			job_prefix = arg
		elif opt in ("-i", "--initnum"):
			initnum = arg
			#will be integered in the class
		elif opt in ("-n", "--no_of_jobs"):
			no_of_jobs = int(arg)
		elif opt in ("-p", "--parameter"):
			parameter = arg
		elif opt in ("-o", "--oprefix"):
			oprefix = arg
		elif opt in ("-a", "--argument_order"):
			argument_order = arg
		elif opt in ("-s", "--shell"):
			shell = 1 
		elif opt in ("-f", "--file_form"):
			file_form = 1


	if call and no_of_jobs and len(args) == 1:
		instance = batch_qsub(args[0], call, job_prefix, initnum, no_of_jobs, parameter, oprefix, argument_order, shell, file_form)
		instance.submit()
	else:
		print __doc__
		sys.exit(2)
