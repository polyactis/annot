#!/usr/bin/env python
import sys,os,cStringIO


class mcl_input_iterator:
	'''looping over a mcl result file, generate a block of clusters related to one splat pattern'''
	def __init__(self, inf):
		self.inf = inf
		self.mcl_input_block = ''
	def __iter__(self):
		return self		
	def next(self):
		self.read()
		return cStringIO.StringIO(self.mcl_input_block)
	def read(self):
		self.mcl_input_block = ''
		mcl_begin = 0
		line = self.inf.readline()
		while line != '\n':
			if line == '':
				raise StopIteration
				break
			if mcl_begin == 1:
				self.mcl_input_block += line
			if line.find('(splat_id') == 0:
				mcl_begin = 1
				self.mcl_input_block += line
			line = self.inf.readline()
		self.mcl_input_block += line

class mcl_batch_run:
	'''
	read a concatenated mcl input file. Iteratively run mcl on the mcl_input_block.
	And concatenate the mcl output into one big file.
	'''
	def __init__(self, infname, outfname, parameter=''):
		self.inf = open(infname,'r')
		self.outf = open(outfname, 'w')
		self.parameter = parameter
		basename_of_infname = os.path.basename(infname)
		self.mclinput_fname = '/tmp/%s_in'%basename_of_infname
		self.mcloutput_fname = '/tmp/%s_out'%basename_of_infname
	
	def run(self):
		iter = mcl_input_iterator(self.inf)
		for mcl_input_block in iter:
			out_block = mcl_input_block.readline()
			mclinf = open(self.mclinput_fname, 'w')
			mclinf.write(mcl_input_block.read())
			mclinf.close()
			parameter = self.parameter.split()
			wl = ['mcl', self.mclinput_fname, '-o',self.mcloutput_fname] + parameter
			try:
				os.spawnvp(os.P_WAIT, 'mcl', wl)
			except:
				sys.stderr.write('MCL running error.\n')
				os.remove(self.mclinput_fname)
				if os.isfile(self.mcloutput_fname):
					os.remove(self.mcloutput_fname)
				sys.exit(1)
			out_block += '(parameter %s )\n'%self.parameter
			mcloutf = open(self.mcloutput_fname, 'r')
			out_block += mcloutf.read()
			mcloutf.close()
			self.outf.write(out_block)
			
		os.remove(self.mclinput_fname)
		os.remove(self.mcloutput_fname)
		
if __name__ == '__main__':
	def helper():
		sys.stderr.write('\
	argv[1] is the input file(concatenated mcl input).\n\
	argv[2] is the output file(concatenated mcl output).\n\
	argv[3] is the parameter to be passed to mcl. Ignore to use default.\n')
		
	if len(sys.argv) == 4:
		instance = mcl_batch_run(sys.argv[1], sys.argv[2], sys.argv[3])
	elif len(sys.argv) == 3:
		instance = mcl_batch_run(sys.argv[1], sys.argv[2])
	else:
		helper()
		sys.exit(1)
	instance.run()
