#!/usr/bin/env python
import pickle,sys,os,gzip


class triplet_counting:
	def __init__(self, triplet_file_path, dir):
		self.counting_dir = os.path.abspath(dir)
		if not os.path.isdir(self.counting_dir):
			os.makedirs(self.counting_dir)
		self.ancestor_fname = triplet_file_path
		self.source_fname = '/tmp/father'
		self.target_fname = '/tmp/son'
		self.triplet_file_dict = {}
		self.triplet_dict = {}
		self.no_of_runs = 0
	def triplet_file_dict_initialize(self):
		for i in range(1,40):
			self.triplet_file_dict[i] = open(os.path.join(self.counting_dir,'triplet_%d'%i), 'w')	

	def run(self):
		self.triplet_file_dict_initialize()
		self.atom_run(self.ancestor_fname, self.source_fname)
		source_f = open(self.source_fname, 'r')
		line = source_f.readline()
		while line:
			source_f.close()
			self.triplet_dict = {}	#ready for another 10million
			self.atom_run(self.source_fname, self.target_fname)
			
			tmp = self.source_fname
			self.source_fname = self.target_fname
			self.target_fname = tmp
			source_f = open(self.source_fname, 'r')
			line = source_f.readline()
			
		
	def atom_run(self, source, target):
		source_f = open(source, 'r')
		target_f = open(target, 'w')
		line = source_f.readline()
		while line:
			if line in self.triplet_dict:
				self.triplet_dict[line] += 1
			elif len(self.triplet_dict) <2.3e7:
				self.triplet_dict[line] = 1
			else:
				target_f.write(line)
			line = source_f.readline()
		source_f.close()
		target_f.close()
		for item in self.triplet_dict.iteritems():
			self.triplet_file_dict[item[1]].write(item[0])
		self.no_of_runs += 1
		sys.stderr.write('\ttransfer %d done\n'%self.no_of_runs)

if __name__ == '__main__':
	def helper():
		sys.stderr.write('\
	argv[1] is the gzipped triplet file(from triplet_construct.py).\n\
	argv[2] is the directory to hold the 39 files.\n')
	
	if len(sys.argv) == 3:
		instance = triplet_counting(sys.argv[1], sys.argv[2])
	else:
		helper()
		sys.exit(1)
	instance.run()
