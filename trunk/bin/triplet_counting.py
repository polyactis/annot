#!/usr/bin/env python
"""
Usage: triplet_counting.py [OPTION] FILE TRIPLET_DIR

Option:
	FILE	the file outputted by triplet_construct.py. '-z' specifies gzipped or not.
	TRIPLET_DIR	the directory to hold the 39 files.
	-t ..., --threshhold=...,	Memory threshhold. Specifies #triplets to be processed
		in one run. Default is 1.6e7
	-b ..., --bufferdir=...,	the directory to store the bucketing files
		at least allow space 2 times as large as FILE, default is '/tmp'
	-z, --gzipped	FILE is in gzipped(to save space) format
	-h, --help	show this help
	
Examples:
	triplet_counting.py -t 2.0e7 pickle/triplets gph_result/dir_39
	triplet_counting.py -z -b /usr/local/src pickle/quadruplets.gz gph_result/dir_39
	
Description:
.	It doesn't matter FILE is triplets or quadruplets.
"""

import pickle,sys,os,gzip,getopt


class triplet_counting:
	def __init__(self, triplet_file_path, dir, threshhold, gzipped, bufferdir):
		self.counting_dir = os.path.abspath(dir)
		if not os.path.isdir(self.counting_dir):
			os.makedirs(self.counting_dir)
		self.ancestor_fname = triplet_file_path
		self.threshhold = float(threshhold)
		self.gzipped = int(gzipped)
		if not os.path.isdir(bufferdir):
			sys.stderr.write("%s doesn't exist\n"%bufferdir)
			sys.exit(2)
		self.bufferdir = bufferdir
		self.source_fname = os.path.join(self.bufferdir, 'father')
		self.target_fname = os.path.join(self.bufferdir, 'son')
		self.triplet_file_dict = {}
		self.triplet_dict = {}
		self.no_of_runs = 0
	
	def triplet_file_dict_initialize(self):
		for i in range(1,40):
			self.triplet_file_dict[i] = open(os.path.join(self.counting_dir,'triplet_%d'%i), 'w')	

	def run(self):
		self.triplet_file_dict_initialize()
		self.atom_run(self.ancestor_fname, self.source_fname)
		if self.gzipped == 1:
			source_f = gzip.open(self.source_fname, 'r')
		else:
			source_f = open(self.source_fname, 'r')
		line = source_f.readline()
		while line:
			source_f.close()
			self.triplet_dict = {}	#ready for another 10million
			self.atom_run(self.source_fname, self.target_fname)
			
			tmp = self.source_fname
			self.source_fname = self.target_fname
			self.target_fname = tmp
			if self.gzipped == 1:
				source_f = gzip.open(self.source_fname, 'r')
			else:
				source_f = open(self.source_fname, 'r')
			line = source_f.readline()
			
		
	def atom_run(self, source, target):
		if self.gzipped:
			source_f = gzip.open(source, 'r')
		else:
			source_f = open(source, 'r')
		if self.gzipped:
			target_f = gzip.open(target, 'w')
		else:
			target_f = open(target, 'w')
		line = source_f.readline()
		while line:
			if self.triplet_dict.has_key(line):
				self.triplet_dict[line] += 1
			elif len(self.triplet_dict) <self.threshhold:
				self.triplet_dict[line] = 1
			else:
				target_f.write(line)
			line = source_f.readline()
		source_f.close()
		target_f.close()
		for i in range(len(self.triplet_dict)):
			(triplet,freq) = self.triplet_dict.popitem()
			self.triplet_file_dict[freq].write(triplet)
		self.no_of_runs += 1
		sys.stderr.write('\ttransfer %d done\n'%self.no_of_runs)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "ht:zb:", ["help", "threshhold=", "gzipped", "bufferdir="])
	except:
		print __doc__
		sys.exit(2)
	
	threshhold = 1.6e7
	gzipped = 0
	bufferdir = '/tmp'
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-t", "--threshhold"):
			threshhold = float(arg)
		elif opt in ("-z", "--gzipped"):
			gzipped = 1
		elif opt in ("-b", "--bufferdir"):
			bufferdir = arg


	if len(args) == 2:
		instance = triplet_counting(args[0], args[1], threshhold, gzipped, bufferdir)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
