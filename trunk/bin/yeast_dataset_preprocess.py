#!/usr/bin/env python
import pickle,sys,os

class yeast_two_way_filter:
	'''
	This class is different from the one in known_genes_filter.py. That one preprocess the yeast graph_construction results from graph.cc.
	It is specific for kangyu's project.
	Preprocess yeast raw dataset files.
	first remove the single quotes around the orfname(normal_outf).
	secondly, tag the transcription factor orfname(trans_factor_tagged_outf)
	'''
	def __init__(self):
		self.trans_factor_dict = {}
	
	def trans_factor_dict_construct(self, inf):
		line = inf.readline()
		while line:
			orf = line[:-1]
			if orf not in self.trans_factor_dict:
				self.trans_factor_dict[orf] = 1
			line = inf.readline()

	def filter_based_on_trans_factor_dict(self, inf, outf1, outf2):
		self.normal_outf = outf1
		self.trans_factor_tagged_outf = outf2
		line = inf.readline()
		while line:
			list = line[:-1].split('\t')
			if list[0][0]==list[0][-1] and list[0][0] == "'":
				list[0] = list[0][1:-1]
			vertex = list[0]
			self.normal_outf.write('%s\n'%'\t'.join(list))
			if vertex in self.trans_factor_dict:
				list[0] = 't' + list[0]
				self.trans_factor_tagged_outf.write('%s\n'%'\t'.join(list))
			else:
				self.trans_factor_tagged_outf.write('%s\n'%'\t'.join(list))
			line = inf.readline()


def two_way_filter_batch(dir, output_dir, trans_fname):
	files = os.listdir(dir)
	sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
	if dir[-1] == '/':
		dir = dir[:-1]
		#the last '/' will affect the behavior of os.path.basename.
	basedir = os.path.basename(dir)
	normal_dir = os.path.join(output_dir,basedir+'_normal')
	if not os.path.isdir(normal_dir):
		os.makedirs(normal_dir)
	trans_factor_tagged_dir = os.path.join(output_dir,basedir+'_trans_factor_tagged')
	if not os.path.isdir(trans_factor_tagged_dir):
		os.makedirs(trans_factor_tagged_dir)
	#pickle_file = os.path.join(os.path.expanduser('~'), 'pickle/known_genes_dict')
	#known_genes_dict = pickle.load(open(pickle_file,'r'))
	instance = yeast_two_way_filter()
	trans_f = open(trans_fname,'r')
	instance.trans_factor_dict_construct(trans_f)
	
	for f in files:
		pathname = os.path.join(dir, f)
		sys.stderr.write("%d/%d:\t%s\n"%(files.index(f)+1,len(files),f))
		inf = open(pathname, 'r')
		normal_f = open(os.path.join(normal_dir,f), 'w')
		trans_factor_tagged_f = open(os.path.join(trans_factor_tagged_dir, f), 'w')
		instance.filter_based_on_trans_factor_dict(inf, normal_f, trans_factor_tagged_f)
		inf.close()
		normal_f.close()
		trans_factor_tagged_f.close()


if __name__ == '__main__':
	two_way_filter_batch(sys.argv[1],sys.argv[2], sys.argv[3])
	# argv[1] specifies which directory contains raw yeast data
	# argv[2] specifies the directory to store two directories of datafiles

