#!/usr/bin/env python


import os,sys

'''
Usage:	'program_name' arg1 arg2
	
	arg1 is the name of the directory containing datasets
	arg2 is the name of the directory to hold the results
'''


def batch(dataset_dir, result_dir):
	files = os.listdir(dataset_dir)
	if not os.path.isdir(result_dir):
		os.makedirs(result_dir)
	sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
	homedir = os.path.expanduser('~')
	bin_dir = os.path.join(homedir, 'script/annot/bin/')
	
	for f in files:
		sys.stderr.write("%d/%d:\t%s\n"%(files.index(f)+1,len(files),f))
		src_file = os.path.join(dataset_dir, f)
		output_file = os.path.join(result_dir, f+'.gph')
		wl = ['graph.bin', src_file, output_file, f]
		os.spawnvp(os.P_WAIT, bin_dir + 'graph.bin', wl)
		

if __name__ == '__main__':
	#argv[1] is the source directory which stores all the datasets.
	#argv[2] is the output directory which stores all the output files.
	batch(sys.argv[1],sys.argv[2])
