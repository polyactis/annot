#!/usr/bin/python

import os,sys
from graph import *

'''
Usage:	'program_name' arg1 arg2

	arg1 is the name of the directory containing the datasets
	arg2 is the name of the file to hold the results of ALL the datasets.
		only the filename. it'll be put into user's home directory.
		arg2 is different from the arg2 of batch_graph_cc.py.
'''


def batch(dir, of_name):
	files = os.listdir(dir)
	sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
	of_pathname = os.path.join(os.path.expanduser('~'), of_name)
	of = open(of_pathname, 'w')
	
	for f in files:
		pathname = os.path.join(dir, f)
		sys.stderr.write("%d/%d:\t%s\n"%(files.index(f)+1,len(files),f))
		instance = graph_construct(pathname)
		instance.mask_array_construct()
		instance.edge_construct()
		instance.cleanup()
		instance.output(out=of)
		del instance
		
if __name__ == '__main__':
	batch(sys.argv[1],sys.argv[2])
	# argv[1] specifies which directory contains data
	# argv[2] specifies the outputfile name in the user's home directory
