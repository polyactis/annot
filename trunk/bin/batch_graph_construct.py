#!/usr/bin/python

import os,sys
from graph import *

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
