#!/usr/bin/python

import os,sys
from graph import *

def batch(dir):
	files = os.listdir(dir)
	sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
	of_pathname = os.path.join(os.getcwd(), 'output.gph')
	of = open(of_pathname, 'w')
	
	for f in files:
		pathname = os.path.join(dir, f)
		sys.stderr.write("%d/%d:\t%s\n"%(files.index(f)+1,len(files),f))
		instance = graph_construct(pathname)
		instance.mask_array_construct()
		instance.edge_construct()
		instance.cleanup()
		instance.output(out=of)
		
if __name__ == '__main__':
	batch(sys.argv[1])
