#!/usr/bin/env python
import sys,os,re

class rename:
	def __init__(self):
		pass
	
	def run(self, dir):
		files = os.listdir(dir)
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
		
		for f in files:
			sys.stderr.write("%d/%d:\t%s\n"%(files.index(f)+1,len(files),f))
			self.lower_case(dir,f)
	
	def prefix_ext_swap(self, dir, f):
		src_pathname = os.path.join(dir, f)
		prefix,ext = os.path.splitext(f)
		if ext == ext[1:]:
			sys.stderr.write('the filename has no extension, ignore\n')
			return
		ext = ext[1:]
		#remove the dot(.)
		dst_pathname = os.path.join(dir, ext+'_'+prefix)
		os.rename(src_pathname, dst_pathname)
	
	def lower_case(self, dir, f):
		src_pathname = os.path.join(dir, f)
		dst_pathname = os.path.join(dir,f.lower())
		os.rename(src_pathname, dst_pathname)

if __name__ == '__main__':
	instance = rename()
	instance.run(sys.argv[1])
	#argv[1] specifies the directory containing the files to be renamed.
