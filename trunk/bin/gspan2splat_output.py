#!/usr/bin/env python
"""
Usage: gspan2splat_output.py -n NO_OF_DATASETS [OPTION] INPUTFILE OUTPUTFILE

Option:
	INPUTFILE is in gspan format.
	OUTPUTFILE will be in splat output format.
	-n ..., --no_of_datasets=...	specify the total number of datasets
	-h, --help              show this help
	
Examples:
	gspan2splat_output.py -n 38 gph_result/sc/gph_dataset1 gph_result/sc_splat/splat_gph_dataset1

Description:
	This program converts gspan formatted file generated by graph_reorganize.py
	to file in splat output format.
	Use batch.qsub.py to iterate over files in a directory.
"""


import sys, os, re, getopt

class gspan2splat_output:
	'''
	'''
	def __init__(self, infname, ofname, no_of_datasets):
		self.infname = infname
		self.inf = open(self.infname, 'r')
		self.of = open(ofname, 'w')
		#used to construct the recurrence_pattern bit string
		self.no_of_datasets = int(no_of_datasets)
		self.no_of_edges = 0
		self.recurrence_pattern = ''
		self.edge_string = ''
	
	def run(self):
		#construct the recurrence_pattern from the filename
		p_no = re.compile(r'\d+$')
		try:
			no = int(p_no.search(self.infname).group())	#integer conversion
		except AttributeError, error:
			sys.stderr.write('%s\n'%error)
			sys.stderr.write("can't construct recurrence_pattern from %s\n"%self.infname)
			sys.exit(2)
		recurrence_array = ['0']*self.no_of_datasets
		recurrence_array[no-1] = '1'
		self.recurrence_pattern = ''.join(recurrence_array)
		#I don't care about no_of_edges, 0.
		self.of.write('%d\n%s\n'%(self.no_of_edges, self.recurrence_pattern))
		for line in self.inf:
			if line[0] == 'e':
				#edge here, like 'e 3807 3859 0.804645'
				line_list = line[:-1].split()
				self.no_of_edges += 1
				self.of.write('(%s %s)'%(line_list[1], line_list[2]))
				if (self.no_of_edges%6) == 0:
					#every six edges, one line
					self.of.write('\n')
		if (self.no_of_edges%6) != 0:
			#add '\n' to the last edge line
			self.of.write('\n')
		#ensure the file is ended with '\n'
		self.of.write('\n')

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hn:", ["help", "no_of_datasets="])
	except:
		print __doc__
		sys.exit(2)
	
	no_of_datasets = None
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-n", "--no_of_datasets"):
			no_of_datasets = int(arg)

			
	if no_of_datasets and len(args) == 2:
		instance = gspan2splat_output(args[0], args[1], no_of_datasets)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
