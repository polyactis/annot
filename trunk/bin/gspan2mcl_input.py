#!/usr/bin/env python
"""
Usage: gspan2mcl_input.py [OPTION] INPUTFILE OUTPUTFILE

Option:
	INPUTFILE is in gspan format.
	OUTPUTFILE will be in mcl_input format.
	-h, --help              show this help
	
Examples:
	gspan2mcl_input.py gph_result/sc/gph_dataset1 gph_result/sc_mcl/mcl_gph_dataset1

Description:
	This program converts gspan formatted file generated by graph_reorganize.py
	to file in mcl_input format.
	Use batch.qsub.py to iterate over files in a directory.
"""


import sys, os, re, getopt

class gspan2mcl_input:
	'''
	'''
	def __init__(self, infname, ofname):
		self.infname = infname
		self.inf = open(self.infname, 'r')
		self.of = open(ofname, 'w')
		self.vertex_dict = {}
	
	def run(self):
		p_no = re.compile(r'\d+$')
		try:
			no = int(p_no.search(self.infname).group())	#integer conversion
		except AttributeError, error:
			sys.stderr.write('%s\n'%error)
			sys.stderr.write("can't find the number of this dataset from %s\n"%self.infname)
			sys.exit(2)
		for line in self.inf:
			if line[0] == 'e':
				#edge here, like 'e 3807 3859 0.804645'
				line_list = line[:-1].split()
				vertex1 = line_list[1]
				vertex2 = line_list[2]
				if vertex1 in self.vertex_dict:
					self.vertex_dict[vertex1].append(vertex2)
				else:
					self.vertex_dict[vertex1] = [vertex2]
				if vertex2 in self.vertex_dict:
					self.vertex_dict[vertex2].append(vertex1)
				else:
					self.vertex_dict[vertex2] = [vertex1]
		dim = len(self.vertex_dict)
		out_block = '(splat_id %s )\n'%no		# here it is '=' not '+='
		out_block += '(mclheader\n'
		out_block += 'mcltype matrix\n'
		out_block += 'dimensions %dx%d\n)\n'%(dim,dim)
		out_block += '(mcldoms\n'
		vertex_list = self.vertex_dict.keys()
		vertex_list.sort()
		out_block += '%s $\n)\n'%' '.join(vertex_list)
		out_block += '(mclmatrix\nbegin\n'
		self.of.write(out_block)
		for vertex in vertex_list:
			self.of.write('%s '%vertex)
			self.of.write('%s $\n'%' '.join(self.vertex_dict[vertex]))
		self.of.write(')\n\n')
			
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "h", ["help"])
	except:
		print __doc__
		sys.exit(2)
	
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)

			
	if len(args) == 2:
		instance = gspan2mcl_input(args[0], args[1])
		instance.run()
	else:
		print __doc__
		sys.exit(2)