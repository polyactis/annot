#!/usr/bin/env python
"""
Usage: triplet_construct.py [OPTION] DATA_DIR FILE

Option:
	DATA_DIR specifies which directory contains results from
		graph construction. Those files should be gspan format.
	FILE is the file to store the resultant triplets or quadruplets.
		The file will be gzipped.
	-q, --quadruplets	search for the quadruplets, default is triplets
	-h, --help              show this help
	
Examples:
	triplet_construct.py gph_result/sc_yh60_type0 pickle/triplets.gz
	triplet_construct.py -q gph_result/sc_yh60_type0 pickle/quadruplets.gz
	
Description:
.
"""

import pickle,sys,os,gzip,getopt

class triplet_construct:
	def __init__(self, of, quadruplets=0):
		self.local_vertex_dict = {}
		self.local_duplet_dict = {}
		self.outf = of
		self.quadruplets = int(quadruplets)

	def init(self):
		self.local_vertex_dict = {}
		#sort of a linked list(see log)
		self.local_duplet_dict = {}
		#a graph dictionary
		
	def parse(self, inf):
		#loop below initilizes local_duplet_dict & local_vertex_dict.
		line = inf.readline()
		while line:
			list = line[:-1].split()
			if line[0] == 'e':
				vertex1 = int(list[1])
				vertex2 = int(list[2])
				if vertex1 < vertex2:
					#ONly add the bigger vertex to the smaller vertex's adjacent list.
					if vertex1 in self.local_vertex_dict:
						self.local_vertex_dict[vertex1].append(vertex2)
					else:
						self.local_vertex_dict[vertex1] = [vertex2]
					#two vertices is stored in ascending order in the key
					self.local_duplet_dict[(vertex1,vertex2)] = 1
				else:
					if vertex2 in self.local_vertex_dict:
						self.local_vertex_dict[vertex2].append(vertex1)
					else:
						self.local_vertex_dict[vertex2] = [vertex1]
					self.local_duplet_dict[(vertex2,vertex1)] = 1
			line = inf.readline()
			
	def local_triplet_construct_to_file(self):
		triplets_count = 0
		for vertex in self.local_vertex_dict:
			#secure that vertices in triplet and duplet are in ascending order.
			self.local_vertex_dict[vertex].sort()
			neighbour_list = self.local_vertex_dict[vertex]
			#sys.stderr.write('\t%d\n'%vertex)
			no_of_neighbours = len(neighbour_list)
			for i in range(no_of_neighbours):
				for j in range(i+1, no_of_neighbours):
					triplet = [vertex, neighbour_list[i], neighbour_list[j] ]
					duplet = [neighbour_list[i], neighbour_list[j]]
					if self.local_duplet_dict.has_key((duplet[0], duplet[1])):
						#sys.stderr.write('triplet: %s\n'%repr(triplet))
						#sys.stderr.write('duplet: %s\n'%repr(duplet))
						self.outf.write('%d,%d,%d\n'%(triplet[0],triplet[1],triplet[2],))
						triplets_count += 1
		sys.stderr.write('\tTotal triplets: %d\n'%triplets_count)

	def local_quadruplet_construct_to_file(self):
		quadruplets_count = 0
		for vertex in self.local_vertex_dict:
			#secure that vertices are in ascending order.
			self.local_vertex_dict[vertex].sort()
			neighbour_list = self.local_vertex_dict[vertex]
			#sys.stderr.write('\t%d\n'%vertex)
			no_of_neighbours = len(neighbour_list)
			for i in range(no_of_neighbours):
				for j in range(i+1, no_of_neighbours):
					if self.local_duplet_dict.has_key((neighbour_list[i], neighbour_list[j])):
						#sys.stderr.write('duplet: %s\n'%repr(duplet))
						for k in range(j+1, no_of_neighbours):
							if (neighbour_list[i], neighbour_list[k]) in self.local_duplet_dict and (neighbour_list[j], neighbour_list[k]) in self.local_duplet_dict:
								self.outf.write('%d,%d,%d,%d\n'%(vertex, neighbour_list[i], neighbour_list[j], neighbour_list[k]))
								quadruplets_count += 1
		sys.stderr.write('\tTotal quadruplets: %d\n'%quadruplets_count)

	def run(self, inf):
		self.init()
		self.parse(inf)
		if self.quadruplets==1:
			self.local_quadruplet_construct_to_file()
		else:
			self.local_triplet_construct_to_file()

def triplet_batch(dir, ofname, quadruplets):
	files = os.listdir(dir)
	sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
	of = gzip.open(ofname, 'w')
	instance = triplet_construct(of, quadruplets)

	for f in files:
		pathname = os.path.join(dir, f)
		sys.stderr.write("%d/%d:\t%s\n"%(files.index(f)+1,len(files),f))
		inf = open(pathname, 'r')
		instance.run(inf)
		inf.close()
		


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hq", ["help", "quadruplets"])
	except:
		print __doc__
		sys.exit(2)

	quadruplets = 0
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-q", "--quadruplets"):
			quadruplets = 1


	if len(args) == 2:
		triplet_batch(args[0], args[1], quadruplets)
	else:
		print __doc__
		sys.exit(2)
