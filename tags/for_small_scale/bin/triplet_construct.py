#!/usr/bin/env python
import pickle,sys,os,gzip

class triplet_construct:
	def __init__(self, of):
		self.local_vertex_dict = {}
		self.local_duplet_dict = {}
		self.outf = of

	def init(self):
		self.local_vertex_dict = {}
		self.local_duplet_dict = {}
		
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
			sys.stderr.write('\t%d\n'%vertex)
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
	
	def run(self, inf):
		self.init()
		self.parse(inf)
		self.local_triplet_construct_to_file()

def triplet_batch(dir, ofname):
	files = os.listdir(dir)
	sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
	of = gzip.open(ofname, 'w')
	instance = triplet_construct(of)

	for f in files:
		pathname = os.path.join(dir, f)
		sys.stderr.write("%d/%d:\t%s\n"%(files.index(f)+1,len(files),f))
		inf = open(pathname, 'r')
		instance.run(inf)
		inf.close()
		


if __name__ == '__main__':
	def helper():
		sys.stderr.write("\
	argv[1] specifies which directory contains results from \n\
		graph construction. Those files should be gspan format.\n\
	argv[2] the file to store the resultant triplets.\n\
		The file will be gzipped.\n")
	
	if len(sys.argv) == 3:
		triplet_batch(sys.argv[1], sys.argv[2])
	else:
		helper()
		sys.exit(1)
