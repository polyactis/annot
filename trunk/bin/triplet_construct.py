#!/usr/bin/env python
import pickle,sys,os,random

class triplet_construct:
	def __init__(self):
		self.global_triplet_dict = {}
		self.local_triplet_dict = {}
		self.local_vertex_dict = {}
		self.local_duplet_dict = {}
		
	def init(self):
		self.local_triplet_dict = {}
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
					#two vertices is stored in ascending order in the key
					self.local_duplet_dict[(vertex1,vertex2)] = 1
				else:
					self.local_duplet_dict[(vertex2,vertex1)] = 1
					
				if vertex1 in self.local_vertex_dict:
					self.local_vertex_dict[vertex1].append(vertex2)
				else:
					self.local_vertex_dict[vertex1] = [vertex2]
				if vertex2 in self.local_vertex_dict:
					self.local_vertex_dict[vertex2].append(vertex1)
				else:
					self.local_vertex_dict[vertex2] = [vertex1]
			line = inf.readline()
			
	def local_triplet_dict_construct(self):
		for vertex in self.local_vertex_dict:
			neighbour_list = self.local_vertex_dict[vertex]
			no_of_neighbours = len(neighbour_list)
			for i in range(no_of_neighbours):
				for j in range(i+1, no_of_neighbours):
					triplet = [vertex, neighbour_list[i], neighbour_list[j] ]
					triplet.sort()
					tuple = (triplet[0],triplet[1],triplet[2])	#only tuple can be hashed.
					if tuple not in self.local_triplet_dict:
						duplet = [neighbour_list[i], neighbour_list[j]]
						duplet.sort()
						if self.local_duplet_dict.has_key((duplet[0], duplet[1])):
							self.local_triplet_dict[tuple] = 1
			
	def local_to_global(self, global_dict):
		for triplet in self.local_triplet_dict:
			if global_dict.has_key(triplet):
				global_dict[triplet] += 1
			else:
				global_dict[triplet] = 1
	
	def run(self, inf):
		self.init()
		self.parse(inf)
		self.local_triplet_dict_construct()
		self.local_to_global()

def triplet_batch(dir, ofname):
	files = os.listdir(dir)
	sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
	pickle_filename = os.path.join(os.path.expanduser('~'), ofname)
	of = open(pickle_filename, 'w')
	global_dict = {}

	for f in files:
		instance = triplet_construct()
		pathname = os.path.join(dir, f)
		sys.stderr.write("%d/%d:\t%s\n"%(files.index(f)+1,len(files),f))
		inf = open(pathname, 'r')
		instance.parse(inf)
		sys.stderr.write('parsing done\n')
		instance.local_triplet_dict_construct()
		sys.stderr.write('local triplet done\n')
		instance.local_to_global(global_dict)
		sys.stderr.write('local to global transfer done\n')
		inf.close()
		del instance
		
	pickle.dump(global_dict, of)
	'''
	# this block is for triplet output.
	triplet_list = instance.global_triplet_dict.keys()
	triplet_list.sort()
	pickle_file2 = os.path.join(os.path.expanduser('~'), 'pickle/yeast_global_struc')
	global_struc = pickle.load(open(pickle_file2, 'r'))
	vertex_list = global_struc['vertex_list']
	for triplet in triplet_list:
		sys.stderr.write("%s\t%d\n"%(triplet, instance.global_triplet_dict[triplet],) )
		#sys.stderr.write("('%s', '%s', '%s')\t%d\n"%(vertex_list[triplet[0]-1], vertex_list[triplet[1]-1], vertex_list[triplet[2]-1],instance.global_triplet_dict[triplet],))
	'''


if __name__ == '__main__':
	triplet_batch(sys.argv[1], sys.argv[2])
	# argv[1] specifies which directory contains results from graph construction
	# the resultant triplet dictionary will be stored in '~/pickle/yeast_triplet'.
