#!/usr/bin/env python
import pickle,sys,os

class triplet_construct:
	def __init__(self):
		self.global_tirplet_dict = {}
		self.global_recurrence_list = []
		self.local_triplet_dict = {}
		self.local_vertex_dict = {}
		self.local_duplet_dict = {}
		
	def init(self):
		self.local_triplet_dict = {}
		self.local_vertex_dict = {}
		self.local_duplet_dict = {}
		
	def parse(self, inf):
		#loop below initilizes graph_list_dict.
		line = inf.readline()
		while line:
			list = line[:-1].split('\t')
			if line[0] == 'e':
				vertex1 = list[1]
				vertex2 = list[2]
				self.local_duplet_dict[(vertex1,vertex2)] = 1
				if vertex1 in self.local_vertex_dict:
					self.local_vertex_dict[vertex1].append(vertex2)
				else:
					self.local_vertex_dict[vertex1] = []
				if vertex2 in self.local_vertex_dict:
					self.local_vertex_dict[vertex2].append(vertex1)
				else:
					self.local_vertex_dict[vertex2] = []
			line = inf.readline()
			
	def local_triplet_dict_construct(self):
		for vertex in self.local_vertex_dict:
			neighbour_list = self.local_vertex_dict[vertex]
			no_of_neighbours = len(neighbour_list)
			for i in range(no_of_neighbours):
				for j in range(i+1, no_of_neighbours):
					triplet = [vertex, neighbour_list[i], neighbour_list[j] ]
					triplet.sort()
					tuple = (triplet[0],triplet[1],triplet[2])
					if tuple not in self.local_triplet_dict:
						if self.local_duplet_dict.has_key((neighbour_list[i],neighbour_list[j])):
							self.local_triplet_dict[tuple] = 1
			
	def local_to_global(self):
		for triplet in self.local_triplet_dict:
			if self.global_tirplet_dict.has_key(triplet):
				self.global_tirplet_dict[triplet] += 1
			else:
				self.global_tirplet_dict[triplet] = 1
	
	def run(self, inf):
		self.init()
		self.parse(inf)
		self.local_triplet_dict_construct()
		self.local_to_global()

def batch(dir, of_name):
	files = os.listdir(dir)
	sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
	of_pathname = os.path.join(os.path.expanduser('~'), of_name)
	of = open(of_pathname, 'w')
	instance = triplet_construct()
	
	for f in files:
		pathname = os.path.join(dir, f)
		sys.stderr.write("%d/%d:\t%s\n"%(files.index(f)+1,len(files),f))
		inf = open(pathname, 'r')
		instance.run(inf)
		inf.close()
		
	pickle.dump ( instance.global_tirplet_dict, of )
	triplet_list = instance.global_tirplet_dict.keys()
	triplet_list.sort()
	for triplet in triplet_list:
		sys.stderr.write('%s\t%d\n'%(triplet, instance.global_tirplet_dict[triplet],))
	
if __name__ == '__main__':
	batch(sys.argv[1],sys.argv[2])
	# argv[1] specifies which directory contains results from graph construction
	# argv[2] specifies the outputfile name in the user's home directory
