#!/usr/bin/env python
import pickle,sys,os,random


class triplet_stat:
	def __init__(self):
		self.recurrence_triplet_list = []
		self.recurrence_stat_func_list = []
		self.recurrence_stat_trans_list = []
		self.triplet_dict = {}
		self.triplet_pickle_fname = os.path.join(os.path.expanduser('~'),'pickle/yeast_triplet')
		self.transfac_pickle_fname = os.path.join(os.path.expanduser('~'),'pickle/yeast_transfac_dict')
		self.known_genes_pickle_fname = os.path.join(os.path.expanduser('~'),'pickle/known_genes_dict')
		self.recurrence_stat_fname = os.path.join(os.path.expanduser('~'),'pickle/recurrence_stat')
		self.global_struc_fname = os.path.join(os.path.expanduser('~'), 'pickle/yeast_global_struc')
		self.vertex_list = []
		self.transfac_dict = {}
		self.known_genes_dict = {}
	
	def dstruc_loadin(self):
		if not os.path.isfile(self.transfac_pickle_fname):
			sys.stderr.write('You need to construct the transfac_dict first.\n')
			sys.exit(1)
		if not os.path.isfile(self.known_genes_pickle_fname):
			sys.stderr.write('You need to construct known_genes_dict first.\n')
			sys.exit(1)
		if not os.path.isfile(self.triplet_pickle_fname):
			sys.stderr.write('You need to construct triplet dict first.\n')
			sys.exit(1)
		if not os.path.isfile(self.global_struc_fname):
			sys.stderr.write('You need to construct global_struc(vertex labeling) first.\n')
			sys.exit(1)
		self.transfac_dict = pickle.load(open(self.transfac_pickle_fname,'r'))
		self.known_genes_dict = pickle.load(open(self.known_genes_pickle_fname,'r'))		
		self.triplet_dict = pickle.load(open(self.triplet_pickle_fname, 'r'))
		global_struc = pickle.load(open(self.global_struc_fname, 'r'))
		self.vertex_list = global_struc['vertex_list']
		
	def transfac_dict_construct(self, inf):
		line = inf.readline()
		while line:
			list = line[:-1].split()
			if self.transfac_dict.has_key(list[0]):
				self.transfac_dict[list[0]].append(list[1])
			else:
				self.transfac_dict[list[0]] = [list[1]]
			line = inf.readline()
		pickle.dump(self.transfac_dict, open(self.transfac_pickle_fname, 'w'))
		
	def recurrence_triplet_list_construct(self):
		for i in range(10):
			self.recurrence_triplet_list.append([])
		while self.triplet_dict:
			item = self.triplet_dict.popitem()
			vertex = item[0]
			freq = item[1]
			if freq>10 or len(self.recurrence_triplet_list[freq])==100000:
				continue
			self.recurrence_triplet_list[freq-1].append(vertex)
				
	def recurrence_stat_list_construct(self):
		for i in xrange(10):
			self.recurrence_stat_func_list.append([])
			self.recurrence_stat_trans_list.append([])
			for j in xrange(100):
				functional_homo = 0
				transcriptional_homo = 0
				len = len(self.recurrence_triplet_list[i])
				if len < 100000:
					sys.stderr.write('\tonly %d triplets with recurrence %d\n'%(len, i+1,))
				if len <1000:
					sys.stderr.write('\ttriplets with recurrence %d is below 1000(%d). Ignore.\n'%( i+1,len,))
					continue
				index_list = random.sample(xrange(len),1000)
				for k in xrange(1000):
					triplet = self.recurrence_triplet_list[i][index_list[k]]
					if is_homogenious(triplet, self.known_genes_dict):
						functional_homo += 1
					if is_homogenious(triplet, self.transfac_dict):
						transcriptional_homo += 1
				functional_homo_ratio = functional_homo/1000.00
				transcriptional_homo_ratio = transcriptional_homo/1000.00
				self.recurrence_stat_func_list[i].append(functional_homo_ratio)
				self.recurrence_stat_trans_list[i].append(transcriptional_homo_ratio)
		
		object['func_list'] = self.recurrence_stat_func_list
		object['trans_list'] = self.recurrence_stat_trans_list
		pickle.dump(object, open(self.recurrence_stat_fname,'w'))
		
	def is_homogenious(self, triplet, dict):
		vertex1_label = self.vertex_list[triplet[0]-1]
		vertex2_label = self.vertex_list[triplet[1]-1]
		vertex3_label = self.vertex_list[triplet[2]-1]
		list = dict[vertex1_label] + dict[vertex2_label] + dict[vertex3_label]
		judge_dict = {}
		for item in list:
			if judge_dict.has_key(item):
				judge_dict[item] += 1
			else:
				judge_dict[item] = 1
			if judge_dict[item] == 3:
				return 1
		return 0
		
		
if __name__ == '__main__':
	instance = triplet_stat()
	instance.dstruc_loadin()
	instance.recurrence_triplet_list_construct()
	instance.recurrence_stat_list_construct()
	'''
	# this block is for transfac_dict construction.
	inf = open(sys.argv[1])
	instance.transfac_dict_construct(inf)
	key_list = instance.transfac_dict.keys()
	key_list.sort()
	for item in key_list:
		print '%s\t%s'%(item,instance.transfac_dict[item],)
	'''
