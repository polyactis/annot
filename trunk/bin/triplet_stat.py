#!/usr/bin/env python
import pickle,sys,os,random


class triplet_stat:
	def __init__(self):
		self.recurrence_triplet_list = []
		self.recurrence_stat_func_list = []
		self.recurrence_stat_trans_list = []
		self.transfac_pickle_fname = os.path.join(os.path.expanduser('~'),'pickle/yeast_transfac_dict')
		self.known_genes_pickle_fname = os.path.join(os.path.expanduser('~'),'pickle/known_genes_dict')
		self.global_struc_fname = os.path.join(os.path.expanduser('~'), 'pickle/yeast_global_struc')
		self.vertex_dict = {}
		self.transfac_dict = {}
		self.known_genes_dict = {}
	
	def dstruc_loadin(self):
		if not os.path.isfile(self.transfac_pickle_fname):
			sys.stderr.write('You need to construct the transfac_dict first.\n')
			sys.exit(1)
		if not os.path.isfile(self.known_genes_pickle_fname):
			sys.stderr.write('You need to construct known_genes_dict first.\n')
			sys.exit(1)
		if not os.path.isfile(self.global_struc_fname):
			sys.stderr.write('You need to construct global_struc(vertex labeling) first.\n')
			sys.exit(1)
		global_struc = pickle.load(open(self.global_struc_fname, 'r'))
		self.vertex_dict = global_struc['vertex_dict']

		transfac_dict_orf = pickle.load(open(self.transfac_pickle_fname,'r'))

		#replace the orfname with no. because the triplets are in no. form.
		for item in transfac_dict_orf.iteritems():
			if self.vertex_dict.has_key(item[0]):
				no = self.vertex_dict[item[0]]
				self.transfac_dict[no] = item[1]
			
		known_genes_dict_orf = pickle.load(open(self.known_genes_pickle_fname,'r'))		
		for item in known_genes_dict_orf.iteritems():
			if self.vertex_dict.has_key(item[0]):
				no = self.vertex_dict[item[0]]
				self.known_genes_dict[no] = item[1]

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
		
	def recurrence_triplet_list_construct(self, triplet_fname):
		inf = open(triplet_fname, 'r')
		line = inf.readline()
		while line:
			vertex_list = line[:-1].split(',')
			vertex_list[0] = int(vertex_list[0])
			vertex_list[1] = int(vertex_list[1])
			vertex_list[2] = int(vertex_list[2])
			self.recurrence_triplet_list.append(vertex_list)
			line = inf.readline()
		inf.close()	
		
	def recurrence_stat_list_construct(self):
		no_of_triplets = len(self.recurrence_triplet_list)
		if no_of_triplets < 100000:
			sys.stderr.write('\tonly %d triplets\n'%no_of_triplets)
		if no_of_triplets <1000:
			sys.stderr.write('\tthe number of triplets is below 1000(%d). Aborted.\n'%no_of_triplets)
			#sys.exit(1)
		sys.stderr.write('func\ttrans\n')
		for j in xrange(100):
			functional_homo = 0
			transcriptional_homo = 0
			index_list = random.sample(xrange(no_of_triplets),100)
			for k in xrange(100):
				triplet = self.recurrence_triplet_list[index_list[k]]
				if self.is_homogenious(triplet, self.known_genes_dict):
					functional_homo += 1
				if self.is_homogenious(triplet, self.transfac_dict):
					transcriptional_homo += 1
			functional_homo_ratio = functional_homo/100.00
			transcriptional_homo_ratio = transcriptional_homo/100.00
			sys.stderr.write('%f\t%f\n'%\
				(functional_homo_ratio,transcriptional_homo_ratio))
			self.recurrence_stat_func_list.append(functional_homo_ratio)
			self.recurrence_stat_trans_list.append(transcriptional_homo_ratio)
		
		
	def is_homogenious(self, triplet, dict):
		list = dict[triplet[0]] + dict[triplet[1]] + dict[triplet[2]]
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
	def helper():
		sys.stderr.write('\
	argv[1] is the file containing triplets with specific frequency.\n')
	if len(sys.argv) == 2:
		instance = triplet_stat()
		instance.dstruc_loadin()
		instance.recurrence_triplet_list_construct(sys.argv[1])
		instance.recurrence_stat_list_construct()
	else:
		helper()
		sys.exit(1)
	'''
	# this block is for transfac_dict construction.
	inf = open(sys.argv[1])
	instance.transfac_dict_construct(inf)
	key_list = instance.transfac_dict.keys()
	key_list.sort()
	for item in key_list:
		print '%s\t%s'%(item,instance.transfac_dict[item],)
	'''
