#!/usr/bin/env python
"""
Usage: triplet_stat.py -k SCHEMA [OPTION] TRIPLET_FILE

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-h, --help              show this help
	
Examples:
	triplet_stat.py -k ming2 gph_result/triplet_6

Description:
	Triplet never means triplet. Triplet can be any large(>3).
	The name is history relics.
"""

import pickle,sys,os,random,getopt,psycopg,csv


class triplet_stat:
	def __init__(self, hostname, dbname, schema):
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to %s"%schema)
		self.recurrence_triplet_list = []
		self.recurrence_stat_func_list = []
		self.recurrence_stat_trans_list = []
		self.transfac_pickle_fname = os.path.join(os.path.expanduser('~'),'pickle/yeast_transfac_dict')
		self.vertex_dict = {}
		self.transfac_dict = {}
		self.known_genes_dict = {}
		self.logfile = open('/tmp/triplet_stat.log','w')
	
	def dstruc_loadin(self):
		if not os.path.isfile(self.transfac_pickle_fname):
			sys.stderr.write('You need to construct the transfac_dict first.\n')
			sys.exit(1)
		
		self.curs.execute("select gene_id, gene_no from gene")
		rows = self.curs.fetchall()
		for row in rows:
			self.vertex_dict[row[0]] = row[1]
		
		transfac_dict_orf = pickle.load(open(self.transfac_pickle_fname,'r'))
		#replace the orfname with no. because the triplets are in no. form.
		for item in transfac_dict_orf.iteritems():
			if self.vertex_dict.has_key(item[0]):
				no = self.vertex_dict[item[0]]
				self.transfac_dict[no] = item[1]
		
		self.curs.execute("select gene_no,go_functions from gene where known=TRUE")
		rows = self.curs.fetchall()
		for row in rows:
			if row[1] == '{}':
				sys.stderr.write('No function for gene: %d\n'%row[0])
				continue
			go_list = row[1][1:-1].split(',')
			self.known_genes_dict[row[0]] = go_list


	
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
		reader = csv.reader(inf)
		for row in reader:
			row = map(int, row)
			self.recurrence_triplet_list.append(row)
		inf.close()
		
	def recurrence_stat_list_construct(self):
		no_of_triplets = len(self.recurrence_triplet_list)
		if no_of_triplets < 100000:
			sys.stderr.write('\tonly %d triplets\n'%no_of_triplets)
		if no_of_triplets <1000:
			sys.stderr.write('\tthe number of triplets is below 1000(%d). Aborted.\n'%no_of_triplets)
			#sys.exit(1)
		sys.stdout.write('func\ttrans\n')
		for j in xrange(100):
			functional_homo = 0
			transcriptional_homo = 0
			index_list = random.sample(xrange(no_of_triplets),1000)
			for k in xrange(1000):
				triplet = self.recurrence_triplet_list[index_list[k]]
				if self.is_homogenious(triplet, self.known_genes_dict):
					self.logfile.write('f %s\n'%repr(triplet))
					functional_homo += 1
				if self.is_homogenious(triplet, self.transfac_dict):
					self.logfile.write('t %s\n'%repr(triplet))
					transcriptional_homo += 1
			functional_homo_ratio = functional_homo/1000.00
			transcriptional_homo_ratio = transcriptional_homo/1000.00
			sys.stdout.write('%f\t%f\n'%\
				(functional_homo_ratio,transcriptional_homo_ratio))
			self.recurrence_stat_func_list.append(functional_homo_ratio)
			self.recurrence_stat_trans_list.append(transcriptional_homo_ratio)
		
		
	def is_homogenious(self, triplet, dict):
		list = []
		for vertex in triplet:
			if vertex not in dict:
			#for the transfac_dict, some known genes' trans-factors are unknown.
				return 0
			else:
				list += dict[vertex]
		no_of_vertices = len(triplet)
		judge_dict = {}
		for item in list:
			if judge_dict.has_key(item):
				judge_dict[item] += 1
			else:
				judge_dict[item] = 1
			if judge_dict[item] == no_of_vertices:
				return 1
		return 0


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:", ["help", "hostname=", "dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-z", "--hostname"):
			hostname = arg
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-k", "--schema"):
			schema = arg
			
	if schema and len(args) ==1:
		instance = triplet_stat(hostname, dbname, schema)
		instance.dstruc_loadin()
		instance.recurrence_triplet_list_construct(args[0])
		instance.recurrence_stat_list_construct()	
	else:
		print __doc__
		sys.exit(2)
	'''
	# this block is for transfac_dict construction.
	inf = open(sys.argv[1])
	instance.transfac_dict_construct(inf)
	key_list = instance.transfac_dict.keys()
	key_list.sort()
	for item in key_list:
		print '%s\t%s'%(item,instance.transfac_dict[item],)
	'''
