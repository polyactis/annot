#!/usr/bin/env python
import sys,os,psycopg,cStringIO


class mcl_result_iterator:
	'''looping over a mcl result file, generate a block of clusters related to one splat pattern'''
	def __init__(self, inf):
		self.inf = inf
		self.mcl_cluster_block = ''
	def __iter__(self):
		return self		
	def next(self):
		self.read()
		return cStringIO.StringIO(self.mcl_cluster_block)
	def read(self):
		self.mcl_cluster_block = ''
		mcl_begin = 0
		line = self.inf.readline()
		while line != '\n':
			'''
			\n
			(mclruninfo
			...
			is discarded.
			'''
			if line == '':
				raise StopIteration
				break
			if mcl_begin == 1:
				self.mcl_cluster_block += line
			if line.find('(splat_id') == 0:
				mcl_begin = 1
				self.mcl_cluster_block += line
			line = self.inf.readline()
		self.mcl_cluster_block += line

class mcl_to_db:
	def __init__(self, dir, dbname, needcommit=0):
		self.dir = os.path.abspath(dir)
		self.conn = psycopg.connect('dbname=%s'%dbname)
		self.needcommit = int(needcommit)
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to graph")
		self.splat_id = ''
		self.cluster_set = []
		self.parameter = ''
				
	def submit(self):
		f_list = os.listdir(self.dir)
		no = 0
		for fname in f_list:
			if fname.find('out') == 0:
				path = os.path.join(self.dir, fname)
				inf = open(path, 'r')
				iter = mcl_result_iterator(inf)
				for mcl_cluster_block in iter:
					self.parse(mcl_cluster_block)
					for i in range(len(self.cluster_set)):

						vertex_set = self.cluster_set[i]
						self.mcl_id = '%s_%s_%d'%(self.splat_id,self.parameter,(i+1))
						string_vertex_set = repr(vertex_set)
						string_vertex_set = string_vertex_set.replace('[','{')
						string_vertex_set = string_vertex_set.replace(']','}')
						try:
							self.curs.execute("insert into mcl_result(splat_id, vertex_set, parameter) \
							values ('%s','%s','%s')"%(self.splat_id, string_vertex_set, self.parameter))
						except:
							sys.stderr.write('Error occured when inserting pattern. Aborted.\n')
							self.conn.rollback()
							sys.exit(1)
						no+=1
					string = repr(no)
					sys.stderr.write('%s%s'%('\x08'*(len(string)+1), string))
		if self.needcommit:
			self.conn.commit()
		sys.stderr.write('\n\tTotal clusters: %d\n'%no)

	def parse(self,inf):
		line = inf.readline()
		self.cluster_set = []
		cluster_begin = 0
		while line:
			if line.find('(splat_id') == 0:
				self.splat_id = line[10:-3]
			if line.find('(parameter') == 0:
				self.parameter = line[11:-3]
			if line == 'begin\n':
				cluster_begin = 1
			if cluster_begin and (line == ')\n'):
				break;
			elif cluster_begin:
				vertex_list = line[1:-2].split()	#remove the first column and $\n.
				if len(vertex_list) >2:
					for i in range(len(vertex_list)):
						vertex_list[i] = int(vertex_list[i])
					self.cluster_set.append(vertex_list)
			line = inf.readline()
				
if __name__ == '__main__':
	def helper():
		sys.stderr.write('\
	argv[1] is the directory containing mcl results.\n\
	argv[2] is the database name.\n\
	argv[3] is 1 or 0 indicating whether to commit or not. Default is 0.\n')
		
	if len(sys.argv) == 4:
		instance = mcl_to_db(sys.argv[1], sys.argv[2], sys.argv[3])
	elif len(sys.argv) == 3:
		instance = mcl_to_db(sys.argv[1], sys.argv[2])
	else:
		helper()
		sys.exit(1)
	instance.submit()
