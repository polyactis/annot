#!/usr/bin/env python
import sys,os,psycopg

class mcl_to_db:
	def __init__(self, dir, dbname, needcommit=0):
		self.dir = os.path.abspath(dir)
		self.conn = psycopg.connect('dbname=%s'%dbname)
		self.needcommit = int(needcommit)
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to graph")
		self.mcl_id = ''
		self.splat_id = ''
		self.cluster_set = []
		self.parameter = ''
				
	def submit(self):
		f_list = os.listdir(self.dir)
		no = 0
		for fname in f_list:
			if fname.find('out') == 0:
				fname_list = fname.split('.')
				self.splat_id = fname_list[1]
				self.parameter = fname_list[2]
				path = os.path.join(self.dir, fname)
				inf = open(path, 'r')
				self.cluster_set_construct(inf)
				for i in range(len(self.cluster_set)):
					#try:
					vertex_set = self.cluster_set[i]
					self.mcl_id = '%s_%s_%d'%(self.splat_id,self.parameter,(i+1))
					string_vertex_set = repr(vertex_set)
					string_vertex_set = string_vertex_set.replace('[','{')
					string_vertex_set = string_vertex_set.replace(']','}')
					self.curs.execute("insert into mcl_result(mcl_id, splat_id, vertex_set, parameter) \
						values ('%s','%s','%s','%s')"%(self.mcl_id, self.splat_id, string_vertex_set, self.parameter))
					#except:
						#sys.stderr.write('Error occured when inserting pattern. Aborted.\n')
						#self.conn.rollback()
						#sys.exit(1)
					no+=1
					sys.stderr.write('.')
		if self.needcommit:
			self.conn.commit()
		sys.stderr.write('\n\tTotal clusters: %d\n'%no)

	def cluster_set_construct(self,inf):
		line = inf.readline()
		self.cluster_set = []
		cluster_begin = 0
		while line:
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
