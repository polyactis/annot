#!/usr/bin/python
import sys,os,cStringIO,psycopg

class splat_result_iterator:
	'''looping over a splat result file, generate a single pattern result'''
	def __init__(self, inf):
		self.inf = inf
		self.pattern = ''
	def __iter__(self):
		return self		
	def next(self):
		self.read()
		return cStringIO.StringIO(self.pattern)
	def read(self):
		self.pattern = ''
		line = self.inf.readline()
		while line != '\n':
			if line == '':
				raise StopIteration
				break
			self.pattern += line
			line = self.inf.readline()
		self.pattern += line
	
class splat_to_db:
	'''
	'''
	def __init__(self, infname, dbname, organism, needcommit=0):
		self.splat_id = ''
		self.organism = organism
		self.no_of_edges = ''
		self.recurrence_pattern = ''
		self.edge_set = []
		self.inf = open(infname, 'r')
		self.needcommit = needcommit
		self.conn = psycopg.connect('dbname=%s'%dbname)
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to graph")
		self.org_short2long = {'at':'Arabidopsis thaliana',
					'ce':'Caenorhabditis elegans',
					'dm':'Drosophila melanogaster',
					'hs':'Homo sapiens',
					'mm':'Mus musculus',
					'sc':'Saccharomyces cerevisiae',
					'Arabidopsis thaliana':'Arabidopsis thaliana',
					'Caenorhabditis elegans':'Caenorhabditis elegans',
					'Drosophila melanogaster':'Drosophila melanogaster',
					'Homo sapiens':'Homo sapiens',
					'Mus musculus':'Mus musculus',
					'Gorilla gorilla Pan paniscus Homo sapiens':'Homo sapiens',
					'Saccharomyces cerevisiae':'Saccharomyces cerevisiae'}
					
	def parse(self, pattern):
		line = pattern.readline()
		no_in_string = line[:-1]
		self.no_of_edges = no_in_string
		line = pattern.readline()
		self.recurrence_pattern = line[:-1]
		line = pattern.readline()
		self.edge_set = []	#initialize the edge_set structure
		while line != '\n':
			list = line.split(')')
			list.pop()	#throw away the last '\n'
			for item in list:
				l_list = item[1:].split()
				edge = [int(l_list[0]), int(l_list[1])]
				self.edge_set.append(edge)
			line = pattern.readline()
			
	def run(self):
		iter = splat_result_iterator(self.inf)
		no = 1
		for pattern in iter:
			self.parse(pattern)
			self.splat_id = '%s_%d'%(self.organism,no)
			sys.stderr.write(self.splat_id + '\t' + self.no_of_edges + '\n' + self.recurrence_pattern +'\n')
			try:
				self.curs.execute("insert into splat_result(splat_id, organism,no_of_edges, \
							recurrence_pattern,edge_set)values ('%s','%s','%s',B'%s',ARRAY%s)"%\
							(self.splat_id, self.org_short2long[self.organism], self.no_of_edges, \
							self.recurrence_pattern,repr(self.edge_set) ))
			except:
				sys.stderr.write('Error occured when inserting pattern. Aborted.\n')
				self.conn.rollback()
				sys.exit(1)
			no+=1
		if self.needcommit:
			self.conn.commit()


if __name__ == '__main__':
	def helper():
		sys.stderr.write('\
		argv[1] is the splat result file.\n\
		argv[2] is the database name.\n\
		argv[3] is the two abbreviation letters for organism.\n\
		argv[4] is 0 or 1 indicating whether to commit or not.\n')
		
	if len(sys.argv) ==5:
		instance = splat_to_db(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
	elif len(sys.argv) == 4:
		instance = splat_to_db(sys.argv[1], sys.argv[2], sys.argv[3], 0)
	else:
		helper()
		sys.exit(1)
	instance.run()
