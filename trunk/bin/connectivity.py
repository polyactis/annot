#!/usr/bin/env python
import sys,os,cStringIO,psycopg

class compute_connectivity:
	'''
	a class for computing the connectivity of splat and mcl results.
	1000 records are done in one time. Reduce the memory usage.
	This is done through postgresql's CURSOR mechanism.
	Computing connectivity of mcl results is based on splat_id.
	Because edge_dict is required in advance.
	'''
	def __init__(self, table, dbname, needcommit=0):
		self.table = table
		self.conn = psycopg.connect('dbname=%s'%dbname)
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to graph")
		self.needcommit = int(needcommit)
		self.vertex_dict = {}
		self.edge_dict = {}
		self.no_of_splat_records = 0
		self.no_of_mcl_records = 0
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
					
	def run(self, organism):
		self.curs.execute("begin")
		self.curs.execute("DECLARE crs CURSOR FOR select splat_id,edge_set from splat_result \
			where organism='%s'"%(self.org_short2long[organism]))
		self.curs.execute("fetch 1000 from crs")
		rows = self.curs.fetchall()
		while rows:
			for row in rows:
				splat_id = row[0]
				edge_set = row[1]
				if self.table == 'splat':
					self.splat_atom_update(splat_id, edge_set)
				elif self.table == 'mcl':
					self.mcl_atom_update(splat_id, edge_set)
			
			if self.table == 'splat':
				sys.stderr.write("%s%s"%("\x08"*80,self.no_of_splat_records))
			elif self.table == 'mcl':
				sys.stderr.write("%s%s"%("\x08"*80,self.no_of_mcl_records))
			self.curs.execute("fetch 1000 from crs")
			rows = self.curs.fetchall()
		if self.needcommit:
			self.curs.execute("end")
			sys.stderr.write('\n\tTotal %d splat records updated\n'%self.no_of_splat_records)
			sys.stderr.write('\tTotal %d mcl records updated\n'%self.no_of_mcl_records)
		else:
			self.conn.rollback()
			sys.stderr.write('\n\tNo real updates\n')
	
	def mcl_atom_update(self, splat_id, edge_set):
		'''
		Compute the connectivity of all mcl clusters pertaining to the same splat_id.
		'''
		self.dstruc_from_edge_set(edge_set)
		self.curs.execute("select mcl_id, vertex_set from mcl_result where splat_id='%s'"%\
			(splat_id))
		rows = self.curs.fetchall()
		for row in rows:
			mcl_id = row[0]
			vertex_set = row[1][1:-1]
			vertex_list = vertex_set.split(',')
			no_of_vertices = len(vertex_list)
			no_of_edges = 0
			for i in xrange(no_of_vertices):
				for j in xrange(i+1, no_of_vertices):
					tuple = (vertex_list[i],vertex_list[j])
					if self.edge_dict.has_key(tuple) or self.edge_dict.has_key(tuple):
						no_of_edges += 1
			connectivity = 2.0*no_of_edges/((no_of_vertices-1)*no_of_vertices)
			try:
				self.curs.execute("update mcl_result set connectivity=%f where mcl_id='%s'"% \
				(connectivity, mcl_id))
			except:
				sys.stderr.write('Error occurred while setting mcl connectivity\n')
				sys.exit(1)
			self.no_of_mcl_records += 1

	def splat_atom_update(self, splat_id, edge_set):
		self.dstruc_from_edge_set(edge_set)
		no_of_edges = len(self.edge_dict)
		no_of_vertices = len(self.vertex_dict)
		connectivity = 2.0*no_of_edges/((no_of_vertices-1)*no_of_vertices)
		try:
			self.curs.execute("update splat_result set connectivity=%f where splat_id='%s'"% \
			(connectivity, splat_id))
		except:
			sys.stderr.write('Error occurred while setting splat connectivity\n')
			sys.exit(1)
		self.no_of_splat_records += 1
	
	def dstruc_from_edge_set(self, edge_set):
		self.edge_dict = {}
		self.vertex_dict = {}
		edge_list = edge_set[2:-2].split('},{')
		for edge in edge_list:
			vertex_list = edge.split(',')
			vertex_list = (vertex_list[0], vertex_list[1])
			self.edge_dict[vertex_list] = 1
			vertex1 = vertex_list[0]
			vertex2 = vertex_list[1]
			if vertex1 not in self.vertex_dict:
				self.vertex_dict[vertex1] = 1
			if vertex2 not in self.vertex_dict:
				self.vertex_dict[vertex2] = 1
		

if __name__ == '__main__':
	def helper():
		sys.stderr.write('\
	argv[1] could be splat or mcl, specifying which connectivity to compute.\n\
	argv[2] is the database name\n\
	argv[2] is two-letter abbreviation for an organism.\n\
	argv[4] is 1 or 0 indicating whether to commit or not. Default is 0.\n')
	
	if len(sys.argv) == 4:
		instance = compute_connectivity(sys.argv[1], sys.argv[2])
		instance.run(sys.argv[3])
	elif len(sys.argv) == 5:
		instance = compute_connectivity(sys.argv[1], sys.argv[2],sys.argv[4])
		instance.run(sys.argv[3])
	else:
		helper()
