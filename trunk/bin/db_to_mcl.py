#!/usr/bin/env python
import sys,os,cStringIO,psycopg


class db_to_mcl:
	def __init__(self, dbname):
		self.vertex_dict = {}
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
					
	def get_from_db(self, organism, outdir):
		self.curs.execute("select splat_id,edge_set from splat_result where \
			organism='%s'"%(self.org_short2long[organism]))
		if not os.path.isdir(outdir):
			os.makedirs(outdir)
		row = self.curs.fetchone()
		while row:
			splat_id = row[0]
			outfname = splat_id
			path = os.path.join(outdir, outfname)
			outf = open(path, 'w')
			self.vertex_dict = {}	#intialize the vertex_dict
			self.mcl_singleton_write(row[1], outf)
			row = self.curs.fetchone()
			
	def mcl_singleton_write(self, edge_set, outf):
		edge_list = edge_set[2:-2].split('},{')
		for edge in edge_list:
			vertex_list = edge.split(',')
			vertex1 = vertex_list[0]
			vertex2 = vertex_list[1]
			if vertex1 in self.vertex_dict:
				self.vertex_dict[vertex1].append(vertex2)
			else:
				self.vertex_dict[vertex1] = [vertex2]
			if vertex2 in self.vertex_dict:
				self.vertex_dict[vertex2].append(vertex1)
			else:
				self.vertex_dict[vertex2] = [vertex1]
		dim = len(self.vertex_dict)
		outf.write('(mclheader\n')
		outf.write('mcltype matrix\n')
		outf.write('dimensions %dx%d\n)\n'%(dim,dim,))
		outf.write('(mcldoms\n')
		vertex_list = self.vertex_dict.keys()
		outf.write('%s $\n)\n'%' '.join(vertex_list))
		outf.write('(mclmatrix\nbegin\n')
		for vertex in self.vertex_dict:
			outf.write('%s '%vertex)
			outf.write('%s $\n'%' '.join(self.vertex_dict[vertex]))
		outf.write(')\n')
		
if __name__ == '__main__':
	def helper():
		sys.stderr.write('\
		argv[1] is the output directory.\n\
		argv[2] is the database name.\n\
		argv[3] is two-letter organism abbreviation.\n')
		
	if len(sys.argv) == 4:
		instance = db_to_mcl(sys.argv[2])

		instance.get_from_db(sys.argv[3], sys.argv[1])
	else:
		helper()
