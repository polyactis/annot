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
					
	def get_from_db(self, outdir, organism, package_size):
		if not os.path.isdir(outdir):
			os.makedirs(outdir)
		self.curs.execute("begin")
		self.curs.execute("DECLARE crs CURSOR FOR select splat_id,edge_set from splat_result \
			where organism='%s'"%(self.org_short2long[organism]))
		self.curs.execute("fetch %d from crs"%(package_size))
		rows = self.curs.fetchall()
		no_of_splat_records = 0
		no_of_packages = 0
		while rows:
			outfname = 'mcl_input%d'%no_of_packages
			path = os.path.join(outdir, outfname)
			outf = open(path, 'w')
			for row in rows:
				self.mcl_singleton_write(row, outf)
				no_of_splat_records += 1
			outf.close()
			no_of_packages += 1
			string = 'Packages: %d\t Splat_patterns: %d'%(no_of_packages, no_of_splat_records)
			sys.stderr.write('%s%s'%('\x08'*80,string))
			self.curs.execute("fetch %d from crs"%(package_size))
			rows = self.curs.fetchall()
		self.curs.execute("end")
		sys.stderr.write('\n\tTotal patterns: %d\n'%no_of_splat_records)
		
	def mcl_singleton_write(self, row, outf):
		self.vertex_dict = {}	#intialize the vertex_dict
		splat_id = row[0]
		edge_set = row[1]
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
		out_block = '(splat_id %s )\n'%splat_id		# here it is '=' not '+='
		out_block += '(mclheader\n'
		out_block += 'mcltype matrix\n'
		out_block += 'dimensions %dx%d\n)\n'%(dim,dim)
		out_block += '(mcldoms\n'
		vertex_list = self.vertex_dict.keys()
		out_block += '%s $\n)\n'%' '.join(vertex_list)
		out_block += '(mclmatrix\nbegin\n'
		for vertex in self.vertex_dict:
			out_block += '%s '%vertex
			out_block += '%s $\n'%' '.join(self.vertex_dict[vertex])
		out_block += ')\n\n'
		outf.write(out_block)
		
if __name__ == '__main__':
	def helper():
		sys.stderr.write('\
		argv[1] is the output directory.\n\
		argv[2] is the database name.\n\
		argv[3] is two-letter organism abbreviation.\n\
		argv[4] specifies the number of splat records in one output file.\n')
		
	if len(sys.argv) == 5:
		instance = db_to_mcl(sys.argv[2])
		instance.get_from_db(sys.argv[1], sys.argv[3], int(sys.argv[4]))
	else:
		helper()
