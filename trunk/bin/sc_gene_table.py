#!/usr/bin/env python
import pickle,sys,os,psycopg


class sc_gene:
	'''
	Initialize the gene_table of yeast from two pickled data structures,
		known_genes_dict
		global_struc['vertex_list']
	'''
	def __init__(self, dbname, needcommit=0):
		self.conn = psycopg.connect('dbname=%s'%dbname)
		self.curs = self.conn.cursor()
		self.needcommit = int(needcommit)
		self.known_genes_dict_fname = os.path.join(os.path.expanduser('~'), 'pickle/known_genes_dict')
		self.global_struc_fname = os.path.join(os.path.expanduser('~'),'pickle/yeast_global_struc')
		self.go_id_to_no_dict = {}

	def dstruc_loadin(self):
		self.curs.execute("select go_id, go_no from graph.sc_go")
		rows = self.curs.fetchall()
		for row in rows:
			self.go_id_to_no_dict[row[0]] = row[1]
		if os.path.isfile(self.known_genes_dict_fname):
			self.known_genes_dict = pickle.load(open(self.known_genes_dict_fname,'r'))
		else:
			sys.stderr.write('No known_genes_dict pickle file\n')
			sys.exit(1)
		
		if os.path.isfile(self.global_struc_fname):
			global_struc = pickle.load(open(self.global_struc_fname,'r'))
			self.vertex_list = global_struc['vertex_list']
		else:
			sys.stderr.write('No yeast_global_struc pickle file\n')
			sys.exit(1)

	def run(self):
		for i in range(len(self.vertex_list)):
			orfname = self.vertex_list[i]
			gene_no = i+1
			known = '0'
			go_functions_list = []
			if orfname in self.known_genes_dict:
				known = '1'
				for go_id in self.known_genes_dict[orfname]:
					no = self.go_id_to_no_dict[go_id]
					go_functions_list.append(no)
			if go_functions_list:
				self.curs.execute("insert into graph.sc_gene(orfname, gene_no, known, \
				go_functions) values ('%s', %d, '%s', ARRAY%s)"%\
				(orfname, gene_no, known,repr(go_functions_list)))
			else:
				self.curs.execute("insert into graph.sc_gene(orfname, gene_no, known)\
				values ('%s', %d, '%s')"%\
				(orfname, gene_no, known))
		if self.needcommit:
			self.conn.commit()


if __name__ == '__main__':
	def helper():
		sys.stderr.write('\
	argv[1] specifies the database name.\n\
	argv[2] is 1 or 0 meaning to commit or not. Ignore is 0.\n\
		This program sets up the sc_gene table.\n')
	if len(sys.argv) == 2:
		instance = sc_gene(sys.argv[1])
		instance.dstruc_loadin()
		instance.run()
	elif len(sys.argv) == 3:
		instance = sc_gene(sys.argv[1],sys.argv[2])
		instance.dstruc_loadin()
		instance.run()
	else:
		helper()
		sys.exit(1)
