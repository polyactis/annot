#!/usr/bin/env python
import pickle,sys,os
import Martel
from xml.sax import saxutils
from Martel import LAX

termtype = Martel.Str("TERMTYPE:[") + Martel.ToSep("termtype", "]") + Martel.AnyEol()
id = Martel.Str("id:") + Martel.ToEol("id")
name = Martel.Str("name:") + Martel.ToEol("name")
alt_id = Martel.Str("alt_id:") + Martel.ToEol("alt_id")
namespace = Martel.Str("namespace:") + Martel.ToEol("namespace")
definition = Martel.Str("def:") + Martel.ToEol("def")
comment = Martel.Str("comment:") + Martel.ToEol("comment")
is_a = Martel.Str("is_a:") + Martel.ToEol("is_a")
no_of_genes = Martel.Digits("no_of_genes")
genes = Martel.ToEol("genes")
dataline = Martel.Group("dataline", no_of_genes+\
												Martel.Str(";")+\
												genes)
												
term = Martel.Group("term",  Martel.Rep(definition)+\
										Martel.Rep(is_a)+\
										termtype+\
										id+\
										name+\
										Martel.Rep(alt_id)+\
										namespace+\
										Martel.Rep(definition)+\
										Martel.Rep(comment)+\
										Martel.Rep(is_a)+\
										dataline+\
										Martel.Rep(Martel.AnyEol()) )
go_bioprocess = Martel.Group("go_bioprocess",  Martel.Rep(term), {"format":"go"} )
														
class go_term:
	'''
	a structure for holding information related to a GeneOntology term.
	'''
	def __init__(self):
		no = None
		name = None
		no_of_genes = None
		whole_gene_array = None
		gene_array = None


class go_bioprocess_parser:
	'''
	parse the yeast known_gene(GO bioprocess) file,
	construct a known-gene dictionary with gene name as key and GO id as value.
	'''
	def __init__(self, infname):
		self.go_dict = {}
		self.global_struc_fname = os.path.join(os.path.expanduser('~'),'pickle/yeast_global_struc')
		if os.path.isfile(self.global_struc_fname):
			global_struc = pickle.load(open(self.global_struc_fname,'r'))
			self.vertex_dict = global_struc['vertex_dict']
		else:
			sys.stderr.write('No such file: %s.\n'%self.global_struc_fname)
			sys.exit(1)
		self.inf = open(infname, 'r')
		
	def parse(self):
		self.term_iterator = go_bioprocess.make_iterator("term")
		self.lax=LAX.LAX()
		for record in self.term_iterator.iterateFile(self.inf, self.lax):
			gene_array = []
			whole_gene_array = record['genes'][0].split(';')
			for gene in whole_gene_array:
				if gene in self.vertex_dict:
					gene_array.append(self.vertex_dict[gene])
			gene_array.sort()
			info = go_term()
			info.name = record['name'][0]
			info.no_of_genes = int(record['no_of_genes'][0])
			info.whole_gene_array = whole_gene_array
			info.gene_array = gene_array
			self.go_dict[record['id'][0]] = info
			
		key_list = self.go_dict.keys()
		key_list.sort()
		for i in range(len(key_list)):
			#output in tab delimited format for database import
			id = key_list[i]
			string1 = repr(self.go_dict[id].whole_gene_array)
			string1 = string1.replace('[','{')
			string1 = string1.replace(']','}')
			string1 = string1.replace("'",'"')
			string2 = repr(self.go_dict[id].gene_array)
			string2 = string2.replace('[','{')
			string2 = string2.replace(']','}')
			print "%s\t%d\t%d\t%s\t%s\t%s"%\
				(id,\
				i+1,\
				self.go_dict[id].no_of_genes,\
				self.go_dict[id].name,\
				string1,\
				string2)
		#parser = go_bioprocess.make_parser()
		#parser.setContentHandler(saxutils.XMLGenerator())
		#parser.parseFile(self.inf)

if __name__ == '__main__':
	def helper():
		sys.stderr.write('\
	argv[1] specifies the go bioprocess raw file.\n\
		output will be dumped on stdout.\n')

	if len(sys.argv) == 2:
		instance = go_bioprocess_parser(sys.argv[1])
		instance.parse()
	else:
		helper()
