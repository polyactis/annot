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
														
															

class known_genes:
	'''
	parse the yeast known_gene(GO bioprocess) file,
	construct a known-gene dictionary with gene name as key and GO id as value.
	'''
	def __init__(self, inf):
		self.known_genes_dict = {}
		self.inf = inf
		
	def parse(self):
		self.term_iterator = go_bioprocess.make_iterator("term")
		self.lax=LAX.LAX()
		for record in self.term_iterator.iterateFile(inf, self.lax):
			gene_list = record['genes'][0].split(';')
			for gene in gene_list:
				if gene in self.known_genes_dict:
					self.known_genes_dict[gene].append(record['id'][0])
				else:
					self.known_genes_dict[gene] = [ record['id'][0] ]

		for gene in self.known_genes_dict:
			self.known_genes_dict[gene].sort()
			print '%s\t%s\n'%(gene, self.known_genes_dict[gene],)
		#parser = go_bioprocess.make_parser()
		#parser.setContentHandler(saxutils.XMLGenerator())
		#parser.parseFile(self.inf)

class yeast_two_way_filter:
	'''
	Read yeast raw datafile, first remove the single quotes around the orfname.
	secondly, output to two files.
		normal_outf is the same as the raw datafile but without the single quotes.
		known_genes_outf contains only the edges of known genes.
	'''
	def __init__(self, known_genes_dict):
		self.dict = known_genes_dict
	
	def filter_based_on_known_genes_dict(self, inf, outf1, outf2):
		self.normal_outf = outf1
		self.known_genes_outf = outf2
		line = inf.readline()
		while line:
			if line[0] == 't':
				self.normal_outf.write(line)
				self.known_genes_outf.write(line)
			if line[0] == 'e':
				list = line[:-1].split('\t')
				list[1] = list[1][1:-2]
				list[2] = list[2][1:-2]
				vertex1 = list[1]
				vertex2 = list[2]
				self.normal_outf.write('%s\n'%'\t'.join(list))
				if (vertex1 in self.dict) and (vertex2 in self.dict):
					self.known_genes_outf.write('%s\n'%'\t'.join(list))
			line = inf.readline()

def batch(dir, output_dir):
	files = os.listdir(dir)
	sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
	normal_dir = os.path.join(output_dir,'sc_normal')
	if not os.path.isdir(normal_dir):
		os.makedirs(normal_dir)
	known_dir = os.path.join(output_dir,'sc_known')
	if not os.path.isdir(known_dir):
		os.makedirs(known_dir)
	pickle_file = os.path.join(os.path.expanduser('~'), 'pickle/known_genes_dict')
	known_genes_dict = pickle.load(open(pickle_file,'r'))
	instance = yeast_two_way_filter(known_genes_dict)
	
	for f in files:
		pathname = os.path.join(dir, f)
		sys.stderr.write("%d/%d:\t%s\n"%(files.index(f)+1,len(files),f))
		inf = open(pathname, 'r')
		normal_f = open(os.path.join(normal_dir,f), 'w')
		known_f = open(os.path.join(known_dir, f), 'w')
		instance.filter_based_on_known_genes_dict(inf, normal_f, known_f)
		inf.close()
		normal_f.close()
		known_f.close()



if __name__ == '__main__':
	batch(sys.argv[1],sys.argv[2])
	# argv[1] specifies which directory contains raw yeast data
	# argv[2] specifies the directory to store two bunchs of datafiles
	'''
	# this part is for constructing known_genes_dict
	inf = open(sys.argv[1], 'r')
	pickle_file = os.path.join(os.path.expanduser('~'), 'pickle/known_genes_dict')
	instance = known_genes(inf)
	instance.parse()
	pickle.dump ( instance.known_genes_dict, open(pickle_file,'w') )
	'''
