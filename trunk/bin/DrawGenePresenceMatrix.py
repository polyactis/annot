#!/usr/bin/env python
"""
Usage: DrawGenePresenceMatrix.py -i INPUTDIR -o OUTPUTFILE -p OUTPUTPIC [OPTIONS]

Option:
	-i ...	the input directory
	-o ...	the incidence output_fname
	-p ...	the jpg output_fname
	-b	enable debugging, no debug by default
	-r	report the progress(a number)
	-h, --help              show this help

Examples:
	DrawGenePresenceMatrix.py -i ~/datasets/hs_fim_92/
		-o /tmp/yuhuang/hs_fim_92.gim -p /tmp/yuhuang/test.png

Description:
	Files in INPUTDIR are sorted by file_rename.py(datasets_sort).
	Program to draw a picture to show how genes are located among datasets.
	Hint: png is better than jpeg. jpeg has maximum 65500 pixels.
"""

import os, sys, csv, getopt
sys.path += [os.path.join(os.path.expanduser('~/script/annot/bin'))]
from graph.complete_cor_vector import complete_cor_vector
from sets import Set
import Image, ImageDraw

class DrawGenePresenceMatrix:
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None,\
		input_dir=None, output_fname=None, jpg_ofname=None, debug=0, report=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_dir = input_dir
		self.output_fname = output_fname
		self.jpg_ofname = jpg_ofname
		self.debug = int(debug)
		self.report = int(report)
	
	def expand_gene_id2presence_vector(self, gene_id2presence_vector, filename, file_index):
		"""
		09-26-05
			fix a bug, non-present genes' presence_vector should have a trailing 0
		"""
		reader = csv.reader(file(filename), delimiter='\t')
		for row in reader:
			gene_id = row[0]
			if gene_id not in gene_id2presence_vector:
				gene_id2presence_vector[gene_id] = [0]*file_index + [1]
			else:
				gene_id2presence_vector[gene_id].append(1)
		#add the trailing 0
		for gene_id in gene_id2presence_vector:
			if len(gene_id2presence_vector[gene_id])==file_index:
				gene_id2presence_vector[gene_id].append(0)
		del reader
	
	def write_gene_incidence_matrix(self, dir, output_fname):
		sys.stderr.write("Writing gene presence_vector into %s...\n"%output_fname)
		complete_cor_vector_instance = complete_cor_vector()
		gene_id2presence_vector = {}
		
		files = os.listdir(dir)
		files = complete_cor_vector_instance.files_sort(files)
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
		
		for i in range(len(files)):
			f = files[i]
			print f
			f_path = os.path.join(dir, f)
			self.expand_gene_id2presence_vector(gene_id2presence_vector, f_path, i)
			
		frequency_presence_vector_gene_id_ls =[]
		for gene_id, presence_vector in gene_id2presence_vector.iteritems():
			frequency_presence_vector_gene_id_ls.append([sum(presence_vector)] + presence_vector + [gene_id])
		
		del gene_id2presence_vector
		frequency_presence_vector_gene_id_ls.sort()	#sorted according to frequency
		frequency_presence_vector_gene_id_ls.reverse()	#09-26-05 reverse it
		
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		for row in frequency_presence_vector_gene_id_ls:
			writer.writerow(row)
		del writer
		sys.stderr.write("Done.\n")
		
	def get_char_dimension(self):
		im = Image.new('RGB', (50,50))
		draw = ImageDraw.Draw(im)
		char_dimension = draw.textsize('a')
		del im, draw
		return char_dimension
	
	def get_text_region(self, text, dimension, rotate=1):
		text_im = Image.new('RGB', dimension, (255,255,255))
		text_draw = ImageDraw.Draw(text_im)
		text_draw.text((0,0), text, fill=(0,0,255))
		box = (0,0,dimension[0], dimension[1])
		text_reg = text_im.crop(box)
		if rotate:
			text_reg = text_reg.transpose(Image.ROTATE_90)	#90 is anti-clockwise
		return text_reg
	
	
	def draw_incidence_matrix(self, input_fname, output_fname):
		"""
		09-26-05
			xlength is no of datasets
			ylength is no of genes
		09-26-05
			fix sme bugs, frequency is already string, repr causes additional "'"
		09-26-05
			each cell's height reduced to 1, too big to use char_height,
			gene_id is not drawn, frequency drawn every 2000 genes
		"""
		sys.stderr.write("Drawing gene presence_vector...")
		ylength_output = os.popen('wc %s'%input_fname)
		ylength_output = ylength_output.read()
		ylength = int(ylength_output.split()[0])
		
		xlength_output = os.popen('%s %s'%(os.path.expanduser('~/script/shell/count_columns.py'), input_fname))
		xlength_output = xlength_output.read()
		xlength = int(xlength_output.split()[-1])-2	#the first column is frequency and last column is gene_id
		
		char_width, char_height = self.get_char_dimension()
		dataset_no_dimension = (char_width*len(repr(xlength)), char_height)
		x_offset0 = 0
		x_offset1 = x_offset0 + dataset_no_dimension[0]
		y_offset0 = 0
		y_offset1 = y_offset0 + dataset_no_dimension[0]
		
		whole_dimension = (x_offset1+dataset_no_dimension[1]*xlength, y_offset1+ylength)
		
		reader = csv.reader(open(input_fname,'r'), delimiter='\t')
		im = Image.new('RGB',whole_dimension,(255,255,255))
		draw = ImageDraw.Draw(im)
		
		for i in range(xlength):	#write the text to a region and rotate anti-clockwise and paste it back
			text_region = self.get_text_region(repr(i+1), dataset_no_dimension)	#rotate
			box = (x_offset1+i*dataset_no_dimension[1], y_offset0, x_offset1+(i+1)*dataset_no_dimension[1], y_offset1)
			im.paste(text_region, box)
		
		gene_id_max_length = 10
		counter = 0
		for row in reader:
			frequency = row[0]
			presence_vector = row[1:-1]
			gene_id = row[-1]	#not used
			if len(gene_id)>gene_id_max_length:
				gene_id = gene_id[:gene_id_max_length]
			
			y_offset_upper = y_offset1+counter
			y_offset_lower = y_offset1+counter
			#draw the presence_vector
			for i in range(len(presence_vector)):
				if presence_vector[i]=='1':
					draw.rectangle((x_offset1+i*dataset_no_dimension[1], y_offset_upper, \
						x_offset1+(i+1)*dataset_no_dimension[1], y_offset_lower), fill=(0,255,0))
			if counter%2000 == 0:
				#draw the frequency every 2000 genes
				text_region = self.get_text_region(frequency, dataset_no_dimension, rotate=0)	#no rotate
				box = (x_offset0, y_offset_upper, x_offset1, y_offset_upper+dataset_no_dimension[1])
				im.paste(text_region, box)
			counter += 1
		del reader
		im.save(output_fname)
		del im
		
		sys.stderr.write('Done.\n')
	
	def run(self):
		"""
		09-26-05
			--write_gene_incidence_matrix()
				--expand_gene_id2presence_vector()
			--draw_incidence_matrix()
				--get_char_dimension()
				--get_text_region()
		"""
		self.write_gene_incidence_matrix(self.input_dir, self.output_fname)
		self.draw_incidence_matrix(self.output_fname, self.jpg_ofname)
	
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema="]
	opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:o:p:br", long_options_list)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	input_dir = None
	output_fname = None
	jpg_ofname = None
	debug = 0
	report = 0

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
		elif opt in ("-i"):
			input_dir = arg
		elif opt in ("-o"):
			output_fname = arg
		elif opt in ("-p"):
			jpg_ofname = arg
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r"):
			report = 1
		
	if input_dir and output_fname and jpg_ofname:
		instance = DrawGenePresenceMatrix(hostname, dbname, schema, input_dir, output_fname, \
			jpg_ofname, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
