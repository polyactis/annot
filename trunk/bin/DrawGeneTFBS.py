#!/usr/bin/env python
"""
Usage: DrawGeneTFBS.py -i -p [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-i ...,	the input gene_id
	-p ...,	output prefix
	-b	enable debugging
	-r	report the progress(a number)
	-h, --help	show this help

Examples:
	DrawGeneTFBS.py -i 11130 -p tf_bs_11130

Description:
	A program to draw TFBSs of a gene. Two pictures would be generated:
		1. tf_legend picture.(could ignore it)
		2. tf_line picture.
	data is from transfac.binding_site and transfac.prom_seq
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import getopt, math
from codense.common import get_char_dimension, get_text_region, db_connect
from sets import Set
import Image, ImageDraw


class DrawGeneTFBS:
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None,\
		input_gene_id=None, output_prefix=None, debug=0, report=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_gene_id = int(input_gene_id)
		self.output_prefix = output_prefix
		self.debug = int(debug)
		self.report = int(report)

	def get_gene_binding_sites(self, curs, gene_id, binding_site_table='transfac.binding_site', \
		prom_seq_table='transfac.prom_seq'):
		sys.stderr.write("Getting binding_sites for %s"%gene_id)
		curs.execute("select b.mt_id, b.strand, b.bs_disp_start, b.bs_disp_end from %s b, %s p \
			where p.id=b.prom_id  and p.prom_acc=%s and p.prom_type_id=1"%(binding_site_table, prom_seq_table, gene_id))
			#only upstream
		rows = curs.fetchall()
		sys.stderr.write(".\n")
		return rows
	
	def get_tf_name2binding_sites(self, tf_info_list):
		sys.stderr.write("Getting tf_name2binding_sites")
		tf_name2binding_sites = {}
		for tf_info_row in tf_info_list:
			tf_name = tf_info_row[0]
			if tf_name not in tf_name2binding_sites:
				tf_name2binding_sites[tf_name] = []
			tf_name2binding_sites[tf_name].append(tf_info_row)
		sys.stderr.write(".\n")
		return tf_name2binding_sites
	
	def get_gene_prom_seq_info(self, curs, gene_id, prom_seq_table='transfac.prom_seq'):
		sys.stderr.write("Getting gene prom_seq info for %s"%gene_id)
		curs.execute("SELECT prom_acc, length(sequence), strand from %s \
			where prom_acc=%s and prom_type_id=1"%(prom_seq_table, gene_id))
		rows = curs.fetchall()
		if rows:
			sys.stderr.write(".\n")
			return rows[0]
		else:
			sys.stderr.write(" NOT FOUND.\n")
			sys.exit(2)
	
	def draw_tf_legend(self, tf_name_list, padding_width=5, max_tf_name_length = 20):
		sys.stderr.write("Drawing tf_legend")
		#transform tf_name_list into a distinct list
		tf_name_set = Set(tf_name_list)
		tf_name_list = list(tf_name_set)
		tf_name_list.sort()
		
		#how many cuts for 0 to 255 (no_of_partitions-1)^3 <= len(tf_name_list) <= no_of_partitions^3
		no_of_partitions = int(math.ceil(math.pow(len(tf_name_list), 1/3.0)))
		step = 256/no_of_partitions
		#RGB three channels
		color_list = []
		for i in range(no_of_partitions):
			for j in range(no_of_partitions):
				for k in range(no_of_partitions):
					color_list.append((step*i, step*j, step*k))
		
		char_width, char_height = get_char_dimension()
		tf_name_dimension = (char_width*max_tf_name_length, char_height)
		legend_dimension = (5, char_height)
		x_offset0 = 0
		x_offset1 = x_offset0 + padding_width	#tf_name begin
		x_offset2 = x_offset1 + tf_name_dimension[0]	#tf_name end
		x_offset3 = x_offset2 + padding_width	#legend begin
		x_offset4 = x_offset3 + legend_dimension[0]	#legend end
		whole_dimension = (x_offset4 + padding_width, char_height*2*len(tf_name_list))
		tf_name2color = {}
		
		im = Image.new('RGB',(whole_dimension[0],whole_dimension[1]),(255,255,255))
		draw = ImageDraw.Draw(im)
		
		for i in range(len(tf_name_list)):
			tf_name = tf_name_list[i]
			tf_name2color[tf_name] = color_list[i]
			if len(tf_name) >max_tf_name_length:
				tf_name = tf_name[:max_tf_name_length]
			text_region = get_text_region(tf_name, tf_name_dimension, rotate=0)	#no rotate
			y_offset_upper = i*2*tf_name_dimension[1]
			y_offset_bottom = (i*2+1)*tf_name_dimension[1]
			box = (x_offset1, y_offset_upper, x_offset2, y_offset_bottom)
			im.paste(text_region, box)
			draw.rectangle((x_offset3, y_offset_upper, x_offset4, y_offset_bottom), fill=color_list[i])
		sys.stderr.write(".\n")
		return im, tf_name2color
	
	def draw_tf_line(self, seq_name, seq_length, seq_strand, tf_info_list, tf_name2color, padding_width=5,\
		padding_height=10, max_seq_name_length = 20, bs_width=1, im_10kb_length=800):
		
		sys.stderr.write("\t Drawing tf_line")
		char_width, char_height = get_char_dimension()
		#the idea is 10kb is treated as 800 pixels
		ratio = im_10kb_length/10000.0
		im_seq_length = int(ratio*seq_length)
		
		seq_dimension = (im_seq_length, char_height)
		seq_length_dimension = ((len(repr(seq_length))+3)*char_width, char_height)
		seq_name_dimension = (max_seq_name_length*char_width, char_height)
		x_offset0 = 0
		x_offset1 = x_offset0 + padding_width #seq starts
		x_offset2 = x_offset1 + seq_dimension[0] # seq ends
		x_offset3 = x_offset2 + padding_width #seq_length starts
		x_offset4 = x_offset3 + seq_length_dimension[0]	#seq_length ends
		x_offset5 = x_offset4 + padding_width	#seq_name starts
		x_offset6 = x_offset5 + seq_name_dimension[0]	#seq_name ends
		y_offset0 = 0
		y_offset1 = y_offset0 + padding_height	#the row starts
		y_offset2 = y_offset1 + seq_dimension[1]	#the row ends
		whole_dimension = (x_offset6 + padding_width, y_offset2 + padding_height)
		
		#1st draw the seq line,
		seq_im = Image.new('RGB', seq_dimension, (255, 255, 255))
		seq_draw = ImageDraw.Draw(seq_im)
		#draw the seq baseline in the middle
		seq_draw.rectangle((0, seq_dimension[1]/2, seq_dimension[0], seq_dimension[1]/2+1), fill=(0,0,0))
		for tf_info_row in tf_info_list:
			tf_name, bs_strand, bs_disp_start, bs_disp_end = tf_info_row
			im_bs_disp_start = int(ratio*bs_disp_start)
			seq_draw.rectangle((im_bs_disp_start, 0, im_bs_disp_start+bs_width, seq_dimension[1]), fill=tf_name2color[tf_name])
		if seq_strand == '-':	#to rotate 180 if seq_strand == '-'
			seq_im = seq_im.rotate(180)
		box = (0,0,seq_dimension[0], seq_dimension[1])
		seq_reg = seq_im.crop(box)
		
		#2nd paste seq line into main image
		im = Image.new('RGB',(whole_dimension[0],whole_dimension[1]),(255,255,255))
		draw = ImageDraw.Draw(im)
		box = (x_offset1, y_offset1, x_offset2, y_offset2)
		im.paste(seq_reg, box)
		
		#3rd draw seq_length
		text_region = get_text_region(str(seq_length)+'(%s)'%seq_strand, seq_length_dimension, rotate=0)	#no rotate
		box = (x_offset3, y_offset1, x_offset4, y_offset2)
		im.paste(text_region, box)
		
		#4gh draw seq_name
		seq_name = str(seq_name)
		if len(seq_name)>max_seq_name_length:
			seq_name = seq_name[:max_seq_name_length]
		text_region = get_text_region(seq_name, seq_name_dimension, rotate=0)	#no rotate
		box = (x_offset5, y_offset1, x_offset6, y_offset2)
		im.paste(text_region, box)
		sys.stderr.write(".\n")
		return im
	
	def get_composite_and_individual_tf_line(self, composite_tf_im, tf_name2binding_sites, \
		seq_name, seq_length, seq_strand, tf_name2color, padding_width,\
		padding_height, max_seq_name_length, bs_width, im_10kb_length):
		
		sys.stderr.write("Drawing composite_and_individual_tf_line\n")
		char_width, char_height = get_char_dimension()
		tf_name_dimension = (char_width*max_seq_name_length, char_height)
		tf_line_dimension = composite_tf_im.size
		
		x_offset0 = 0
		x_offset1 = x_offset0 + padding_width	#tf_name starts
		x_offset2 = x_offset1 + tf_name_dimension[0]	#tf_name ends
		x_offset3 = x_offset2 + padding_width	#tf_line starts
		x_offset4 = x_offset3 + tf_line_dimension[0]	#tf_line ends
		
		whole_dimension = (x_offset4, (len(tf_name2binding_sites)+1)*tf_line_dimension[1])
		im = Image.new('RGB',(whole_dimension[0],whole_dimension[1]),(255,255,255))
		draw = ImageDraw.Draw(im)
		
		#1st paste the composite_tf_im
		box = (x_offset3, 0, x_offset4, tf_line_dimension[1])
		im.paste(composite_tf_im, box)
		
		counter = 0
		tf_name_list = tf_name2binding_sites.keys()
		tf_name_list.sort()
		for tf_name in tf_name_list:
			tf_info_list = tf_name2binding_sites[tf_name]
			counter += 1	#starting from 1 because composite_tf_im is 0
			if len(tf_name) >max_seq_name_length:
				tf_name = tf_name[:max_seq_name_length]
			text_region = get_text_region(tf_name, tf_name_dimension, rotate=0)	#no rotate
			y_offset_upper = counter*tf_line_dimension[1]
			y_offset_bottom = (counter+1)*tf_line_dimension[1]
			box = (x_offset1, y_offset_upper+padding_height, x_offset2, y_offset_bottom-padding_height)
			im.paste(text_region, box)
			tf_line_im = self.draw_tf_line(seq_name, seq_length, seq_strand, tf_info_list, tf_name2color, \
				padding_width, padding_height, max_seq_name_length, bs_width, im_10kb_length)
			box = (x_offset3, y_offset_upper, x_offset4, y_offset_bottom)
			im.paste(tf_line_im, box)
		sys.stderr.write(".\n")
		return im
	
	def run(self):
		"""
		01-03-06
		"""
		padding_width=5
		padding_height=10
		max_seq_name_length = 20
		bs_width=1
		im_10kb_length=800
		
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		tf_info_list = self.get_gene_binding_sites(curs, self.input_gene_id)
		tf_name2binding_sites = self.get_tf_name2binding_sites(tf_info_list)
		seq_name, seq_length, seq_strand = self.get_gene_prom_seq_info(curs, self.input_gene_id)
		
		tf_name_list = tf_name2binding_sites.keys()
		tf_legend_im, tf_name2color = self.draw_tf_legend(tf_name_list, padding_width, max_seq_name_length)
		composite_tf_im = self.draw_tf_line(seq_name, seq_length, seq_strand, tf_info_list, tf_name2color, \
			padding_width, padding_height, max_seq_name_length, bs_width, im_10kb_length)
		
		im = self.get_composite_and_individual_tf_line(composite_tf_im, tf_name2binding_sites, \
			seq_name, seq_length, seq_strand, tf_name2color, padding_width,\
			padding_height, max_seq_name_length, bs_width, im_10kb_length)
		
		tf_legend_output_fname = '%s_tf_legend.png'%self.output_prefix
		tf_legend_im.save(tf_legend_output_fname)
		tf_line_output_fname = '%s_tf_line.png'%self.output_prefix
		im.save(tf_line_output_fname)
		

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema="]
	opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:p:br", long_options_list)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = None
	input_gene_id = None
	output_prefix = None
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
		elif opt in ("-i",):
			input_gene_id = arg
		elif opt in ("-p",):
			output_prefix = arg
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r",):
			report = 1
		
	if input_gene_id and output_prefix:
		instance = DrawGeneTFBS(hostname, dbname, schema, input_gene_id, output_prefix, \
			debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
