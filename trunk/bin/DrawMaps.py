#!/usr/bin/env python
"""
Usage: DrawMaps.py -k -i -p -s -o [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-p ...,	pattern_table
	-s ...,	cluster_bs_table
	-i ...	the inputfname, filtered prediction file(output of parse_haifeng_markov())
	-l ...,	lm_bit of setting1, ('00001', default)
	-a ...	acc_cutoff of setting1, (0.6 default)
	-o ...	output prefix
	-y ...	type of -i(inputfname), 1(schema prefix, for p_gene_table and gene_p_table, default), 2(filtered prediction file)
	-x ...	the font size, 20 (default)
	-f ...	the font path, /usr/share/fonts/truetype/msttcorefonts/arial.ttf (default)
	-b	enable debugging, no debug by default
	-r	report the progress(a number)
	-h, --help              show this help

Examples:
	DrawMaps.py -k hs_fim_65 -i ~/script/haifeng_annot/8-pred/human.n2.s200.pred.1.2.max 
	 -p pattern_hs_fim_65_m5x65s4l5 -s cluster_bs_hs_fim_65_m5x65s4l5e0p001geneid
	 -o /tmpcmb-01/yuhuang/map
	
	DrawMaps.py -k hs_fim_65 -i hs_fim_65_n2s175_m5x65s4l5_ft2 -p pattern_hs_fim_65_n2s175_m5x65s4l5 
		-s bs_hs_fim_65_n2s175_m5x65s4l5_ft2_e5_000001a60p01y1 
		-o ~/tmp/hs_fim_65_n2s175_m5x65s4l5_ft2_map -l 000001
Description:
	Draw function_map and gene_function_map and tf_map.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import os, sys, csv, getopt
from codense.common import get_char_dimension, get_text_region, db_connect,\
	form_schema_tables, p_gene_id_set_from_gene_p_table, get_go_no2name, \
	cluster_bs_id_set_from_good_bs_table, get_mt_no2tf_name, draw_grid, \
	get_gene_id2gene_symbol, get_go_id2name, get_go_no2go_id
from sets import Set
import Image, ImageDraw, ImageFont
from MpiFromDatasetSignatureToPattern import encodeOccurrenceBv, decodeOccurrenceToBv, decodeOccurrence

class DrawMaps:
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, pattern_table='', cluster_bs_table='',\
		inputfname=None, lm_bit='00001', acc_cutoff=0.6, output_prefix=None, type=1, font_size=20, \
		font_path='/usr/share/fonts/truetype/msttcorefonts/arial.ttf', debug=0, report=0):
		"""
		2006-09-26
		2006-11-06
			add type
		2006-12-13
			add font_size and font_path
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.pattern_table = pattern_table
		self.cluster_bs_table = cluster_bs_table
		self.inputfname = inputfname
		self.lm_bit = lm_bit
		self.acc_cutoff = float(acc_cutoff)
		self.output_prefix = output_prefix
		self.type = int(type)
		self.font_size = int(font_size)
		self.font_path = font_path
		self.function_name_length = 44	#truncate if exceed, no public interface
		self.debug = int(debug)
		self.report = int(report)
	
	
	def get_no_of_p_funcs_gene_no_go_no_list_from_file(self, inputfname):
		"""
		2006-09-26
			from a file
			go_no is actually go_id (GeneOntology acc, GO:0004283)
			mcl_id2go_no_set is for  get_recurrence_go_no_rec_array_cluster_id_ls()
		2006-11-06
			rename it from get_no_of_p_funcs_gene_no_go_no_list to get_no_of_p_funcs_gene_no_go_no_list_from_file
		"""
		sys.stderr.write("Getting no_of_p_funcs_gene_no_go_no_list...\n")
		gene_no2go_no_set = {}
		mcl_id2go_no_set = {}
		reader = csv.reader(open(inputfname), delimiter='\t')
		counter = 0
		for row in reader:
			gene_id, pattern_id, go_id, is_known, is_correct, markov_score, p_value = row
			gene_no = int(gene_id)
			pattern_id = int(pattern_id)
			if gene_no not in gene_no2go_no_set:
				gene_no2go_no_set[gene_no] = Set()
			gene_no2go_no_set[gene_no].add(go_id)
			if pattern_id not in mcl_id2go_no_set:
				mcl_id2go_no_set[pattern_id] = Set()
			mcl_id2go_no_set[pattern_id].add(go_id)
			counter += 1
		if self.report:
			sys.stderr.write("%s%s\n"%('\x08'*20, counter))
		#transform the dict into list
		no_of_p_funcs_gene_no_go_no_list = []
		for gene_no, go_no_set in gene_no2go_no_set.iteritems():
			no_of_p_funcs_gene_no_go_no_list.append([len(go_no_set), gene_no, go_no_set])
		no_of_p_funcs_gene_no_go_no_list.sort()
		sys.stderr.write("End getting no_of_p_funcs_gene_no_go_no_list.\n")
		return no_of_p_funcs_gene_no_go_no_list, mcl_id2go_no_set
	
	def get_no_of_p_funcs_gene_no_go_no_list_from_db(self, curs, p_gene_table, given_p_gene_set, go_no2go_id):
		"""
		2006-11-06
			modified after the old version of get_no_of_p_funcs_gene_no_go_no_list
		"""
		sys.stderr.write("Getting no_of_p_funcs_gene_no_go_no_list...\n")
		gene_no2go_no_set = {}
		mcl_id2go_no_set = {}
		curs.execute("DECLARE crs CURSOR FOR SELECT p_gene_id, gene_no, go_no, mcl_id \
			from %s"%p_gene_table)
		curs.execute("fetch 10000 from crs")
		rows = curs.fetchall()
		counter = 0
		real_counter = 0
		while rows:
			for row in rows:
				p_gene_id, gene_no, go_no, mcl_id = row
				if p_gene_id in given_p_gene_set:
					if gene_no not in gene_no2go_no_set:
						gene_no2go_no_set[gene_no] = Set()
					gene_no2go_no_set[gene_no].add(go_no2go_id[go_no])
					if mcl_id not in mcl_id2go_no_set:
						mcl_id2go_no_set[mcl_id] = Set()
					mcl_id2go_no_set[mcl_id].add(go_no2go_id[go_no])
					real_counter += 1
			counter += 1
			if self.report:
				sys.stderr.write("%s%s/%s"%('\x08'*20, counter, real_counter))
			curs.execute("fetch 10000 from crs")
			rows = curs.fetchall()
		curs.execute("close crs")
		#transform the dict into list
		no_of_p_funcs_gene_no_go_no_list = []
		for gene_no, go_no_set in gene_no2go_no_set.iteritems():
			no_of_p_funcs_gene_no_go_no_list.append([len(go_no_set), gene_no, go_no_set])
		no_of_p_funcs_gene_no_go_no_list.sort()
		sys.stderr.write("End getting no_of_p_funcs_gene_no_go_no_list.\n")
		return no_of_p_funcs_gene_no_go_no_list, mcl_id2go_no_set
	
	def draw_gene_function_map(self, no_of_p_funcs_gene_no_go_no_list, go_no2index, function_name_region,\
		output_fname, function_name_length, char_dimension, no_of_functions, font):
		"""
		10-31-05
		11-28-05 draw_grid
		2006-12-13 add font
		"""
		sys.stderr.write("Drawing gene_function_map...\n")
		char_width, char_height = char_dimension
		
		gene_no_repr_length = 0
		for row in no_of_p_funcs_gene_no_go_no_list:
			if len(repr(row[1]))>gene_no_repr_length:
				gene_no_repr_length = len(repr(row[1]))
		gene_no_repr_dimension = (char_width*gene_no_repr_length, char_height)
		function_name_dimension = (char_width*function_name_length, char_height)	#will rotate
		no_of_p_funcs_repr_length = len(repr(no_of_p_funcs_gene_no_go_no_list[-1][0]))	#-1 is the largest number
		no_of_p_funcs_repr_dimension = (char_width*no_of_p_funcs_repr_length, char_height)
		
		x_offset0 = 0
		x_offset1 = gene_no_repr_dimension[0]
		x_offset2 = x_offset1 + no_of_functions*function_name_dimension[1]
		y_offset0 = 0
		y_offset1 = function_name_dimension[0]
		
		whole_dimension = (x_offset2+no_of_p_funcs_repr_dimension[0], \
			y_offset1+len(no_of_p_funcs_gene_no_go_no_list)*gene_no_repr_dimension[1])
		
		im = Image.new('RGB',(whole_dimension[0],whole_dimension[1]),(255,255,255))
		draw = ImageDraw.Draw(im)
		
		box = (x_offset1, y_offset0, x_offset2, y_offset1)
		im.paste(function_name_region, box)
		
		for i in range(len(no_of_p_funcs_gene_no_go_no_list)):
			no_of_p_funcs, gene_no, go_no_set = no_of_p_funcs_gene_no_go_no_list[i]
			y_offset_upper = y_offset1 + i*gene_no_repr_dimension[1]
			y_offset_lower = y_offset_upper + gene_no_repr_dimension[1]
			#draw gene_no
			text_region = get_text_region(repr(gene_no), gene_no_repr_dimension, rotate=0, font=font)
			box = (x_offset0, y_offset_upper, x_offset1, y_offset_lower)
			im.paste(text_region, box)
			#draw p_function
			for go_no in go_no_set:
				index = go_no2index[go_no]
				x_offset_left = x_offset1 + index*function_name_dimension[1]
				x_offset_right = x_offset_left + function_name_dimension[1]
				draw.rectangle((x_offset_left, y_offset_upper, x_offset_right, y_offset_lower), fill=(0,255,0))
			#draw no_of_p_funcs
			text_region = get_text_region(repr(no_of_p_funcs), no_of_p_funcs_repr_dimension, rotate=0, font=font)
			box = (x_offset2, y_offset_upper, whole_dimension[0], y_offset_lower)
			im.paste(text_region, box)
		#11-28-05
		draw_grid(im, draw, [x_offset1, y_offset1, x_offset2, whole_dimension[1]], char_height, char_height)
		
		im = im.rotate(270)
		im.save(output_fname)
		del im
		sys.stderr.write("End drawing gene_function_map.\n")
	
	
	def get_recurrence_go_no_rec_array_cluster_id_ls(self, curs, pattern_table, mcl_id2go_no_set):
		"""
		2006-09-26
			from pattern_table and use mcl_id2go_no_set
			go_no_list is the go_id Set
			mcl_id2enc_recurrence is for get_recurrence_rec_array_bs_no_list()
		"""
		sys.stderr.write("Getting recurrence_go_no_rec_array_cluster_id_ls...\n")
		no_of_datasets = 0
		go_no2recurrence_cluster_id = {}
		mcl_id2enc_recurrence = {}
		curs.execute("DECLARE crs CURSOR FOR SELECT id, recurrence_array from %s"\
			%pattern_table)
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		counter = 0
		real_counter = 0
		while rows:
			for row in rows:
				mcl_id, recurrence_array = row
				if mcl_id in mcl_id2go_no_set:
					#if this pattern has functions predicted
					recurrence_array = recurrence_array[1:-1].split(',')
					recurrence_array = map(float, recurrence_array)	#this is not a binary 0/1 array
					occurrence_cutoff_func = lambda x: int(x>=0.8)	#0.8 is arbitrary
					recurrence_array = map(occurrence_cutoff_func, recurrence_array)
					if no_of_datasets == 0:
						no_of_datasets = len(recurrence_array)
					go_no_list = mcl_id2go_no_set[mcl_id]
					encoded_recurrence = encodeOccurrenceBv(recurrence_array)
					mcl_id2enc_recurrence[mcl_id] = encoded_recurrence	#2006-09-26
					for go_no in go_no_list:
						if go_no not in go_no2recurrence_cluster_id:
							go_no2recurrence_cluster_id[go_no] = [encoded_recurrence, Set([mcl_id])]
								#use Set() because mcl_id has duplicates due to different p-values
						else:
							go_no2recurrence_cluster_id[go_no][0] = \
								go_no2recurrence_cluster_id[go_no][0] | encoded_recurrence
							go_no2recurrence_cluster_id[go_no][1].add(mcl_id)
					real_counter += 1
				counter += 1
			if self.report:
				sys.stderr.write("%s%s\t%s"%('\x08'*20, counter, real_counter))
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		curs.execute("close crs")
		
		recurrence_go_no_rec_array_cluster_id_ls = []
		for go_no in go_no2recurrence_cluster_id:
			encoded_recurrence, mcl_id_set = go_no2recurrence_cluster_id[go_no]
			recurrence_array = decodeOccurrence(encoded_recurrence)	#not binary vector
			recurrence = len(recurrence_array)
			recurrence_go_no_rec_array_cluster_id_ls.append([recurrence, go_no, recurrence_array, mcl_id_set])
		
		recurrence_go_no_rec_array_cluster_id_ls.sort()
		sys.stderr.write("End getting recurrence_go_no_rec_array_cluster_id_ls.\n")
		return recurrence_go_no_rec_array_cluster_id_ls, no_of_datasets, mcl_id2enc_recurrence
	
	def draw_function_map(self, recurrence_go_no_rec_array_cluster_id_ls, no_of_datasets, go_no2name, \
		output_fname, function_name_length, char_dimension, no_of_functions, font):
		"""
		10-31-05 from misc.py. make it independent and return something for draw_gene_function_map()
		11-28-05 add grid
		11-28-05 include go_no in go_name
		2006-12-13 add font
		"""
		sys.stderr.write("Drawing function_map...\n")
		
		char_width, char_height = char_dimension
		dataset_no_length = len(repr(no_of_datasets))
		dataset_no_dimension = (char_width*dataset_no_length, char_height)	#one is not rotated, one is rotated
		no_of_clusters_dimension = (char_width*7, char_height)	#will rotate
		function_name_dimension = (char_width*function_name_length, char_height)	#will rotate
		
		x_offset0 = 0
		x_offset1 = dataset_no_dimension[0]
		y_offset0 = 0
		y_offset1 = function_name_dimension[0]
		y_offset2 = y_offset1 + no_of_datasets*dataset_no_dimension[1]
		y_offset3 = y_offset2 + dataset_no_dimension[0]
		whole_dimension = (x_offset1+no_of_functions*char_height, \
			y_offset3+no_of_clusters_dimension[0])
		
		im = Image.new('RGB',(whole_dimension[0],whole_dimension[1]),(255,255,255))
		draw = ImageDraw.Draw(im)
		#dataset_no section
		for i in range(no_of_datasets):
			text_region = get_text_region(repr(i+1), dataset_no_dimension, rotate=0, font=font)	#no rotate
			box = (x_offset0, y_offset1+i*dataset_no_dimension[1], x_offset1, y_offset1+(i+1)*dataset_no_dimension[1])
			im.paste(text_region, box)
		#10-31-05 following for draw_gene_function_map() to ensure correspondence.
		go_no2index = {}
		function_name_box = (x_offset1,y_offset0, whole_dimension[0], y_offset1)
		
		for i in range(len(recurrence_go_no_rec_array_cluster_id_ls)):
			recurrence, go_no, recurrence_array, mcl_id_set = recurrence_go_no_rec_array_cluster_id_ls[i]
			go_no2index[go_no] = i	#10-31-05
			x_offset_left = x_offset1+i*function_name_dimension[1]
			x_offset_right = x_offset1+(i+1)*function_name_dimension[1]
			#function_name
			go_name = '%s %s'%(go_no, go_no2name[go_no])	#11-28-05 include go_no
			if len(go_name)>function_name_length:
				go_name = go_name[:function_name_length]
			text_region = get_text_region(go_name, function_name_dimension, font=font)	#rotate
			box = (x_offset_left, y_offset0, x_offset_right, y_offset1)
			im.paste(text_region, box)
			
			#fill in a cell for each dataset_no
			for dataset_no in recurrence_array:
				draw.rectangle((x_offset_left, y_offset1+(dataset_no-1)*dataset_no_dimension[1], \
					x_offset_right, y_offset1+dataset_no*dataset_no_dimension[1]), fill=(0,255,0))
			# write down the recurrence
			text_region = get_text_region(repr(recurrence), dataset_no_dimension, font=font)	#rotate
			box = (x_offset_left, y_offset2, x_offset_right, y_offset3)
			im.paste(text_region, box)
			#write down the no_of_clusters
			text_region = get_text_region(repr(len(mcl_id_set)), no_of_clusters_dimension, font=font)	#rotate
			box = (x_offset_left, y_offset3, x_offset_right, whole_dimension[1])
			im.paste(text_region, box)
		function_name_region = im.crop(function_name_box)	#10-31-05
		
		#11-28-05
		draw_grid(im, draw, [x_offset1, y_offset1, whole_dimension[0], y_offset2], char_height, char_height)
		
		im = im.rotate(270)
		im.save(output_fname)
		del im
		sys.stderr.write("End drawing function_map.\n")
		return go_no2index, function_name_region
	
	def get_recurrence_rec_array_bs_no_list(self, curs, cluster_bs_table, mcl_id2enc_recurrence):
		"""
		11-01-05
		""" 
		sys.stderr.write("Getting recurrence_rec_array_bs_no_list...\n")
		bs_no2enc_recurrence = {}
		curs.execute("DECLARE crs CURSOR FOR select c.mcl_id, c.bs_no_list from %s c"%(cluster_bs_table))
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		counter = 0
		real_counter =0
		while rows:
			for row in rows:
				mcl_id, bs_no_list = row
				if mcl_id in mcl_id2enc_recurrence:
					encoded_recurrence = mcl_id2enc_recurrence[mcl_id]
					bs_no_list = bs_no_list[1:-1].split(',')
					bs_no_list = map(int, bs_no_list)
					for bs_no in bs_no_list:
						if bs_no not in bs_no2enc_recurrence:
							bs_no2enc_recurrence[bs_no] = encoded_recurrence
						else:
							bs_no2enc_recurrence[bs_no] |= encoded_recurrence
					real_counter += 1
				counter += 1
			if self.report:
				sys.stderr.write("%s%s/%s"%('\x08'*20, counter, real_counter))
			curs.execute("fetch 10000 from crs")
			rows = curs.fetchall()
		curs.execute("close crs")
		recurrence_rec_array_bs_no_list = []
		for bs_no, enc_recurrence in bs_no2enc_recurrence.iteritems():
			recurrence_array = decodeOccurrence(enc_recurrence)	#not binary vector
			recurrence = len(recurrence_array)
			recurrence_rec_array_bs_no_list.append([recurrence, recurrence_array, bs_no])
		recurrence_rec_array_bs_no_list.sort()
		sys.stderr.write("End getting recurrence_rec_array_bs_no_list.\n")
		return recurrence_rec_array_bs_no_list
	
	def draw_tf_map(self, recurrence_rec_array_bs_no_list, no_of_datasets, mt_no2tf_name, \
		output_fname, function_name_length, char_dimension, font):
		"""
		11-28-05 add grid
		2006-12-13 add font
		"""
		sys.stderr.write("Drawing tf_map...\n")
		
		no_of_tfs = len(recurrence_rec_array_bs_no_list)
		char_width, char_height = char_dimension
		dataset_no_length = len(repr(no_of_datasets))
		dataset_no_dimension = (char_width*dataset_no_length, char_height)
			#one is not rotated(dataset_no), one is rotated(recurrence)
		function_name_dimension = (char_width*function_name_length, char_height)	#will rotate
		
		x_offset0 = 0
		x_offset1 = dataset_no_dimension[0]
		y_offset0 = 0
		y_offset1 = function_name_dimension[0]
		y_offset2 = y_offset1 + no_of_datasets*dataset_no_dimension[1]
		y_offset3 = y_offset2 + dataset_no_dimension[0]
		whole_dimension = (x_offset1+no_of_tfs*char_height, y_offset3)
		
		im = Image.new('RGB',(whole_dimension[0],whole_dimension[1]),(255,255,255))
		draw = ImageDraw.Draw(im)
		#dataset_no section
		for i in range(no_of_datasets):
			text_region = get_text_region(repr(i+1), dataset_no_dimension, rotate=0, font=font)	#no rotate
			box = (x_offset0, y_offset1+i*dataset_no_dimension[1], x_offset1, y_offset1+(i+1)*dataset_no_dimension[1])
			im.paste(text_region, box)
			
		for i in range(len(recurrence_rec_array_bs_no_list)):
			recurrence, recurrence_array, bs_no = recurrence_rec_array_bs_no_list[i]
			x_offset_left = x_offset1+i*function_name_dimension[1]
			x_offset_right = x_offset1+(i+1)*function_name_dimension[1]
			#function_name
			tf_name = '%s %s'%(bs_no, mt_no2tf_name[bs_no])
			if len(tf_name)>function_name_length:
				tf_name = tf_name[:function_name_length]
			text_region = get_text_region(tf_name, function_name_dimension, font=font)	#rotate
			box = (x_offset_left, y_offset0, x_offset_right, y_offset1)
			im.paste(text_region, box)
			
			#fill in a cell for each dataset_no
			for dataset_no in recurrence_array:
				draw.rectangle((x_offset_left, y_offset1+(dataset_no-1)*dataset_no_dimension[1], \
					x_offset_right, y_offset1+dataset_no*dataset_no_dimension[1]), fill=(0,255,0))
			# write down the recurrence
			text_region = get_text_region(repr(recurrence), dataset_no_dimension, font=font)	#rotate
			box = (x_offset_left, y_offset2, x_offset_right, y_offset3)
			im.paste(text_region, box)
		#11-28-05
		draw_grid(im, draw, [x_offset1, y_offset1, whole_dimension[0], y_offset2], char_height, char_height)
		
		im = im.rotate(270)
		im.save(output_fname)
		del im
		sys.stderr.write("End drawing tf_map.\n")
	
	def run(self):
		"""
		10-31-05
		2006-09-26
			modify it to be compatible with the modified pipeline from haifeng
		2006-11-06
			add type
		2006-12-13
			use font_path and font_size
			
			--form_schema_tables()
			--db_connect()
			--get_char_dimension()
			
			--get_no_of_p_funcs_gene_no_go_no_list()
			--get_recurrence_go_no_rec_array_cluster_id_ls()
			--get_go_no2name()
			--draw_function_map()
			
			--draw_gene_function_map()

			--get_recurrence_rec_array_bs_no_list()
			--get_mt_no2tf_name()
			--draw_tf_map()
		"""
		schema_instance = form_schema_tables(self.inputfname, self.acc_cutoff, self.lm_bit)
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		font = ImageFont.truetype(self.font_path, self.font_size)
		char_dimension = font.getsize('a')
		#char_dimension = get_char_dimension()
		
		#go_no2name = get_go_no2name(curs)
		go_no2name = get_go_id2name(curs)
		if self.type==1:
			go_no2go_id = get_go_no2go_id(curs)
			given_p_gene_set = p_gene_id_set_from_gene_p_table(curs, schema_instance.gene_p_table)
			no_of_p_funcs_gene_no_go_no_list, mcl_id2go_no_set = self.get_no_of_p_funcs_gene_no_go_no_list_from_db(curs, \
				schema_instance.p_gene_table, given_p_gene_set, go_no2go_id)
		elif self.type==2:
			no_of_p_funcs_gene_no_go_no_list, mcl_id2go_no_set = self.get_no_of_p_funcs_gene_no_go_no_list_from_file(self.inputfname)
		
		
		recurrence_go_no_rec_array_cluster_id_ls, no_of_datasets, mcl_id2enc_recurrence = \
			self.get_recurrence_go_no_rec_array_cluster_id_ls(curs, self.pattern_table, mcl_id2go_no_set)
		
		no_of_functions = len(recurrence_go_no_rec_array_cluster_id_ls)
		function_map_output_fname = '%s.function_map.png'%self.output_prefix
		go_no2index, function_name_region = self.draw_function_map(recurrence_go_no_rec_array_cluster_id_ls, no_of_datasets,\
			go_no2name, function_map_output_fname, self.function_name_length, char_dimension, no_of_functions, font)				
		
		gene_function_map_output_fname = '%s.gene_function_map.png'%self.output_prefix
		self.draw_gene_function_map(no_of_p_funcs_gene_no_go_no_list, go_no2index, function_name_region,\
			gene_function_map_output_fname, self.function_name_length, char_dimension, no_of_functions, font)
		
		
		#tf_map requires mcl_id2enc_recurrence and no_of_datasets from above
		recurrence_rec_array_bs_no_list = self.get_recurrence_rec_array_bs_no_list(curs, self.cluster_bs_table, mcl_id2enc_recurrence)
		mt_no2tf_name = get_gene_id2gene_symbol(curs, tax_id=9606)
		#mt_no2tf_name = get_mt_no2tf_name()
		tf_map_output_fname = '%s.tf_map.png'%self.output_prefix
		self.draw_tf_map(recurrence_rec_array_bs_no_list, no_of_datasets, mt_no2tf_name, \
			tf_map_output_fname, self.function_name_length, char_dimension, font)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema="]
	opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:p:s:i:l:a:o:y:x:f:br", long_options_list)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	pattern_table = ''
	cluster_bs_table = ''
	inputfname = None
	lm_bit = '00001'
	acc_cutoff = 0.6
	output_prefix = None
	type = 1
	font_size = 20
	font_path = '/usr/share/fonts/truetype/msttcorefonts/arial.ttf'
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
		elif opt in ("-p",):
			pattern_table = arg
		elif opt in ("-s",):
			cluster_bs_table = arg
		elif opt in ("-i", ):
			inputfname = arg
		elif opt in ("-l", ):
			lm_bit = arg
		elif opt in ("-a", ):
			acc_cutoff = float(arg)
		elif opt in ("-o",):
			output_prefix = arg
		elif opt in ("-y",):
			type = int(arg)
		elif opt in ("-x",):
			font_size = int(arg)
		elif opt in ("-f",):
			font_path = arg
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r", ):
			report = 1
		
	if pattern_table and cluster_bs_table and inputfname and output_prefix:
		instance = DrawMaps(hostname, dbname, schema, pattern_table, cluster_bs_table, inputfname, lm_bit, acc_cutoff, output_prefix, \
			type, font_size, font_path, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
