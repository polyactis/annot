#!/usr/bin/env python
"""
Usage: DrawMaps.py -k -i -o [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-i ...	the inputfname
	-l ...,	lm_bit of setting1, ('00001', default)
	-a ...	acc_cutoff of setting1, (0.6 default)
	-p ...	output prefix
	-b	enable debugging, no debug by default
	-r	report the progress(a number)
	-h, --help              show this help

Examples:
	DrawMaps.py -k rnfim33 -i rnfim33m3x33bfsdfg0p0 -p /tmp/rnfim33

Description:
	Draw function_map and gene_function_map.
"""

import os, sys, csv, getopt
sys.path += [os.path.join(os.path.expanduser('~/script/annot/bin'))]
from codense.common import get_char_dimension, get_text_region, db_connect,\
	form_schema_tables, p_gene_id_set_from_gene_p_table, get_go_no2name, \
	cluster_bs_id_set_from_good_bs_table, get_mt_no2tf_name
from sets import Set
import Image, ImageDraw
from MpiFromDatasetSignatureToPattern import encodeOccurrenceBv, decodeOccurrenceToBv, decodeOccurrence

class DrawMaps:
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None,\
		inputfname=None, lm_bit='00001', acc_cutoff=0.6, output_prefix=None, debug=0, report=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.inputfname = inputfname
		self.lm_bit = lm_bit
		self.acc_cutoff = float(acc_cutoff)
		self.output_prefix = output_prefix
		self.function_name_length = 40	#truncate if exceed, no public interface
		self.debug = int(debug)
		self.report = int(report)
	
	def get_no_of_p_funcs_gene_no_go_no_list(self, curs, p_gene_table, given_p_gene_set):
		sys.stderr.write("Getting no_of_p_funcs_gene_no_go_no_list...\n")
		gene_no2go_no_set = {}
		curs.execute("DECLARE crs CURSOR FOR SELECT p_gene_id, gene_no, go_no \
			from %s"%p_gene_table)
		curs.execute("fetch 10000 from crs")
		rows = curs.fetchall()
		counter = 0
		real_counter = 0
		while rows:
			for row in rows:
				p_gene_id, gene_no, go_no = row
				if p_gene_id in given_p_gene_set:
					if gene_no not in gene_no2go_no_set:
						gene_no2go_no_set[gene_no] = Set()
					gene_no2go_no_set[gene_no].add(go_no)
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
		return no_of_p_funcs_gene_no_go_no_list
		
	def draw_gene_function_map(self, no_of_p_funcs_gene_no_go_no_list, go_no2index, function_name_region,\
		output_fname, function_name_length, char_dimension, no_of_functions):
		"""
		10-31-05
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
			text_region = get_text_region(repr(gene_no), gene_no_repr_dimension, rotate=0)
			box = (x_offset0, y_offset_upper, x_offset1, y_offset_lower)
			im.paste(text_region, box)
			#draw p_function
			for go_no in go_no_set:
				index = go_no2index[go_no]
				x_offset_left = x_offset1 + index*function_name_dimension[1]
				x_offset_right = x_offset_left + function_name_dimension[1]
				draw.rectangle((x_offset_left, y_offset_upper, x_offset_right, y_offset_lower), fill=(0,255,0))
			#draw no_of_p_funcs
			text_region = get_text_region(repr(no_of_p_funcs), no_of_p_funcs_repr_dimension, rotate=0)
			box = (x_offset2, y_offset_upper, whole_dimension[0], y_offset_lower)
			im.paste(text_region, box)
		im = im.rotate(270)
		im.save(output_fname)
		del im
		sys.stderr.write("End drawing gene_function_map.\n")
	
	
	def get_recurrence_go_no_rec_array_cluster_id_ls(self, curs, good_cluster_table):
		"""
		"""
		sys.stderr.write("Getting recurrence_go_no_rec_array_cluster_id_ls...\n")
		no_of_datasets = 0
		go_no2recurrence_cluster_id = {}
		curs.execute("DECLARE crs CURSOR FOR SELECT mcl_id, recurrence_array, go_no_list from %s"\
			%good_cluster_table)
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		counter = 0
		while rows:
			for row in rows:
				mcl_id, recurrence_array, go_no_list = row
				recurrence_array = recurrence_array[1:-1].split(',')
				recurrence_array = map(int, recurrence_array)
				if no_of_datasets == 0:
					no_of_datasets = len(recurrence_array)
				go_no_list = go_no_list[1:-1].split(',')
				go_no_list = map(int, go_no_list)
				for go_no in go_no_list:
					if go_no not in go_no2recurrence_cluster_id:
						go_no2recurrence_cluster_id[go_no] = [encodeOccurrenceBv(recurrence_array), Set([mcl_id])]
							#use Set() because mcl_id has duplicates due to different p-values
					else:
						go_no2recurrence_cluster_id[go_no][0] = \
							go_no2recurrence_cluster_id[go_no][0] | encodeOccurrenceBv(recurrence_array)
						go_no2recurrence_cluster_id[go_no][1].add(mcl_id)
				counter += 1
			if self.report:
				sys.stderr.write("%s%s"%('\x08'*20, counter))
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
		return recurrence_go_no_rec_array_cluster_id_ls, no_of_datasets
	
	def draw_function_map(self, recurrence_go_no_rec_array_cluster_id_ls, no_of_datasets, go_no2name, \
		output_fname, function_name_length, char_dimension, no_of_functions):
		"""
		10-31-05 from misc.py. make it independent and return something for draw_gene_function_map()
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
			text_region = get_text_region(repr(i+1), dataset_no_dimension, rotate=0)	#no rotate
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
			go_name = go_no2name[go_no]
			if len(go_name)>function_name_length:
				go_name = go_name[:function_name_length]
			text_region = get_text_region(go_name, function_name_dimension)	#rotate
			box = (x_offset_left, y_offset0, x_offset_right, y_offset1)
			im.paste(text_region, box)
			
			#fill in a cell for each dataset_no
			for dataset_no in recurrence_array:
				draw.rectangle((x_offset_left, y_offset1+(dataset_no-1)*dataset_no_dimension[1], \
					x_offset_right, y_offset1+dataset_no*dataset_no_dimension[1]), fill=(0,255,0))
			# write down the recurrence
			text_region = get_text_region(repr(recurrence), dataset_no_dimension)	#rotate
			box = (x_offset_left, y_offset2, x_offset_right, y_offset3)
			im.paste(text_region, box)
			#write down the no_of_clusters
			text_region = get_text_region(repr(len(mcl_id_set)), no_of_clusters_dimension)	#rotate
			box = (x_offset_left, y_offset3, x_offset_right, whole_dimension[1])
			im.paste(text_region, box)
		function_name_region = im.crop(function_name_box)	#10-31-05
		im = im.rotate(270)
		im.save(output_fname)
		del im
		sys.stderr.write("End drawing function_map.\n")
		return go_no2index, function_name_region
	
	def get_recurrence_rec_array_bs_no_list(self, curs, schema_instance):
		"""
		11-01-05
		""" 
		cluster_bs_id_set = cluster_bs_id_set_from_good_bs_table(curs, schema_instance.good_bs_table)
		sys.stderr.write("Getting recurrence_rec_array_bs_no_list...\n")
		bs_no2enc_recurrence = {}
		no_of_datasets = 0
		curs.execute("DECLARE crs CURSOR FOR select g.recurrence_array, c.bs_no_list,c.id from %s g, %s c\
				where c.mcl_id=g.mcl_id"%(schema_instance.good_cluster_table,\
				schema_instance.cluster_bs_table))
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		counter = 0
		real_counter =0
		while rows:
			for row in rows:
				recurrence_array, bs_no_list, id = row
				if id in cluster_bs_id_set:
					recurrence_array = recurrence_array[1:-1].split(',')
					recurrence_array = map(int, recurrence_array)
					bs_no_list = bs_no_list[1:-1].split(',')
					bs_no_list = map(int, bs_no_list)
					if no_of_datasets == 0:
						no_of_datasets = len(recurrence_array)
					for bs_no in bs_no_list:
						if bs_no not in bs_no2enc_recurrence:
							bs_no2enc_recurrence[bs_no] = encodeOccurrenceBv(recurrence_array)
						else:
							bs_no2enc_recurrence[bs_no] |= encodeOccurrenceBv(recurrence_array)
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
		return recurrence_rec_array_bs_no_list, no_of_datasets
	
	def draw_tf_map(self, recurrence_rec_array_bs_no_list, no_of_datasets, mt_no2tf_name, \
		output_fname, function_name_length, char_dimension):
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
			text_region = get_text_region(repr(i+1), dataset_no_dimension, rotate=0)	#no rotate
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
			text_region = get_text_region(tf_name, function_name_dimension)	#rotate
			box = (x_offset_left, y_offset0, x_offset_right, y_offset1)
			im.paste(text_region, box)
			
			#fill in a cell for each dataset_no
			for dataset_no in recurrence_array:
				draw.rectangle((x_offset_left, y_offset1+(dataset_no-1)*dataset_no_dimension[1], \
					x_offset_right, y_offset1+dataset_no*dataset_no_dimension[1]), fill=(0,255,0))
			# write down the recurrence
			text_region = get_text_region(repr(recurrence), dataset_no_dimension)	#rotate
			box = (x_offset_left, y_offset2, x_offset_right, y_offset3)
			im.paste(text_region, box)
		im = im.rotate(270)
		im.save(output_fname)
		del im
		sys.stderr.write("End drawing tf_map.\n")
	
	def run(self):
		"""
		10-31-05
			--form_schema_tables()
			--db_connect()
			--get_char_dimension()
			
			--get_recurrence_rec_array_bs_no_list()
			--get_mt_no2tf_name()
			--draw_tf_map()
			
			--get_recurrence_go_no_rec_array_cluster_id_ls()
			--get_go_no2name()
			--draw_function_map()
			
			--p_gene_id_set_from_gene_p_table()
			--get_no_of_p_funcs_gene_no_go_no_list()
			--draw_gene_function_map()
		"""
		schema_instance = form_schema_tables(self.inputfname, self.acc_cutoff, self.lm_bit)
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		char_dimension = get_char_dimension()
		
		recurrence_rec_array_bs_no_list, no_of_datasets = self.get_recurrence_rec_array_bs_no_list(curs, schema_instance)
		mt_no2tf_name = get_mt_no2tf_name()
		tf_map_output_fname = '%s.tf_map.png'%self.output_prefix
		self.draw_tf_map(recurrence_rec_array_bs_no_list, no_of_datasets, mt_no2tf_name, \
			tf_map_output_fname, self.function_name_length, char_dimension)
		
		recurrence_go_no_rec_array_cluster_id_ls, no_of_datasets = self.get_recurrence_go_no_rec_array_cluster_id_ls(curs, \
			schema_instance.good_cluster_table)
		go_no2name = get_go_no2name(curs)
		no_of_functions = len(recurrence_go_no_rec_array_cluster_id_ls)
		function_map_output_fname = '%s.function_map.png'%self.output_prefix
		go_no2index, function_name_region = self.draw_function_map(recurrence_go_no_rec_array_cluster_id_ls, no_of_datasets,\
			go_no2name, function_map_output_fname, self.function_name_length, char_dimension, no_of_functions)
		
		given_p_gene_set = p_gene_id_set_from_gene_p_table(curs, schema_instance.gene_p_table)
		no_of_p_funcs_gene_no_go_no_list = self.get_no_of_p_funcs_gene_no_go_no_list(curs, \
			schema_instance.p_gene_table, given_p_gene_set)
		gene_function_map_output_fname = '%s.gene_function_map.png'%self.output_prefix
		self.draw_gene_function_map(no_of_p_funcs_gene_no_go_no_list, go_no2index, function_name_region,\
		gene_function_map_output_fname, self.function_name_length, char_dimension, no_of_functions)
		

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema="]
	opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:l:a:p:br", long_options_list)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	inputfname = None
	lm_bit = '00001'
	acc_cutoff = 0.6
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
		elif opt in ("-i"):
			inputfname = arg
		elif opt in ("-l"):
			lm_bit = arg
		elif opt in ("-a"):
			acc_cutoff = float(arg)
		elif opt in ("-p"):
			output_prefix = arg
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r"):
			report = 1
		
	if inputfname and output_prefix:
		instance = DrawMaps(hostname, dbname, schema, inputfname, lm_bit, acc_cutoff, output_prefix, \
			debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
