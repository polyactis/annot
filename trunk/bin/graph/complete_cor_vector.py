#!/usr/bin/env python
"""
Usage: complete_cor_vector.py  -i INPUT_FILE -o OUTPUT_FILE [OPTION] DATADIR

Option:
	DATADIR is the directory of the datasets to be choosen
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	the splat_result table (edge_set)
	-i ..., --input_file=...	gspan format(by graph_merge.py)
	-o ..., --output_file=...	the edge correlation vector file, for CODENSE and COPATH
	-b ..., --label=...	use gene 1(index, default) or 2(no) or 3(id) to label
	-f ..., --dir_files=...	the directory to store the temporary files, /tmp/yh(default)
	-r, --report	report the progress(a number)
	-h, --help              show this help
	
Examples:
	#create a file to dump into database by haiyan_cor_vector2db.py
	complete_cor_vector.py -i sc_merge_gspan -o /tmp/cor_vector ~/datasets/sc_54
	#output cor_vector from database
	complete_cor_vector.py -k sc_54 -i sc_merge_gspan -o /tmp/cor_vector

Description:
	given a summary graph, get all the edges, find the correlations in all
	the datasets, whether significant or insignificant, for each edge,
	output in haiyan's format with gene_index based on the order of nodes
	appearing in the 'v' part of the gspan file.

"""

import sys, os, getopt, csv, re, math, psycopg
import graph_modeling

	
def string_convert_to_three_int(str):
	'''
	a function to convert a string to a float, then get the first three digits after the point
	'''	
	return int(math.modf(float(str)*1000)[1])

class complete_cor_vector:
	'''
	run
		--cor_vector_from_db
		or
		--cor_vector_from_files
	cor_vector_from_files
		--dstruc_loadin
		--edge_tuple_list_output
		--files_sort
		(loop begins)
			--gene_id2expr_array_setup
			--expr_array2expr_list
			--python_call
			--file_combine
	'''
	def __init__(self, hostname, dbname, schema, table, input_file, output_file, label, dir_files, report, dir):
		"""
		run	--dstruc_loadin
			--edge_tuple_list_output
			--files_sort
			(loop)
				gene_index2expr_array_setup
				cor_calculate
				file_combine(self.output_file, self.tmp_file1, self.tmp_file2)
		"""
		self.table = table
		self.input_file = input_file
		self.output_file = output_file
		self.label = int(label)
		self.dir_files = dir_files
		self.dir = dir
		if schema:
			self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
			self.curs = self.conn.cursor()
			self.schema = schema
			self.curs.execute("set search_path to %s"%schema)
		else:
			if self.dir == None:
				sys.stderr.write("Either schema or data directory should be specified\n")
				sys.exit(2)
			else:
				self.files = os.listdir(self.dir)
		
			#check if the temporary directory already exists
			if not os.path.isdir(self.dir_files):
				os.makedirs(self.dir_files)
			else:
				sys.stderr.write("Warning, directory %s already exists.\n"%(self.dir_files))

		self.tmp_file1 = os.path.join(self.dir_files, 'tmp_complete_cor_vector1')
		self.tmp_file2 = os.path.join(self.dir_files, 'tmp_complete_cor_vector2')
		#mapping between gene_no and gene_id
		self.gene_index2no = {}
		self.gene_no2index = {}
		self.gene_index2id = {}
		self.gene_id2index = {}
		#the label dictionary
		self.label_dict = {}
		#edge_tuple_list, constructed based on gene_index
		self.edge_tuple_list = []
		#to setup the 'NA' vector if gene_id not present in some dataset
		self.no_of_cols = 0
		#this stores all the correlations for these edges across different datasets
		self.edge_cor_across_datasets = []
		#to get the number out of the dataset name
		self.p_no = re.compile(r'\d+$')

	def dstruc_loadin(self):
		'''
		This method loads in the data structures from database.
		'''
		sys.stderr.write("Loading Data STructure...")
		#setup three mapping stuff and setup the self.edge_tuple_list
		reader = csv.reader(file(self.input_file), delimiter=' ')
		self.no_of_genes = 0
		for row in reader:
			if row[0] == 'v':
				gene_no = int(row[1])
				gene_id = row[2]
				self.gene_index2no[self.no_of_genes] = gene_no
				self.gene_no2index[gene_no] = self.no_of_genes
				self.gene_index2id[self.no_of_genes] = gene_id
				self.gene_id2index[gene_id] = self.no_of_genes
				self.no_of_genes += 1
			elif row[0] == 'e':
				gene_no1 = int (row[1])
				gene_no2 = int(row[2])
				gene_index1 = self.gene_no2index[gene_no1]
				gene_index2 = self.gene_no2index[gene_no2]
				self.edge_tuple_list.append(gene_index1)
				self.edge_tuple_list.append(gene_index2)
		del reader
		
		#setup up the label dictionary
		gene_index_list = range(self.no_of_genes)
		self.label_dict = {1: gene_index_list,
			2: self.gene_index2no,
			3: self.gene_index2id}
		sys.stderr.write("Done\n")

	def gene_index2expr_array_setup(self, f_path):
		gene_index2expr_array = {}
		reader = csv.reader(file(f_path), delimiter='\t')
		for row in reader:
			gene_id = row[0]
			if gene_id in self.gene_id2index:
				gene_index = self.gene_id2index[gene_id]
				#shorten the float part
				for i in range(1,len(row)):
					if row[i] == 'NA':
						#100000000 is regarded as NAN
						row[i] = float(100000000)
					elif row[i]!='':
						row[i] = float(row[i])
				if row[-1] == '':
					new_row = row[1:-1]
				else:
					new_row = row[1:]
				gene_index2expr_array[gene_index] = new_row
				if self.no_of_cols == 0:
					#hasn't been setup yet
					self.no_of_cols = len(gene_index2expr_array[gene_index])
		del reader
		#makeup those unpresent genese
		for i in range(self.no_of_genes):
			if i not in gene_index2expr_array:
				gene_index2expr_array[i] = [100000000]*self.no_of_cols
		
		return gene_index2expr_array

	def expr_array2expr_list(self, gene_index2expr_array):
		expr_list = []
		for i in range(self.no_of_genes):
			if i in gene_index2expr_array:
				expr_list += gene_index2expr_array[i]
			else:
				expr_list += ['NA']*self.no_of_cols
		return expr_list

	def file_combine(self, output_file, tmp_file1, tmp_file2):
		#use the unix command 'paste' to merge lines of two files
		job_fname = 'paste %s %s > %s'%(output_file, tmp_file1, tmp_file2)
		wl = ['sh', '-c', job_fname]
		os.spawnvp(os.P_WAIT, 'sh', wl)
		#overlapping the first file by 'mv'
		job_fname = 'mv %s %s'%(tmp_file2, output_file)
		wl = ['sh', '-c', job_fname]
		os.spawnvp(os.P_WAIT, 'sh', wl)

	def edge_tuple_list_output(self, output_file):
		of = open(output_file, 'w')
		#start from 0, one step by 2
		for i in range(0, len(self.edge_tuple_list), 2):
			gene_index1 = self.edge_tuple_list[i]
			gene_index2 = self.edge_tuple_list[i+1]
			#haiyan's format requires the first index should be larger than the second index
			if gene_index1>gene_index2:
				of.write("%s\t%s\n"%(self.label_dict[self.label][gene_index1], self.label_dict[self.label][gene_index2]))
			else:
				of.write("%s\t%s\n"%(self.label_dict[self.label][gene_index2], self.label_dict[self.label][gene_index1]))
		of.close()
	
	def files_sort(self, files_list):
		new_files_list = ['']*len(files_list)
		for f in files_list:
			no = int(self.p_no.search(f).group())	#integer conversion
			new_files_list[no-1] = f
		return new_files_list
	
	def cor_calculate(self, output_file, gene_index2expr_array):
		of = open(output_file, 'w')
		#start from 0, one step by 2
		for i in range(0, len(self.edge_tuple_list), 2):
			gene_index1 = self.edge_tuple_list[i]
			gene_index2 = self.edge_tuple_list[i+1]
			y = graph_modeling.ind_min_cor(gene_index2expr_array[gene_index1], gene_index2expr_array[gene_index2])
			#get the integer part of the float*1000
			of.write("%s\n"%(int(math.modf(y.value*1000)[1]) ) )
		of.close()
	
	def cor_vector_from_db(self, output_file):
		writer = csv.writer(open(output_file, 'w'), delimiter='\t')
		#start from 0, one step by 2
		for i in range(0, len(self.edge_tuple_list), 2):
			gene_index1 = self.edge_tuple_list[i]
			gene_index2 = self.edge_tuple_list[i+1]
			gene_no1 = self.gene_index2no[gene_index1]
			gene_no2 = self.gene_index2no[gene_index2]
			if gene_no1 < gene_no2:
				edge_string = '{%d,%d}'%(gene_no1, gene_no2)
			else:
				edge_string = '{%d,%d}'%(gene_no2, gene_no1)
			self.curs.execute("select cor_vector from %s where edge_name='%s'"%(self.table,edge_string))
			rows = self.curs.fetchall()
			if len(rows) == 0:
				sys.stderr.write('%s not found in edge_cor_vector\n'%edge_string)
				sys.exit(1)
			cor_vector = rows[0][0][1:-1].split(',')
			cor_vector = map(string_convert_to_three_int, cor_vector)
			#haiyan's format requires the first index should be larger than the second index
			if gene_index1> gene_index2:
				row = [gene_index1, gene_index2]
			else:
				row = [gene_index2, gene_index1]
			row += cor_vector
			writer.writerow(row)
		del writer

	
	def run(self):
		self.dstruc_loadin()
		if self.schema:
			self.cor_vector_from_db(self.output_file)
		else:
			self.cor_vector_from_files()

	def cor_vector_from_files(self):
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(self.files))
		self.edge_tuple_list_output(self.output_file)
		#sort all the files based on the dataset number, to order the columns of the outputed edge correlation vector
		self.files = self.files_sort(self.files)
		for f in self.files:
			sys.stderr.write("%d/%d:\t%s\n"%(self.files.index(f)+1,len(self.files),f))
			f_path = os.path.join(self.dir, f)
			#get the expression values only for those genes
			gene_index2expr_array = self.gene_index2expr_array_setup(f_path)
			self.cor_calculate(self.tmp_file1, gene_index2expr_array)
			"""#convert the expression values dictionary to a list for python_call
			expr_list = self.expr_array2expr_list(gene_index2expr_array)
			#call the C++ module
			python_call(self.tmp_file1, self.edge_tuple_list, expr_list, self.no_of_genes)
			"""
			#combine the self.output_file and self.tmp_file1 into self.output_file by means of self.tmp_file2
			self.file_combine(self.output_file, self.tmp_file1, self.tmp_file2)
		os.remove(self.tmp_file1)
		#self.tmp_file2 is removed in the last call to file_combine
		os.rmdir(self.dir_files)

	
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:i:o:b:f:r", ["help",  "hostname=", "dbname=", "schema=", "table=", \
			"input_file=", "output_file=", "label=", "dir_files=", "report"])
	except:
		print __doc__
		sys.exit(2)

	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = 'edge_cor_vector'
	input_file = None
	output_file = None
	label = 1
	dir_files = '/tmp/yh'
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
		elif opt in ("-t", "--table"):
			table = arg
		elif opt in ("-i", "--input_file"):
			input_file = arg
		elif opt in ("-o", "--output_file"):
			output_file = arg
		elif opt in ("-b", "--label"):
			label = int(arg)
		elif opt in ("-f", "--dir_files"):
			dir_files = arg
		elif opt in ("-r", "--report"):
			report = 1
			
	if len(args) == 1:
		data_dir = args[0]
	else:
		data_dir = None
	if input_file and output_file:
		instance = complete_cor_vector(hostname, dbname, schema, table, \
			input_file, output_file, label, dir_files, report, data_dir)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
