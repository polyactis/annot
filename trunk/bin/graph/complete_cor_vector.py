#!/usr/bin/env mpipython
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
	-s ..., --significance_file=...	the edge significance flag vector file.
	-b ..., --label=...	use gene 1(index, default) or 2(no) or 3(id) to label
	-g ..., --gph_dir=...	the directory where the graph files are stored
		if p_value_cut_off and cor_cut_off are both 0, get the corCut from gph_dir
	-p ..., --p_value_cut_off=...	the p_value_cut_off for an edge to be significant, 0.01(default)
	-c ..., --cor_cut_off=...	the cor_cut_off for an edge to be significant, 0.6(default)
		priority gives p_value_cut_off, cor_cut_off is used only when p_value_cut_off=0
	-r, --report	report the progress(a number)
	-u, --debug	enable debugging output.
	-h, --help              show this help
	
Examples:
	#create a file to dump into database by haiyan_cor_vector2db.py
	mpirun.lam N complete_cor_vector.py -i sc_merge_gspan -o /tmp/cor_vector
		-s /tmp/edge_significance_vector ~/datasets/sc_54
	
	#output cor_vector from database
	complete_cor_vector.py -k sc_54 -i sc_merge_gspan -o /tmp/cor_vector
	
	#read the corCut from top 1% gph files.
	complete_cor_vector.py -i sc_merge_gspan -o /tmp/cor_vector -g gph_result/sc
		-p 0 -c 0 -s /tmp/edge_significance_vector ~/datasets/sc_54

Description:
	given a summary graph, get all the edges, find the correlations in all
	the datasets, and also its significance flag, for each edge,
	output in haiyan's format with gene_index based on the order of nodes
	appearing in the 'v' part of the gspan file. 03-02-05, significance flag
	vector is output in another file.

"""

import sys, os, getopt, csv, re, math, psycopg
import graph_modeling
from Scientific import MPI
from sets import Set

	
def string_convert_to_three_int(str):
	'''
	a function to convert a string to a float, then get the first three digits after the point
	'''	
	return int(math.modf(float(str)*1000)[1])

class complete_cor_vector:
	'''
	
	03-02-05
		Now, graph_modeling can return the significance flag for an edge.
		Add a significance_file to hold those flags.
		
	05-14-05
		modified to be MPI-capable
		
		cor_vector_from_files() almost overhauled.
		add collect_and_merge_output(), node_fire() and modify file_combine()
	
	'''
	def __init__(self, data_dir, hostname, dbname, schema, table, input_file, output_file, \
		significance_file, label, gph_dir=None, p_value_cut_off=0.01, cor_cut_off=0.6, \
		report=0, debug=0):
		"""
		05-14-05
			add debug flag.
			dir_files is useless.
		06-30-05
			dir_files renamed to be gph_dir
		"""
		self.dir = data_dir
		self.schema = schema
		self.table = table
		self.input_file = input_file
		self.output_file = output_file
		self.significance_file = significance_file
		self.label = int(label)
		self.gph_dir = gph_dir
		self.p_value_cut_off = float(p_value_cut_off)
		self.cor_cut_off = float(cor_cut_off)
		self.debug = int(debug)
		
		if self.schema:
			self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
			self.curs = self.conn.cursor()
			self.curs.execute("set search_path to %s"%self.schema)
		else:
			if self.dir == None:
				sys.stderr.write("Either schema or data directory should be specified\n")
				sys.exit(2)
		
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
		
		#04-03-05, Important bug.
		self.gene_cut_off = 8

	def dstruc_loadin(self, communicator):
		'''
		This method loads in the data structures from database.
		05-14-05
			add some stderr's
		'''
		sys.stderr.write("Node %s Loading Data STructure...\n"%communicator.rank)
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
		sys.stderr.write("Node %s: Structure Loading Done\n"%communicator.rank)

	def gene_index2expr_array_setup(self, f_path):
		gene_index2expr_array = {}
		reader = csv.reader(file(f_path), delimiter='\t')
		for row in reader:
			gene_id = row[0]
			if gene_id in self.gene_id2index:
				gene_index = self.gene_id2index[gene_id]
				no_of_nas = 0
				#shorten the float part
				for i in range(1,len(row)):
					if row[i] == 'NA':
						#100000000 is regarded as NAN
						row[i] = float(100000000)
						no_of_nas += 1
					elif row[i]!='':
						row[i] = float(row[i])
				if row[-1] == '':
					new_row = row[1:-1]
				else:
					new_row = row[1:]
				if self.no_of_cols == 0:
					#hasn't been setup yet
					self.no_of_cols = len(new_row)
				if len(new_row)-no_of_nas >= self.gene_cut_off:
					#04-03-05 important bug. another difference between graph.cc and complete_cor_vector+graph_modeling.
					gene_index2expr_array[gene_index] = new_row

		del reader
		#makeup those unpresent genes
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
		"""
		05-14-05
			concatenate output_file and tmp_file1 to be output_file via tmp_file2
			
			replace os.spawnvp with os.system
		"""
		#use the unix command 'paste' to merge lines of two files
		job_fname = 'paste %s %s > %s'%(output_file, tmp_file1, tmp_file2)
		os.system(job_fname)	#spawnvp gets dead under MPI
		#wl = ['sh', '-c', job_fname]
		#os.spawnvp(os.P_WAIT, 'sh', wl)
		#overlapping the first file by 'mv'
		job_fname = 'mv %s %s'%(tmp_file2, output_file)
		os.system(job_fname)
		#wl = ['sh', '-c', job_fname]
		#os.spawnvp(os.P_WAIT, 'sh', wl)

	def edge_tuple_list_output(self, output_file):
		"""
		05-14-05
			add some stderr's
		"""
		sys.stderr.write("Outputting edge_tuples to %s..."%output_file)
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
		sys.stderr.write("Done.\n")
	
	def files_sort(self, files_list):
		new_files_list = ['']*len(files_list)
		for f in files_list:
			no = int(self.p_no.search(f).group())	#integer conversion
			new_files_list[no-1] = f
		return new_files_list
	
	def cor_calculate(self, output_file, significance_file, gene_index2expr_array, corCut_list, file_index):
		"""
		06-30-05
			add corCut_list and file_index, used to overwrite the internal cutoff of graph_modeling
		"""
		of = open(output_file, 'w')
		sf = open(significance_file, 'w')
		#start from 0, one step by 2
		for i in range(0, len(self.edge_tuple_list), 2):
			gene_index1 = self.edge_tuple_list[i]
			gene_index2 = self.edge_tuple_list[i+1]
			edge_data = graph_modeling.ind_min_cor(gene_index2expr_array[gene_index1], gene_index2expr_array[gene_index2])
			#get the integer part of the float*1000
			of.write("%s\n"%(int(math.modf(edge_data.value*1000)[1]) ) )
			if corCut_list:	#06-30-05	if corCut_list, overwrite the significance flag
				sf.write("%s\n"%int(edge_data.value>=corCut_list[file_index]))
			else:
				sf.write("%s\n"%(edge_data.significance))
		sf.close()
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
	
	def collect_and_merge_output(self, final_fname, no_of_datasets):
		"""
		05-14-05
			concatenate all files(=no_of_datasets) into final_fname
			use file_combine() recursively
		05-16-05
			don't use file_combine(), the 'paste' and 'mv' in file_combine() become slow as final_fname gets larger and larger 
			paste all files at once and delete them one by one.
			
		"""
		sys.stderr.write("Collecting for %s...\n"%final_fname)
		job_word_list = ['paste']	#later concatenate this list to be a command
		dataset_fname_list = []	#the list of files to be removed
		for i in range(no_of_datasets+1):	#0th is edge tuple file, 1-no_of_datasets are data files
			dataset_fname = "%s_%s"%(final_fname, i)
			job_word_list.append(dataset_fname)
			dataset_fname_list.append(dataset_fname)
		job_word_list.append('>')	#this is different from '>>', the final_fname is truncated first.
		job_word_list.append(final_fname)
		exit_code = os.system(' '.join(job_word_list))	#recursive until success
		while exit_code:
			exit_code = os.system(' '.join(job_word_list))
			
		#remove the dataset_fname
		for dataset_fname in dataset_fname_list:
			os.remove(dataset_fname)
		sys.stderr.write("Done.\n")
	
	def node_fire(self, dir, files, file_index, cor_fname, sig_fname, corCut_list):
		"""
		05-14-05
			mpi unit
			
			--gene_index2expr_array_setup()
			--cor_calculate()
		05-16-05
			suffix each unit file with index+1
		
		06-30-05
			add corCut_list and file_index, passed to cor_calculate()
		"""
		f = files[file_index]
		f_path = os.path.join(dir, f)	#the path to the dataset file
		gene_index2expr_array = self.gene_index2expr_array_setup(f_path)
		cor_tmp_file = "%s_%s"%(cor_fname, file_index+1)	#index+1
		sig_tmp_file = "%s_%s"%(sig_fname, file_index+1)
		self.cor_calculate(cor_tmp_file, sig_tmp_file, gene_index2expr_array, corCut_list, file_index)
		
	def mpi_synchronize(self, communicator):
		"""
		05-14-05
			copied from MpiBiclustering.py
		"""
		sys.stdout.flush()
		sys.stderr.flush()
		communicator.barrier()
	
	def get_corCut_list(self, gph_dir):
		"""
		06-30-05
			get the top 1% corCut from the files in gph_dir
			
			--get_corCut_from_gph_file
		"""
		files = os.listdir(gph_dir)
		corCut_list = [0]*len(files)
		for file in files:
			no = int(self.p_no.search(file).group())
			f_path = os.path.join(gph_dir, file)	#the path to the dataset file
			corCut = self.get_corCut_from_gph_file(f_path)
			corCut_list[no-1] = corCut
		return corCut_list
	
	def get_corCut_from_gph_file(self, gph_file):
		"""
		06-30-05
			the first edge in gph_file is the lowest correlation due to the nature of priority_queue
		"""
		reader = csv.reader(open(gph_file, 'r'), delimiter='\t')
		reader.next()	#skip the first line containing gph name
		row = reader.next()	#one line is 'e       Mm.153061       Mm.3295 0.361159'
		del reader
		return float(row[3])
	
	def cor_vector_from_files(self, communicator, dir, gph_dir, cor_fname, sig_fname, p_value_cut_off, cor_cut_off):
		"""
		05-14-05
			modify to be mpi form, feed mode(in MpiBiclustering.py and MpiGraphModeling.py)
		
		05-16-05
			output the edge tuple into 0th file.
			
		06-30-05	if both 0, get the corCut from the top 1% graph files	
		
			--files_sort
			if node_rank==0:
				--edge_tuple_list_output()
				--edge_tuple_list_output()
			else:
				--graph_modeling.cor_cut_off_vector_construct()
				--get_corCut_list()
				
			if node_rank==0:
				(send signal to other nodes)
			else:
				--node_fire()
					--gene_index2expr_array_setup()
					--cor_calculate()
			if node_rank==0:
				--collect_and_merge_output()
		"""
		files = os.listdir(dir)
		#sort all the files based on the dataset number, to order the columns of the outputed edge correlation vector
		files = self.files_sort(files)
		
		file_index_list = range(len(files))
		node_rank = communicator.rank
		
		if p_value_cut_off ==0 and cor_cut_off == 0 and gph_dir==None:
			sys.stderr.write("p_value_cut_off and cor_cut_off both are 0, but no gph_dir. Aborted.\n")
			sys.exit(3)
		if node_rank == 0:
			#output the name first
			self.edge_tuple_list_output("%s_0"%cor_fname)	#05-16-05 output the edge tuple into 0th file.
			self.edge_tuple_list_output("%s_0"%sig_fname)
		else:
			#set the cor_cut_off_vector, internal structure of graph_modeling
			graph_modeling.cor_cut_off_vector_construct(p_value_cut_off, cor_cut_off)
			if p_value_cut_off ==0 and cor_cut_off == 0:	#06-30-05	if both 0, get the corCut from the top 1% graph files
				corCut_list = self.get_corCut_list(gph_dir)
			else:
				corCut_list = []
		
		self.mpi_synchronize(communicator)
		
		if node_rank == 0:
			sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))

			seed_utilized = Set()
			for node in range(1, communicator.size):
				if len(file_index_list)==0:	#if #nodes > #jobs, tell those nodes to break their listening loop.
					stop_signal = "-1"
					communicator.send(stop_signal, node, 0)	#no more jobs, stop that node,
					if self.debug:
						sys.stderr.write("node %s stopped.\n"%node)
				else:
					input_file_index = file_index_list.pop(0)	#the first item poped first.
					communicator.send(repr(input_file_index), node, 0)	#string format
					if self.debug:
						sys.stderr.write("Node %s schedule a job to %s\n"%(node_rank, node))
					seed_utilized.add(node)
			
			received_value, source, tag = communicator.receiveString(None, None)	#listen
			while received_value:		#??check what the received_value is
				if len(file_index_list) == 0:	#first check if there're still files left, otherwise pop(0) raises error.
					stop_signal = "-1"
					communicator.send(stop_signal, source, 0)	#no more jobs, stop that node,
					if self.debug:
						sys.stderr.write("node %s stopped.\n"%source)
					seed_utilized.remove(source)
					if len(seed_utilized) == 0:	#all seed used have finished their jobs
						break
				else:
					input_file_index = file_index_list.pop(0)
					if input_file_index:
						communicator.send(repr(input_file_index), source, 0)	#string format,
						if self.debug:
							sys.stderr.write("Node %s get one more job\n"%source)
				received_value, source, tag = communicator.receiveString(None, None)	#listen
		else:
			received_data, source, tag = communicator.receiveString(0, None)	#get data from node 0,
				#04-24-05 the array is one-dimension no matter what dimension the original array is
			while received_data:
				if received_data=="-1":	#stop signal
					if self.debug:
						sys.stderr.write("node %s breaked.\n"%node_rank)
					break
				else:
					input_file_index = int(received_data)	#convert it to integer
					sys.stderr.write("node %s working on %s...\n"%(node_rank, received_data))
					self.node_fire(dir, files, input_file_index, cor_fname, sig_fname, corCut_list)
					sys.stderr.write("node %s work on %s finished.\n"%(node_rank, received_data))
					communicator.send("finished", 0, node_rank)
					
				received_data, source, tag = communicator.receiveString(0, None)	#get data from node 0
		
		self.mpi_synchronize(communicator)
		
		if node_rank==0:
			self.collect_and_merge_output(cor_fname, len(files))
		elif node_rank==1:
			self.collect_and_merge_output(sig_fname, len(files))
			
	def run(self):
		"""
		05-14-05
		06-30-05
			add gph_dir
		
		--dstruc_loadin()
		if self.schema:
			--cor_vector_from_db
			or
			--cor_vector_from_files

		"""
		communicator = MPI.world.duplicate()
		
		self.dstruc_loadin(communicator)
		
		self.mpi_synchronize(communicator)
		
		if self.schema:
			self.cor_vector_from_db(self.output_file)
		else:
			if self.significance_file==None:
				sys.stderr.write("Error: where's the significance file?\n")
				sys.exit(2)
			self.cor_vector_from_files(communicator, self.dir, self.gph_dir, self.output_file, \
				self.significance_file, self.p_value_cut_off, self.cor_cut_off)

	
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:i:o:s:b:g:p:c:ru", ["help",  "hostname=", \
			"dbname=", "schema=", "table=", "input_file=", "output_file=", "significance_file=", \
			"label=", "gph_dir=", "p_value_cut_off=", "cor_cut_off=", "report", "debug"])
	except:
		print __doc__
		sys.exit(2)

	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = 'edge_cor_vector'
	input_file = None
	output_file = None
	significance_file = None
	label = 1
	gph_dir = None
	p_value_cut_off = 0.01
	cor_cut_off = 0.6
	report = 0
	debug = 0
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
		elif opt in ("-s", "--significance_file"):
			significance_file = arg
		elif opt in ("-b", "--label"):
			label = int(arg)
		elif opt in ("-g", "--gph_dir"):
			gph_dir = arg
		elif opt in ("-p", "--p_value_cut_off"):
			p_value_cut_off = float(arg)
		elif opt in ("-c", "--cor_cut_off"):
			cor_cut_off = float(arg)
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-u", "--debug"):
			debug = 1
			
	if len(args) == 1:
		data_dir = args[0]
	else:
		data_dir = None
	if input_file and output_file:
		instance = complete_cor_vector(data_dir, hostname, dbname, schema, table, \
			input_file, output_file, significance_file, label, gph_dir, p_value_cut_off, \
			cor_cut_off, report, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
