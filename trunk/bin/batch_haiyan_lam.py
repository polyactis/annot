#!/usr/bin/env python
"""
Usage: batch_haiyan_lam.py --mp=FILEPREFIX [OPTIONS]

Option:
	-m ..., --run_mode=...	0(default), 1, 2
	-n ..., --genenum=...	number of genes, 6661(default)
	-p ..., --svnum=...	number of edges, 805939(default)
	-l ..., --sv_length=...	length of edge correlation vector, #datasets, 54(default)
	-t ..., --ttablefile=...	the ttable file, ttableFromMatlabt1p-3.txt in the id(default)
	
	#following parameters can be given in a list form
	-r ..., --cut_loop_num_list=...	cut loop num, only for copath, 2(default)
	-g ..., --min_graph_size_list=...	the minimum cluster(subgraph) size, 5(default)
	-e ..., --min_edge_freq_list=...	6(default)
	-d ..., --first_density_cutoff_list=...	0.4(default)
	-q ..., --second_density_cutoff_list=...	0.4(default)
	-s ..., --max_pre_graph_size_list=...	the maximum subgraph size to do 2nd-order clustering, 80(default)
	-c ..., --conn_perc_list=...	connect perc restoring the condensedcluster, 0.5(default)
	
	-u ..., --max_degree=...	maximum degree of freedom for user's T table, 500(default)
	
	-y ..., --type=...	0(coden, default), 1(copath)
	--id=...	input directory, '~/bin/hhu_clustering/data/input' (default)
	--od=...	output directory, '~/bin/hhu_clustering/data/ouput' (default)
	-bd=...	the binary directory, '~/bin/hhu_clustering/bin' (default)
	--mp=...	the prefix of the input matrix filename, i.e. 'sc_54_6661_merge_6'
		Their suffix should be .matrix. The corresponding correlation vector files are suffixed as .cor_vector.
	--op=...	output file prefix, default is 
		TypeMpGMin_graph_sizeEMin_edge_freqDFirst_density_cutoffQSecond_density_cutoffSMax_pre_graph_sizeCConn_perc
		Note: those floats are taken the tenth and hundredth digit. i.e. 0.4 -> 40.
	--rr=...	rank range, i.e. 0-3,5,7-9 means useing 0,1,2,3,5,7,8,9
		This is used to specify the nodes to use.
	--js=...	job file starting no, 0 (default)
	-h, --help              show this help
	
Examples:
	#simplest:
	batch_haiyan_lam.py --mp='sc_54_6661_merge_6' --rr=0-3

	#normal:
	batch_haiyan_lam.py --mp='sc_54_6661_merge_6' -n 6661 -p 805939 -l 54
		-e 6 -d 0.4,0.6,0.8 -q 0.6 -y 1 --rr=0-3,5,7-9
	
	
Description:
	This program is used to start several codense or copath, with different settings
	on several nodes via lam/mpi.


"""


#from Scientific import MPI
import sys, os, math, getopt, time, csv

class batch_haiyan_lam:
	"""
	03-16-05
		parameters passed to coden or copath should be in string form.
	03-19-05
		It's not a mpipython program anymore. It uses at to schedule job files which
		are like ~/qjob/batch_haiyan_lam#.sh and these files call lam_sub.py(a mpipython program)
		to submit jobs separately. The drawback of mpipython program is that it exits
		once one job stops abnormally.
	"""
	def __init__(self, run_mode='0', genenum='6661', svnum='805939', sv_length='54',\
		ttablefile='ttableFromMatlabt1p-3.txt', cut_loop_num_list=['2'], min_graph_size_list=['5'], \
		min_edge_freq_list=['6'], first_density_cutoff_list=['0.4'], second_density_cutoff_list=['0.4'],\
		max_pre_graph_size_list=['80'], conn_perc_list=['0.5'], max_degree='500', type=0, \
		id=None, od=None, bd=None, mp='sc_54_6661_merge_6', op=None, rank_range=None, js=0):
		"""
		03-16-05
			those None's in the argument list will be given default values
		"""
		self.run_mode = run_mode
		self.genenum = genenum
		self.svnum = svnum
		self.sv_length = sv_length
		self.ttablefile = ttablefile
		self.cut_loop_num_list = cut_loop_num_list
		self.min_graph_size_list = min_graph_size_list
		self.min_edge_freq_list = min_edge_freq_list
		self.first_density_cutoff_list = first_density_cutoff_list
		self.second_density_cutoff_list = second_density_cutoff_list
		self.max_pre_graph_size_list = max_pre_graph_size_list
		self.conn_perc_list = conn_perc_list
		self.max_degree = max_degree
		self.type = int(type)
		self.id = id
		self.od = od
		self.bd = bd
		self.mp = mp
		self.op = op
		self.rank_range = rank_range
		self.js = int(js)

	def parameter_list_init(self, comm):
		"""
		03-16-05
			input: comm
			output: parameter_list
			
			some other class wide parameters, self.input_matrix_file....
			
			initialization for other arguments and parameter list to call codense or copath
		03-17-05
		
		03-19-05
			1. no default for rank_range, must be specified by user
			2. program_path is parameter_list[0]
			
		"""
		#None arguments
		if self.id == None:
			self.id = os.path.join(os.path.expanduser('~'), 'bin/hhu_clustering/data/input')
		if self.od == None:
			self.od = os.path.join(os.path.expanduser('~'), 'bin/hhu_clustering/data/output')
		if self.bd == None:
			self.bd = os.path.join(os.path.expanduser('~'), 'bin/hhu_clustering/bin')
		
		#additional class wide parameters
		self.input_matrix_file = os.path.join(self.id, self.mp+'.matrix')
		self.input_sv_file = os.path.join(self.id, self.mp+'.cor_vector')
		self.input_ttable_file = os.path.join(self.id, self.ttablefile)
		if self.type==0:
			self.binary = 'coden'
		elif self.type==1:
			self.binary = 'copath'
		program_path = os.path.join(self.bd, self.binary)
		
		parameter_list = []
		for cut_loop_num in self.cut_loop_num_list:
			for min_graph_size in self.min_graph_size_list:
				for min_edge_freq in self.min_edge_freq_list:
					for first_density_cutoff in self.first_density_cutoff_list:
						for second_density_cutoff in self.second_density_cutoff_list:
							for max_pre_graph_size in self.max_pre_graph_size_list:
								for conn_perc in self.conn_perc_list:
									if self.op == None:
										#initialize the output prefix
										#get the tenth and hundredth digit
										int_first_density_cutoff = int(float(first_density_cutoff)*100)
										int_second_density_cutoff = int(float(second_density_cutoff)*100)
										int_conn_perc = int(float(conn_perc)*100)
										op = self.binary+self.mp+'G%sE%sD%sQ%sS%sC%s'%\
		(min_graph_size, min_edge_freq, int_first_density_cutoff, int_second_density_cutoff, max_pre_graph_size, int_conn_perc)
									if self.type==0:
										parameter_list.append([program_path,\
		'-m', self.run_mode, '-i', self.input_matrix_file, '-n', self.genenum, '-v', self.input_sv_file,\
		'-p', self.svnum, '-l', self.sv_length, '-t', self.input_ttable_file, '-o', op, '-g', min_graph_size,\
		'-e', min_edge_freq, '-d', first_density_cutoff, '-q', second_density_cutoff, \
		'-s', max_pre_graph_size, '-c', conn_perc, '-u', self.max_degree])
									elif self.type==1:
										#cut_loop_num is the difference
										parameter_list.append([program_path, '-r', cut_loop_num,\
		'-m', self.run_mode, '-i', self.input_matrix_file, '-n', self.genenum, '-v', self.input_sv_file,\
		'-p', self.svnum, '-l', self.sv_length, '-t', self.input_ttable_file, '-o', op, '-g', min_graph_size,\
		'-e', min_edge_freq, '-d', first_density_cutoff, '-q', second_density_cutoff, \
		'-s', max_pre_graph_size, '-c', conn_perc, '-u', self.max_degree])

		return parameter_list
		
	def _node_fire(self, comm, program_path, parameter_list, node_rank_list):
		"""
		03-16-05
			input: comm, program_path, parameter_list, node_rank_list
			output: print process information
			
			each node fires some jobs based on its own rank
		03-19-05
			not a mpipython program anymore
		"""


	def run(self):
		"""
		03-16-05
			distribute coden or copath along the nodes
		
		03-19-05
			no communicator anymore
		"""
		comm = None
		"""
		#a communicator
		comm = MPI.world.duplicate()
		if self.rank_range[0] > comm.size or self.rank_range[1] > comm.size:
			sys.stderr.write("Error: invalid rank range: %s"%repr(self.rank_range))
			sys.exit(2)
		"""
		parameter_list = self.parameter_list_init(comm)
		#hour:minute, minute is current_minute + 2, if time is like 11:61, it's fine.
		time_tuple = time.localtime()
		time_to_run_jobs = "%s:%s"%(time_tuple[3], time_tuple[4]+2)
		#the path of lam_sub.py
		lam_sub_path = '~/script/annot/bin/lam_sub.py'
		
		"""
		testing

		program_path = 'ls'
		parameter_list = [['ls'], ['ls', '/'], ['ls', '/usr/local']]
		"""
		
		#map parameters to each node_rank
		node_rank2parameter_setting = {}
		for i in range(len(parameter_list)):
			#remainder is the node_rank
			index = int(math.fmod(i, len(self.rank_range)))
			node_rank = self.rank_range[index]
			if node_rank not in node_rank2parameter_setting:
				node_rank2parameter_setting[node_rank] = [parameter_list[i]]
			else:
				node_rank2parameter_setting[node_rank].append(parameter_list[i])
		job_number = self.js
		for node_rank in node_rank2parameter_setting:
			#some nodes will be idle if there're more nodes than jobs
			for parameter in node_rank2parameter_setting[node_rank]:
				job_fname = os.path.join(os.path.expanduser('~'), 'qjob/batch_haiyan_lam%s.sh'%job_number)
				job_f = open(job_fname, 'w')
				writer = csv.writer(job_f, delimiter=' ')
				#change to the ouput directory first, because copath and coden output to the current directory
				writer.writerow(['cd', self.od, ' '])	#the endline is attached to self.od and cd says an error
				jobrow = ['mpirun', 'n%s'%node_rank, lam_sub_path] + parameter
				writer.writerow(['echo']+jobrow)	#print the commandline
				writer.writerow(jobrow)
				#close the file
				del writer
				job_f.close()
				
				print "node: %s, at %s, parameter: %s"%(node_rank, time_to_run_jobs, repr(parameter))
				#schedule it
				wl = ['at', '-mf', job_fname, time_to_run_jobs]	#-m, even if no output, mail me.
				os.spawnvp(os.P_WAIT, 'at', wl)
				
				job_number+=1
						

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hm:n:p:l:t:r:g:e:d:q:s:c:u:y:", ["help","run_mode=",\
			"genenum=", "svnum=", "sv_length=", "ttablefile=", "cut_loop_num_list=", "min_graph_size_list=",\
			"min_edge_freq_list=", "first_density_cutoff_list=", "second_density_cutoff_list=", \
			"max_pre_graph_size_list=", "conn_perc_list=", "max_degree=", "type=", "id=", "od=", \
			"bd=", "mp=", "op=", "rr=", "js="])
	except:
		print __doc__
		sys.exit(2)
	
	run_mode = '0'
	genenum = '6661'
	svnum = '805939'
	sv_length = '54'
	ttablefile ='ttableFromMatlabt1p-3.txt'
	cut_loop_num_list = ['2']
	min_graph_size_list = ['5']
	min_edge_freq_list = ['6']
	first_density_cutoff_list = ['0.4']
	second_density_cutoff_list = ['0.4']
	max_pre_graph_size_list = ['80']
	conn_perc_list = ['0.5']
	max_degree = '500'
	type = '0'
	id = None
	od = None
	bd = None
	mp = None
	op = None
	rank_range = []
	js = 0
	
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-m", "--run_mode"):
			run_mode = arg
		elif opt in ("-n", "--genenum"):
			genenum = arg
		elif opt in ("-p", "--svnum"):
			svnum = arg
		elif opt in ("-l", "--sv_length"):
			sv_length = arg
		elif opt in ("-t", "--ttablefile"):
			ttablefile = arg
		elif opt in ("-r", "--cut_loop_num_list"):
			cut_loop_num_list = arg.split(',')

		elif opt in ("-g", "--min_graph_size_list"):
			min_graph_size_list = arg.split(',')

		elif opt in ("-e", "--min_edge_freq_list"):
			min_edge_freq_list = arg.split(",")
		elif opt in ("-d", "--first_density_cutoff_list"):
			first_density_cutoff_list = arg.split(",")
		elif opt in ("-q", "--second_density_cutoff_list"):
			second_density_cutoff_list = arg.split(",")
		elif opt in ("-s", "--max_pre_graph_size_list"):
			max_pre_graph_size_list = arg.split(',')
		elif opt in ("-c", "--conn_perc_list"):
			conn_perc_list = arg.split(',')
		elif opt in ("-u", "--max_degree"):
			max_degree = arg
		elif opt in ("-y", "--type"):
			type = int(arg)
		elif opt in ("--id"):
			id = arg
		elif opt in ("--od"):
			od = arg
		elif opt in ("--bd"):
			bd = arg
		elif opt in ("--mp"):
			mp = arg
		elif opt in ("--op"):
			op = arg
		elif opt in ("--rr"):
			#a sample rr is like '3-5,6,8-10'
			rr = arg.split(',')
			for rank in rr:
				rank = rank.split('-')
				rank = map(int, rank)
				if len(rank)==2:
					rank_range += range(rank[0], rank[1]+1)
				else:
					rank_range += rank			
		elif opt in ("--js"):
			js = int(arg)
			
	if mp and len(rank_range)>0:
		instance = batch_haiyan_lam(run_mode, genenum, svnum, sv_length,\
			ttablefile, cut_loop_num_list, min_graph_size_list, \
			min_edge_freq_list, first_density_cutoff_list, second_density_cutoff_list,\
			max_pre_graph_size_list, conn_perc_list, max_degree, type,\
			id, od, bd, mp, op, rank_range, js)
		
		instance.run()
		
	else:
		print __doc__
		sys.exit(2)
