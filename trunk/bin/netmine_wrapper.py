#!/usr/bin/env mpipython.lam
"""
Usage: netmine_wrapper.py --mp=FILEPREFIX [OPTIONS]

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
	
	-h ..., --match_cut=...	4.0(default)
	-j ..., --intersect2union_cut=...	0.4(default)
	-w ..., --euclidean_ratio=...	0.2(default)
	-z ..., --selection_code=...	0110(default), 4-digit string, the leftmost digit(1st digit) means using match measure, 
		2nd digit is intersect/union measure, 3rd is euler/max measure, 4th is pearson correlation
		for example, if you combine 3rd and 4th. you use -z 0011
	
	-u ..., --max_degree=...	maximum degree of freedom for user's T table, 500(default)
	-k ...,	this argument specifies the selected dataset index like, 37/0,1-3,4,5-36
	
	-y ..., --type=...	0(coden, default), 1(copath) (IGNORE)
	--id=...	input directory, '~/bin/hhu_clustering/data/input' (default)
	--od=...	output directory, '~/bin/hhu_clustering/data/ouput/netmine' (default)
	--bd=...	the binary directory, '~/bin/hhu_clustering/bin' (default)
	--mp=...	the prefix of the input matrix filename, i.e. 'sc_54_6661_merge_6'
		Their suffix should be .matrix. The corresponding correlation vector files are suffixed as .cor_vector.
	--op=...	output file prefix, default is 
		TypeMpGMin_graph_sizeEMin_edge_freqDFirst_density_cutoffQSecond_density_cutoffSMax_pre_graph_sizeCConn_perc
		Note: those floats are taken the tenth and hundredth digit. i.e. 0.4 -> 40.
	--rr=...	rank range, i.e. 0-3,5,7-9 means useing 0,1,2,3,5,7,8,9
		This is used to specify the nodes to use. (IGNORE, mpirun specifies this)
	--js=...	job file starting no, 0 (default)
	--debug	enable debug
	--help	show this help
	
Examples:
	#simplest:
	mpirun N netmine_wrapper.py --mp='sc_54_6661_6'

	#normal:
	mpirun n0-17 netmine_wrapper.py --mp='sc54_5' -n 6661 -p 1342902 -l 54
		-e 6 -q 0.6 -z0001
	
	mpirun n0-17 netmine_wrapper.py --mp='sc54_5' -n 6661 -p 1342902 -l54
		-e6 -q0.4 -w0.4 -z0010
	
	
Description:
	This program wraps netmine and netmine2nd(two parts of the copath).
	This node will run netmine. Nodes in the rank_range will run netmine2nd.
"""


import sys, os, math, getopt, time, csv, Numeric
from Scientific import MPI
from codense.common import system_call, mpi_schedule_jobs, mpi_synchronize
from sets import Set
		
def run_netmine2nd(cluster_no, netmine2nd_parameter_list):
	"""
	05-19-05
		similar to node_fire() of class netmine_wrapper
	"""
	import sys
	#get the output filename prefix
	ofname_prefix = netmine2nd_parameter_list[20]
	jobrow = netmine2nd_parameter_list + ['-f', repr(cluster_no)]
	exit_code = system_call("%s"%' '.join(jobrow))
	ofname = '%s_%sh'%(ofname_prefix, cluster_no)	#04-09-05 haiyan puts 'h' after cluster_no.??
	return ofname

class netmine_wrapper:
	"""
	04-01-05
		start
	04-09-05
		use mpi, os.system solves the problem of calling external program in mpi environment.
	"""
	def __init__(self, run_mode='0', genenum='6661', svnum='805939', sv_length='54',\
		ttablefile='ttableFromMatlabt1p-3.txt', cut_loop_num_list=['2'], min_graph_size_list=['5'], \
		min_edge_freq_list=['6'], first_density_cutoff_list=['0.4'], second_density_cutoff_list=['0.4'],\
		max_pre_graph_size_list=['80'], conn_perc_list=['0.5'], match_cut='4.0', intersect2union_cut='0.4',\
		euclidean_ratio='0.2', selection_code='0110', max_degree='500', dataset_selected=None, type=0, \
		id=None, od=None, bd=None, mp='sc_54_6661_merge_6', op=None, rank_range=None, \
		js=0, debug=0):
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
		
		self.match_cut = match_cut
		self.intersect2union_cut = intersect2union_cut
		self.euclidean_ratio = euclidean_ratio
		self.selection_code = selection_code
		
		self.max_degree = max_degree
		self.dataset_selected = dataset_selected
		
		self.type = int(type)
		self.id = id
		self.od = od
		self.bd = bd
		self.mp = mp
		self.op = op
		self.rank_range = rank_range
		self.js = int(js)
		self.debug = int(debug)

	def parameter_list_init(self):
		"""
		04-01-05
			
		"""
		#None arguments
		self.homedir = os.path.expanduser("~")
		if self.id == None:
			#self.id = os.path.join(os.path.expanduser('~'), 'bin/hhu_clustering/data/input')
			self.id = os.path.join(self.homedir,'bin/hhu_clustering/data/input')
		if self.od == None:
			#self.od = os.path.join(os.path.expanduser('~'), 'bin/hhu_clustering/data/output')
			self.od = os.path.join(self.homedir, 'bin/hhu_clustering/data/output/netmine')
		if self.bd == None:
			#self.bd = os.path.join(os.path.expanduser('~'), 'bin/hhu_clustering/bin')
			self.bd = os.path.join(self.homedir, 'bin/hhu_clustering/bin')
		
		#additional class wide parameters
		self.input_matrix_file = os.path.join(self.id, self.mp+'.matrix')
		self.input_sv_file = os.path.join(self.id, self.mp+'.cor_vector')
		self.input_bv_file = os.path.join(self.id, self.mp+'.sig_vector')
		self.input_ttable_file = os.path.join(self.id, self.ttablefile)
		self.binary = 'netmine'
		self.binary2nd = 'netmine2nd'
		program_path = os.path.join(self.bd, self.binary)
		program_path2nd = os.path.join(self.bd, self.binary2nd)
		
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
										int_match_cut = int(float(self.match_cut))
										int_intersect2union_cut = int(float(self.intersect2union_cut)*100)
										int_euclidean_ratio = int(float(self.euclidean_ratio)*100)
										self.op = self.mp+'G%sE%sD%sQ%sS%sC%sH%sJ%sW%sZ%s'%\
		(min_graph_size, min_edge_freq, int_first_density_cutoff, int_second_density_cutoff, max_pre_graph_size, int_conn_perc,\
		int_match_cut, int_intersect2union_cut, int_euclidean_ratio, self.selection_code)
									#plus the directory
									op = os.path.join(self.od, self.op)
									#first is netmine's parameter
									parameter_list.append([program_path,\
		'-m', self.run_mode, '-i', self.input_matrix_file, '-n', self.genenum, \
		'-o', op, '-g', min_graph_size,\
		'-e', min_edge_freq, '-d', first_density_cutoff, '-q', second_density_cutoff, \
		'-s', max_pre_graph_size, '-c', conn_perc, '-u', self.max_degree, '-y', self.svnum, '-r', cut_loop_num])
									#second is netmine2nd's parameter
									parameter_list.append([program_path2nd, '-r', cut_loop_num,\
		'-m', self.run_mode, '-i', self.input_matrix_file, '-n', self.genenum, '-v', self.input_sv_file, '-a', self.input_bv_file,\
		'-p', self.svnum, '-l', self.sv_length, '-t', self.input_ttable_file, '-o', op, '-g', min_graph_size,\
		'-e', min_edge_freq, '-d', first_density_cutoff, '-q', second_density_cutoff, \
		'-s', max_pre_graph_size, '-c', conn_perc, '-h', self.match_cut, '-j', self.intersect2union_cut,\
		'-w', self.euclidean_ratio, '-z', self.selection_code, '-u', self.max_degree, '-y', self.svnum])
									if self.dataset_selected:	#not None	06-02-05
										parameter_list[1].append('-k')
										parameter_list[1].append(self.dataset_selected)

		return parameter_list


	def run_netmine(self, node, netmine_parameter_list):
		"""
		04-01-05
			input: wl_list of netmine
			output: number of clusters
			
			nothing to do with mpi
		05-16-05
			os.system() might fail on some platforms, so use recursive system_call()
		"""
		sys.stderr.write("Running netmine...")
		#wl = ['ssh', 'node%s'%node, '%s'%' '.join(netmine_parameter_list)]
		"""
		04-08-05 spawnvp gets dead under MPI, use system instead.
		"""
		#return_code = os.spawnvp(os.P_WAIT, netmine_parameter_list[0], netmine_parameter_list)
		commandline = '%s'%' '.join(netmine_parameter_list)
		if self.debug:
			sys.stderr.write("The commandline of netmine is %s\n"%commandline)
		exit_code = system_call(commandline)	#05-16-05 use the recursive one.
		#exit_code = os.system('%s'%' '.join(netmine_parameter_list))
		op = netmine_parameter_list[8]	#the 8th is the output file
		no_of_clusters = 0
		of = open(op, 'r')
		for line in of:
			no_of_clusters += 1
		sys.stderr.write("total clusters:%s. exit_code: %s.\n"%(no_of_clusters, exit_code))
		return no_of_clusters
	
	def schedule_netmine2nd(self, rank_range, netmine2nd_parameter_list, no_of_clusters):
		"""
		04-01-05
			non-mpi
		"""
		#hour:minute, minute is current_minute + 2, if time is like 11:61, it's fine.
		time_tuple = time.localtime()
		time_to_run_jobs = "%s:%s"%(time_tuple[3], time_tuple[4]+2)
		
		"""
		testing

		program_path = 'ls'
		parameter_list = [['ls'], ['ls', '/'], ['ls', '/usr/local']]
		"""
		
		#map clusters to each node
		node_rank2cluster_no = {}
		for i in range(no_of_clusters):
			#remainder is the node_rank
			index = int(math.fmod(i, len(rank_range)))
			node_rank = rank_range[index]
			if node_rank not in node_rank2cluster_no:
				node_rank2cluster_no[node_rank] = [i]
			else:
				node_rank2cluster_no[node_rank].append(i)
		job_number = self.js
		for node_rank in node_rank2cluster_no:
			#some nodes will be idle if there're more nodes than jobs
			job_fname = os.path.join(os.path.expanduser('~'), 'qjob/netmine_wrapper%s.sh'%job_number)
			job_f = open(job_fname, 'w')
			
			for cluster_no in node_rank2cluster_no[node_rank]:
				job_f.write('date\n')	#the beginning time
				jobrow = ['ssh', 'node%s'%node_rank]
				jobrow = jobrow + netmine2nd_parameter_list + ['-f', repr(cluster_no)]
				job_f.write('echo %s\n'%' '.join(jobrow))	#print the commandline
				job_f.write("%s\n"%' '.join(jobrow))	#command here
				job_f.write('date\n')	#the ending time
			
			#close the file
			job_f.close()
				
			print "node: %s, at %s, clusters: %s"%(node_rank, time_to_run_jobs, repr(node_rank2cluster_no[node_rank]))
			#schedule it
			#wl = ['at', '-mf', job_fname, time_to_run_jobs]	#-m, even if no output, mail me.
			#04-01-05 at has a problem. 'at' only executes the first 'ssh nodexx ....' and it replaces the whole 'at' job. return without
			#running the following 'ssh ...'
			sh_wl = ['sh', job_fname]
			os.spawnvp(os.P_NOWAIT, 'sh', sh_wl)
				
			job_number+=1

	def distribute_jobs(self, communicator, no_of_clusters, node_function, netmine2nd_parameter_list):
		"""
		04-09-05
			input: rank_range, no_of_clusters
			output: node_rank2cluster_no
			
			distribute jobs among the nodes in the rank_range based on no_of_clusters
		05-19-05
			(rewritten)
		"""
		job_list = range(no_of_clusters)
		if self.debug:
			sys.stderr.write("The common parameter_list of netmine2nd is %s.\n"%repr(' '.join(netmine2nd_parameter_list)))
		of_name_list = mpi_schedule_jobs(communicator, job_list, node_function, netmine2nd_parameter_list, self.debug)
		return of_name_list
		"""
		#map clusters to each node
		node_rank2cluster_no = {}
		for i in range(no_of_clusters):
			#remainder is the node_rank
			index = int(math.fmod(i, len(rank_range)))
			node_rank = rank_range[index]
			if node_rank not in node_rank2cluster_no:
				node_rank2cluster_no[node_rank] = [i]
			else:
				node_rank2cluster_no[node_rank].append(i)
		return node_rank2cluster_no
		"""
		
		
	def node_fire(self, communicator, node_rank2cluster_no, netmine2nd_parameter_list):
		"""
		04-08-05
			mpi form
		05-16-05
			os.system() might fail on some platforms, so use recursive system_call()
		"""
		#get the output filename prefix
		ofname_prefix = netmine2nd_parameter_list[20]
		#store the finished output files and send it back to node 0 later
		ofname_list = []
		if communicator.rank in node_rank2cluster_no:	#check if it's necessary to ignite the node
			#fire repeatedly
			for cluster_no in node_rank2cluster_no[communicator.rank]:
				jobrow = netmine2nd_parameter_list + ['-f', repr(cluster_no)]
				sys.stderr.write("node %s working on cluster %s...\n"%(communicator.rank, cluster_no))
				exit_code = system_call("%s"%' '.join(jobrow))
				sys.stderr.write("node %s on cluster %s done, exit_code: %s\n"%(communicator.rank, cluster_no, exit_code))
				ofname_list.append('%s_%sh'%(ofname_prefix, cluster_no))	#04-09-05 haiyan puts 'h' after cluster_no.??
		
			#join the list by blank and send it back to node 0
			communicator.send("%s"%' '.join(ofname_list),0,communicator.rank)
	
	def collect_and_merge_output(self, of_name_list, final_ofname):
		"""
		04-08-05
			only for rank 0
			1. receive the concatenated ofname_list from all nodes
			2. cat all files together into final_ofname and delete the intermediary files.
		05-16-05
			os.system() might fail on some platforms, so use recursive system_call()
		05-19-05
			collecting output filenames is NOT necessary, it's done in distribute_jobs()
				(via mpi_schedule_jobs())
		"""
		sys.stderr.write("Catting them together...")
		#05-22-05 remove the final_ofname first to avoid '>>' to add up
		if os.path.isfile(final_ofname):
			os.remove(final_ofname)
		
		for ofname in of_name_list:
			#do it one by one because 'cat * >final_ofname' might cause error 'Argument list too long'
			#which is the case for 'rm' when there're thousands of files represented by '*'.
			exit_code = system_call("cat %s >> %s"%(ofname, final_ofname))	#it's >> not >
			exit_code = system_call("rm %s"%ofname)	#delete it immediately
		sys.stderr.write("Done.\n")
		
	def run(self):
		"""
		04-01-05
		
		05-19-05
			--parameter_list_init()
			--run_netmine()
			
			--distribute_jobs()
				--run_netmine2nd()
			--collect_and_merge_output()
		"""
		communicator = MPI.world.duplicate()
		rank_range = range(communicator.size)
		parameter_list = self.parameter_list_init()
		no_of_clusters_broadcast  = Numeric.zeros((1,), Numeric.Int)
		if communicator.rank == 0:
			sys.stderr.write("this is node %s\n"%communicator.rank)
			no_of_clusters = self.run_netmine(rank_range[0], parameter_list[0])
			no_of_clusters_broadcast[0] = no_of_clusters
		#let every node know no_of_clusters
		communicator.broadcast(no_of_clusters_broadcast, 0)
		
		mpi_synchronize(communicator)
		
		#ignite each node
		of_name_list = self.distribute_jobs(communicator, no_of_clusters_broadcast[0], run_netmine2nd, parameter_list[1])
		
		mpi_synchronize(communicator)
		
		#collecting
		if communicator.rank==0:
			final_ofname = os.path.join(self.od, 'F'+self.op)
			self.collect_and_merge_output(of_name_list, final_ofname)
				#receive only from those fired nodes
		
		#self.schedule_netmine2nd(self.rank_range, parameter_list[1], no_of_clusters)
		
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	try:
		opts, args = getopt.getopt(sys.argv[1:], "m:n:p:l:t:r:g:e:d:q:s:c:h:j:w:z:u:k:y:", ["help","run_mode=",\
			"genenum=", "svnum=", "sv_length=", "ttablefile=", "cut_loop_num_list=", "min_graph_size_list=",\
			"min_edge_freq_list=", "first_density_cutoff_list=", "second_density_cutoff_list=", \
			"max_pre_graph_size_list=", "conn_perc_list=", "match_cut=", "intersect2union_cut=",\
			"euclidean_ratio=", "selection_code=", "max_degree=", "type=", "id=", "od=", \
			"bd=", "mp=", "op=", "rr=", "js=", "debug"])
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
	
	match_cut = '4.0'
	intersect2union_cut = '0.4'
	euclidean_ratio = '0.2'
	selection_code = '0110'
	
	max_degree = '500'
	dataset_selected = None
	
	type = '0'
	id = None
	od = None
	bd = None
	mp = None
	op = None
	rank_range = []
	js = 0
	debug = 0
	
	for opt, arg in opts:
		if opt in ("--help"):
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
		
		elif opt in ("-h", "--match_cut"):
			match_cut = arg
		elif opt in ("-j", "--intersect2union_cut"):
			intersect2union_cut = arg
		elif opt in ("-w", "--euclidean_ratio"):
			euclidean_ratio = arg
		elif opt in ("-z", "--selection_code"):
			selection_code = arg
			
		elif opt in ("-u", "--max_degree"):
			max_degree = arg
		elif opt in ("-k"):
			dataset_selected = arg
			
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
		elif opt in ("--debug"):
			debug = 1
	if mp:
		instance = netmine_wrapper(run_mode, genenum, svnum, sv_length,\
			ttablefile, cut_loop_num_list, min_graph_size_list, \
			min_edge_freq_list, first_density_cutoff_list, second_density_cutoff_list,\
			max_pre_graph_size_list, conn_perc_list, match_cut, intersect2union_cut,\
			euclidean_ratio, selection_code, max_degree, dataset_selected, type,\
			id, od, bd, mp, op, rank_range, js, debug)
		
		instance.run()
		
	else:
		print __doc__
		sys.exit(2)
