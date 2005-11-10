#!/usr/bin/env python
"""
Usage: p_gene_lm.py -k SCHEMA -t TABLE -s SPLAT_TABLE [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-f ..., --netmine_fname=...	used to construct table names
		if netmine_fname is given, no need to given following four tables
	-t ..., --table=...	the p_gene table
	-s ..., --splat_table=...	the corresponding splat_table
	-m ..., --mcl_table=...	the corresponding mcl_table
	-l ,,,, --lm_table=...	the lm_table to store the linear_model results, needed if needcommit
	-a ..., --accuracy_cut_off=...	0.5(default)
	-p ..., --percentage=...	0.3(default) if accuracy_cut_off=0, use percentage to select score_cut_off
	-j ..., --judger_type=...	how to judge predicted functions, 0(default), 1, 2
	-b ..., --bit_string=...	the bit_string control which parameter to be counted into regression
		p_value, recurrence, connectivity, cluster_size, connectivity_2nd;	11111(default)
	-n, --no_redundancy	use the pair accuracy to find score cutoff
	-o, --logistic	logistic regression
	-c, --commit	commit this database transaction
	-r, --report	report flag
	-u, --debug debug flag
	-h, --help              show this help

Examples:
	p_gene_lm.py -k sc_54 -t p_gene_repos_2_e5 -s splat_table  -j 1 -r >/tmp/p_gene_lm.out
	p_gene_lm.py -k sc_54 -t p_gene_repos_2_e5 -s splat_repost
		-l p_gene_repos_2_e5_lm -j 1 -r -c >/tmp/p_gene_lm.out
	p_gene_lm.py -k mm_oxi_stress_7t1 -f fmos_7t1g1e3d40q20s200c50z0001c8 -j 2 >/tmp/p_gene_lm.out
	
Description:
	linear model:
		is_correct ~ p_value + recurrence + connectivity + cluster_size + connectivity_2nd
	
	Output contains the p_value_cut_off and other information for each go-no
	as well as the linear_model fitting results. Output is dumped on sys.stdout.
	
	bit_string can also be less than 4 digits. The leftout digit is regarded as 0.
"""

import sys, os, psycopg, getopt, csv, math
from p_gene_analysis import prediction_space_attr
from codense.common import *
from numarray import *
from Numeric import average
from rpy import r, set_default_mode,NO_CONVERSION,BASIC_CONVERSION

class p_gene_lm:
	"""
	02-28-05
		get the data from p_gene table
		for each go, for each (recurrence, connectivity), compute its p_value_cut_off based on accuracy_cut_off
		for each go, linear-model fitting of p_value_cut_off ~ recurrence + connectivity
		submit the coefficients and chisq of linear-model to a lm_table.
		
		parameter min_data_points (no of known predictions) controls the validity of each
		prediction space(recurrence, connectivity).
		
		Some go-nos don't have enough data to do linear-model fitting, simply ignore. But lm_table
		will contain an average model for these go-nos to look up in the later prediction stage.
	
	03-27-05
		call R to use logistic regression. No floor taking of recurrence and connectivity.
	04-05-05
		use connectivity in the corresponding splat_table to do logistic regression(previously use
		my definition which is in the p_gene table copied from mcl_table)
	06-30-05
		add corresponding p-values, parameter cluster_size,
		add bit_string
	"""
	def __init__(self, hostname=None, dbname=None, schema=None, netmine_fname=None, table=None, splat_table=None,\
		mcl_table=None, lm_table=None, accuracy_cut_off=0, percentage=0.3, judger_type=0, bit_string='11111', min_data_points=5, \
		no_redundancy=0, logistic=0, needcommit=0, report=0, debug=0, valid_space=20, recurrence_gap_size=2, connectivity_gap_size=2):
		"""
		03-08-05
			add two more parameters, recurrence_gap_size and connectivity_gap_size (make them explicit)
		06-30-05
			add bit_string
			add percentage
		07-06-05
			add logistic
		08-08-05
			add no_redundancy
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.netmine_fname = netmine_fname
		self.table = table
		self.splat_table = splat_table
		self.mcl_table = mcl_table
		
		self.lm_table = lm_table
		self.accuracy_cut_off = float(accuracy_cut_off)
		self.percentage = float(percentage)
		self.judger_type = int(judger_type)
		self.bit_string = bit_string
		self.min_data_points = int(min_data_points)
		self.no_redundancy = int(no_redundancy)
		self.logistic = int(logistic)
		self.needcommit = int(needcommit)
		self.report = int(report)		
		self.debug = int(debug)
		self.valid_space = int(valid_space)
		#the gap between two recurrences
		self.recurrence_gap_size = int(recurrence_gap_size)
		self.connectivity_gap_size = int(connectivity_gap_size)
		
		self.go_no2prediction_space = {}
		self.mcl_id2prediction_space = {}
		#an is_correct dictionary used in database fetch
		self.is_correct_dict = {0: 'is_correct',
			1: 'is_correct_L1',
			2: 'is_correct_lca'}
		

		
	def init(self):
		"""
		03-27-05
			use rpy to call r instead, the gsl linear_model is useless.
		"""
		self.lm_instance = None
		self.no_of_zero_p_values = 0
		
	def data_fetch(self, curs, table):
		"""
		02-28-05
			borrowed from p_gene_analysis.py
		06-30-05
			fetch cluster_size_cut_off
		07-04-05
			fetch m.connectivity(2nd-order density)
		08-15-05
			regard 1e7 as the memory threshold, use step to select entries
		10-10-05
			temporarily modify something
		10-26-05 get connectivity from p_gene_table, not splat_table
			
			--prediction_space_setup()
		"""
		sys.stderr.write("Setting Up prediction_space...\n")
		
		#08-15-05	step to control memory usage
		curs.execute("select count(*) from %s"%table)
		rows = curs.fetchall()
		no_of_entries = rows[0][0]
		step = int(no_of_entries/1e7) + 1
		sys.stderr.write("\tstep is %s.\n"%step)
		
		curs.execute("DECLARE crs CURSOR FOR select p.gene_no, p.go_no, p.mcl_id, p.%s, p.avg_p_value, \
			p.recurrence_cut_off,p.connectivity_cut_off, p.depth_cut_off,p.cluster_size_cut_off,p.edge_gradient from %s p"\
			%(self.is_correct_dict[self.judger_type], table))
		no_of_records = 0
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		while rows:
			for row in rows:
				if (no_of_records%step) == 0:
					self.prediction_space_setup(row)
				no_of_records += 1
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, no_of_records))
			
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		sys.stderr.write("\t%d zero p_values converted to 1e-8.\n"%(self.no_of_zero_p_values))
		sys.stderr.write("done.\n")

	def prediction_space_setup(self, row):
		"""
		02-28-05
			borrowed from _p_gene_analysis() of p_gene_analysis.py
		06-30-05
			add cluster_size
		07-04-05
			add connectivity_2nd
		07-05-05
			prediction_space goes to mcl_id2prediction_space
		07-25-05
			prediction_space goes back to go_no2prediction_space
		08-08-05
			add gene_no and go_no to control accuracy pair
		"""
		gene_no = row[0]
		go_no = row[1]
		mcl_id = row[2]
		is_correct = row[3]
		p_value = row[4]
		recurrence = row[5]
		connectivity = row[6]
		depth_cut_off = row[7]
		cluster_size = row[8]
		connectivity_2nd = row[9]
		
		#take the floor of the recurrence
		#recurrence = int(math.floor(recurrence/self.recurrence_gap_size)*self.recurrence_gap_size)
		#take the floor of the connectivity *10
		#connectivity = int(math.floor(connectivity*10/self.connectivity_gap_size)*self.connectivity_gap_size)
		
		if is_correct==-1:
			#unknown genes, not for model training
			return
		"""
		03-27-05
			Don't distinguish go_no.
			Assign go_no = -1.
		"""
		#convert p_value to -log(p_value)
		"""10-23-05 no need to convert p-value
		if p_value == 0:
			self.no_of_zero_p_values += 1
			p_value = -math.log(1e-8)
		else:
			p_value = -math.log(p_value)
		"""
		unit = [p_value, recurrence, connectivity, cluster_size, connectivity_2nd, gene_no, go_no, is_correct]	#add gene_no and go_no to control accuracy pair
		#if mcl_id not in self.mcl_id2prediction_space:
		#	self.mcl_id2prediction_space[mcl_id] = []
		#self.mcl_id2prediction_space[mcl_id].append(unit)
		go_no = -1
		if go_no not in self.go_no2prediction_space:
			self.go_no2prediction_space[go_no] = []
		self.go_no2prediction_space[go_no].append(unit)

	def mcl_id2prediction_spaceTogo_no2prediction_space(self, mcl_id2prediction_space):
		"""
		07-05-05
			transform mcl_id2prediction_space to go_no2prediction_space, each mcl_id's prediction_space
			is averaged
		"""
		sys.stderr.write("Transforming mcl_id2prediction_space...")
		go_no2prediction_space = {}
		go_no = -1
		go_no2prediction_space[go_no] = []
		for (mcl_id, data) in mcl_id2prediction_space.iteritems():
			data = array(data)
			unit = [average(data[:,0]), average(data[:,1]), average(data[:,2]), average(data[:,3]), average(data[:,4]), average(data[:,-1])]
			if self.debug:
				print data
				print "averaged",unit
				is_continue = raw_input("Continue?Y/n")
				if is_continue =='n':
					sys.exit(2)
			go_no2prediction_space[go_no].append(unit)
		sys.stderr.write("%s data points to do fitting. Done.\n"%len(go_no2prediction_space[go_no]))
		return go_no2prediction_space

	def lm_fit(self, lm_instance, go_no2prediction_space, bit_string, curs=None, lm_table=None):
		"""
		02-28-05
			linear model fitting here
		
		03-08-05
			grouping and accumulating before do linear model fitting, see log of 2005, 
			section 'linear model overfitting' for detail.
		03-27-05
			Use glm of R to do logistic regression
		06-30-05
			add cluster_size
			add bit_string to control which parameter should be enabled.
		07-04-05
			add connectivity_2nd
		07-06-05
			add logistic
		11-09-05 extend coeff_list and coeff_p_value_list
			restructure the list, go_no2lm_results[go_no]
			
			--data_prepare
			--submit
		"""
		sys.stderr.write("Linear Model Fitting...\n")
		go_no2lm_results = {}
		
		#06-30-05	setup the formula_list based on bit_string
		coeff_name_list = ['p_value', 'recurrence', 'connectivity', 'cluster_size', 'connectivity_2nd']
		formula_list = []
		for i in range(len(bit_string)):
			if bit_string[i] == '1':
				formula_list.append(coeff_name_list[i])
		
		for (go_no,data) in go_no2prediction_space.iteritems():
			sys.stderr.write("%s prediction entries from %s.\n"%(len(data), go_no))
			#11-09-05 extend coeff_list and coeff_p_value_list
			coeff_list = [0]*7	#intercept, p_value, recurrence, connectivity, cluster_size
			coeff_p_value_list = [1]*7
			index = 0	#06-30-05	the pointer for summary_stat
			
			if len(data)<=50:
				#two few data
				continue
			#convert it to a 2d array
			data = array(data)
			"""
			data_frame = r("d=data.frame(p_value=c(%s),recurrence=c(%s),connectivity=c(%s), is_correct=c(%s))"%(repr(list(data[:,0]))[1:-1], \
				repr(list(data[:,1]))[1:-1], repr(list(data[:,2]))[1:-1], repr(list(data[:,3]))[1:-1]))
			lm_result = r("lm_result=glm(is_correct~p_value+recurrence+connectivity, data=d,family=binomial)")
			significance_dict = r("summary(lm_result)")
			print significance_dict['coefficients']
			"""
			set_default_mode(NO_CONVERSION) #04-07-05
			data_frame = r.as_data_frame({"p_value":data[:,0], "recurrence":data[:,1], "connectivity":data[:,2], \
				"cluster_size":data[:,3], "connectivity_2nd":data[:,4], "is_correct":data[:,-1]})	#06-30-05	-1 denotes is_correct
			if self.logistic:
				lm_result = r.glm(r("is_correct~%s"%'+'.join(formula_list)), data=data_frame, family=r("binomial"))
			else:
				lm_result = r.glm(r("is_correct~%s"%'+'.join(formula_list)), data=data_frame)	#06-30-05 use formula_list
			set_default_mode(BASIC_CONVERSION) #04-07-05
			#04-07-05 r.summary() requires lm_result in NO_CONVERSION state
			summary_stat = r.summary(lm_result)
			if self.debug:
				print "everything about coefficients from function", go_no, "is"
				print summary_stat['coefficients']	#p-values of coefficients
			"""
			#04-07-05 convert to python dictionary form
			lm_result = lm_result.as_py()
			coeff_list = [lm_result["coefficients"]["(Intercept)"], lm_result["coefficients"]["p_value"], \
				lm_result["coefficients"]["recurrence"], lm_result["coefficients"]["connectivity"], \
				lm_result["coefficients"]["cluster_size"], \
				summary_stat['coefficients'][0][-1], summary_stat['coefficients'][1][-1],\
				summary_stat['coefficients'][2][-1], summary_stat['coefficients'][3][-1],\
				summary_stat['coefficients'][4][-1], 1]
				#the last entry is score_cut_off, replaced later in get_score_cut_off()
				#06-30-05	add corresponding p-values
			"""
			#06-30-05	0 in summary_stat['coefficients'] is intercept
			coeff_list[0] = summary_stat['coefficients'][0][0]	#0 is the coefficient
			coeff_p_value_list[0] = summary_stat['coefficients'][0][-1]	#-1 is the corresponding p-value
			#06-30-05	fill in other efficients based on bit_string, NOTE i+1
			for i in range(len(bit_string)):
				if bit_string[i] == '1':
					index+=1
					coeff_list[i+1] = summary_stat['coefficients'][index][0]	#0 is the coefficient
					coeff_p_value_list[i+1] = summary_stat['coefficients'][index][-1]	#-1 is the corresponding p-value
			#11-09-05 restructure the following list
			go_no2lm_results[go_no] = [coeff_list, coeff_p_value_list, 1]	#the last entry is score_cut_off, replaced later in get_score_cut_off()
		sys.stderr.write("done.\n")
		return go_no2lm_results
	
	def lm_results_output(self, outf, go_no2lm_results):
		"""
		02-28-05
			output the go_no2lm_results
		06-30-05
			change the header of the output.
		"""
		sys.stderr.write("Outputting Linear Model parameters...")
		writer = csv.writer(outf, delimiter='\t')
		writer.writerow(['go_no', 'intercept', 'coefficients and p-values', 'score_cut_off'])
		for (go_no, coeff_list) in go_no2lm_results.iteritems():
			writer.writerow([go_no]+coeff_list)
		del writer
		sys.stderr.write("done.\n")
	
	def get_score_cut_off(self, go_no2prediction_space, go_no2lm_results, accuracy_cut_off, percentage):
		"""
		03-27-05
			input: go_no2prediction_space, go_no2lm_results
			output: score_cut_off
			
			After model fitting, every prediction has a score computed from the model.
			Based on the accuracy_cut_off, this function computes the corresponding
			score_cut_off.
		06-30-05
			add cluster_size
			add percentage
		07-04-05
			add connectivity_2nd
		08-08-05
			add gene_no and go_no to score_list
		11-09-05 extend coeff_list and coeff_p_value_list
			restructure the list, go_no2lm_results[go_no]
			this is for OneParameterCutoffSeeker.py. Itself needs testing.
			
			--return_score_cut_off()
		"""
		sys.stderr.write("Getting score cutoff for accuracy_cut_off %s...\n"%(accuracy_cut_off))
		for go_no,data in go_no2prediction_space.iteritems():
			if go_no not in go_no2lm_results:
				#this go_no has too few data, ignored.
				continue
			score_list = []
			#11-09-05 restructure the list, go_no2lm_results[go_no]
			coeff_list, coeff_p_value_list, score_cut_off = go_no2lm_results[go_no]
			for entry in data:
				#intercept + coeff1*p_value + coeff2*recurrence + coeff3*connectivity + coeff4*cluster_size + coeff5*connectivity_2nd
				score = coeff_list[0]+ coeff_list[1]*entry[0] + coeff_list[2]*entry[1] + coeff_list[3]*entry[2] + coeff_list[4]*entry[3] + coeff_list[5]*entry[4]
				#score, is_correct
				score_list.append([score, entry[-1], entry[-3], entry[-2]])	#06-30-05	-1 denotes is_correct, 08-08-05 -3,-2 are gene_no,go_no.
			if accuracy_cut_off==0:
				score_cut_off = self.return_score_cut_off_by_percentage(score_list, percentage, go_no)
			else:
				score_cut_off = self.return_score_cut_off(score_list, accuracy_cut_off, go_no)
			if score_cut_off:
				#found the cutting point, append the score cutoff to the coeff_list
				go_no2lm_results[go_no][-1] = score_cut_off
			else:
				#can't find a score_cut_off to meet the accuracy_cut_off, append a maximum
				go_no2lm_results[go_no][-1] = 1e10
		sys.stderr.write("Done.\n")
		return go_no2lm_results
	
	def return_score_cut_off(self, score_list, accuracy_cut_off, go_no):
		"""
		03-27-05
		08-08-05
			1. use the bicut method to speed up the score_cut_off searching
			2. judge self.no_redundancy
				get_pair_accuracy()
		"""
		sys.stderr.write("\tReturning score_cut_off for go %d..."%(go_no))
		#default no score_cut_off
		score_cut_off = None
		score_list.sort()
		if self.debug:
			print "score\tis_correct"
			for score_tuple in score_list:
				print "%s\t%s"%(score_tuple[0], score_tuple[1])
		#convert to a 2d array, they have the same indices
		score_array = array(score_list)
		previous_score = None
		no_of_known_predictions = len(score_array)
		if self.debug:
			print "score\tis_correct_array\ttotal\taccuracy"
		#08-08-05 	use the bicut method to speed up
		gap = 1
		lower = 0
		upper = len(score_array)
		list_index = upper/2
		while gap>0.001:
			score = score_array[list_index,0]
			previous_list_index = list_index
			if self.no_redundancy:	#08-08-05
				accuracy = self.get_pair_accuracy(score_array[list_index:])
			else:
				correct_array = score_array[:,1][list_index:]
				accuracy = sum(correct_array)/float(len(correct_array))
				if self.debug:
					sys.stderr.write("correct_array(%s, length): %s\n"%(len(correct_array), repr(correct_array)))
			if self.debug:
				print "list_index: %s; %s\t%s"%(list_index, score, accuracy)
				raw_input("Continue:(Y/n)?")
			gap = abs(accuracy-accuracy_cut_off)
			if accuracy>=accuracy_cut_off:
				upper = list_index
				list_index = (list_index + lower)/2
			else:
				lower = list_index
				list_index = (upper + list_index)/2
			if list_index == previous_list_index:	#convergence failure
				break
		
		score_cut_off = score
		"""
		#calculate the accuracy from the low score cutoff to high score cutoff, once >= accuracy_cut_off, break
		for i in range(len(score_array)):
			score = score_array[i,0]
			if score == previous_score:
				#same score, previous one doesn't meet the cutoff, this one won't either.
				continue
			else:
				#this score becomes previous_score for the next score
				previous_score = score
				if self.no_redundancy:	#08-08-05
					accuracy = self.get_pair_accuracy(score_array[i:])
				else:
					#count the is_correct==1 entries from i-th row, and divide it by (no_of_known_predictions-i+1)
					correct_array = score_array[:,1][i:]
					accuracy = sum(correct_array)/float(len(correct_array))
					if self.debug:
						sys.stderr.write("correct_array(%s, length): %s\n"%(len(correct_array), repr(correct_array)))
				if self.debug:
					print "%s\t%s"%(score, accuracy)
					raw_input("Continue:(Y/n)?")
				if accuracy >= accuracy_cut_off:
					score_cut_off = score
					break
		"""
		sys.stderr.write("Done.\n")
		return score_cut_off
	
	def get_pair_accuracy(self, partial_score_array):
		"""
		08-08-05
			get the pair accuracy which is the non-redundant accuracy
		"""
		gene_no_go_no2correct = {}
		for entry in partial_score_array:
			score, is_correct, gene_no, go_no = entry
			gene_no_go_no2correct[(gene_no, go_no)] = is_correct
		accuracy = sum(gene_no_go_no2correct.values())/float(len(gene_no_go_no2correct))
		if self.debug:
			sys.stderr.write("No of gene_no:go_no pairs:%s; correct: %s\n"%\
				(len(gene_no_go_no2correct), sum(gene_no_go_no2correct.values())))
		return accuracy
	
	def return_score_cut_off_by_percentage(self, score_list, percentage, go_no):
		"""
		06-30-05
		"""
		sys.stderr.write("\tReturning score_cut_off for go %s based on percentage %s...\n"%(go_no, percentage))
		#default no score_cut_off
		score_cut_off = None
		score_list.sort()
		score_list.reverse()	#####different from return_score_cut_off. descending order
		cut_off_index = int(len(score_list)*percentage)
		if cut_off_index!=0:
			cut_off_index -= 1	#index starts from 0, not like counting
		score_cut_off = score_list[cut_off_index][0]
		#convert to a 2d array, they have the same indices
		score_array = array(score_list)
		correct_array = score_array[:,1][:cut_off_index]
		accuracy = sum(correct_array)/float(len(correct_array))
		sys.stderr.write("\tThe rough accuracy for percentage %s is %s.\n"%(percentage, accuracy))
		sys.stderr.write("Done.\n")
		return score_cut_off
		
	def lm_table_create(self, curs, lm_table):
		"""
		03-27-05
			coeff_list has changed, intercept+3 coefficients+score_cut_off
		06-30-05
			add  fields for p-values, cluster_size(coeff4)
		07-04-05
			add connectivity_2nd and its p-value
		11-09-05 add coeff6
		"""
		sys.stderr.write("Creating Linear Model table...")
		try:
			curs.execute("create table %s(\
				go_no	integer,\
				intercept	float,\
				coeff1	float,\
				coeff2	float,\
				coeff3	float,\
				coeff4	float,\
				coeff5	float,\
				coeff6	float,\
				intercept_p_value	float,\
				coeff1_p_value	float,\
				coeff2_p_value	float,\
				coeff3_p_value	float,\
				coeff4_p_value	float,\
				coeff5_p_value	float,\
				coeff6_p_value	float,\
				score_cut_off	float)"%lm_table)
		except:
			sys.stderr.write("Error occurred when creating table %s\n"%lm_table)
			sys.exit(3)
		sys.stderr.write("done.\n")

	def submit(self, curs, lm_table, go_no2lm_results):
		"""
		02-28-05
			submit the linear model parameters to the database
		03-27-05
			coeff_list has changed, intercept+3 coefficients+score_cut_off
		06-30-05
			add corresponding p-values and coeff4(cluster_size)
		07-04-05
			add connectivity_2nd and restructure the coeff, p_value, score_cut_off list
		11-09-05 add coeff6
		"""
		sys.stderr.write("Submitting Linear model parameters...")
		for (go_no, coeff_p_value_score_list) in go_no2lm_results.iteritems():
			coeff_list, p_value_list, score_cut_off = coeff_p_value_score_list
			curs.execute("insert into %s(go_no, intercept, coeff1, coeff2, coeff3, coeff4, coeff5, coeff6, \
				intercept_p_value, coeff1_p_value, coeff2_p_value, coeff3_p_value, coeff4_p_value, coeff5_p_value, \
				coeff6_p_value, score_cut_off) values (%d, %s,%s,%s,%s,%s,%s,%s,   %s,%s,%s,%s,%s,%s,%s,  %s)"%\
				(lm_table, go_no, coeff_list[0], coeff_list[1], coeff_list[2], coeff_list[3], coeff_list[4], coeff_list[5], coeff_list[6], \
				p_value_list[0], p_value_list[1], p_value_list[2], p_value_list[3], p_value_list[4], p_value_list[5], p_value_list[6], score_cut_off) )
		sys.stderr.write("done.\n")

	def run(self):
		"""
		02-24-05
		
		07-05-05
			add mcl_id2prediction_spaceTogo_no2prediction_space
			
			--init()
			--db_connect()
			--lm_table_create()
			--data_fetch()
				--prediction_space_setup()
			--lm_fit()
			--lm_results_output()
			--get_score_cut_off()
				--return_score_cut_off()
					--get_pair_accuracy()
			--lm_results_output()
			--submit()
		
		"""
		#some additional initialization
		self.init()
		
		(conn, curs) = db_connect(self.hostname, self.dbname, self.schema)
		curs.execute("begin")
		if self.needcommit:
			if self.lm_table:
				self.lm_table_create(curs, self.lm_table)
			else:
				sys.stderr.write("Please the lm_table to commit.\n")
				sys.exit(127)
			
		self.data_fetch(curs, self.table)
		#self.go_no2prediction_space = self.mcl_id2prediction_spaceTogo_no2prediction_space(self.mcl_id2prediction_space)	#07-25-05 discard this step
		go_no2lm_results = self.lm_fit(self.lm_instance, self.go_no2prediction_space, self.bit_string)
		self.lm_results_output(sys.stdout, go_no2lm_results)
		go_no2lm_results = self.get_score_cut_off(self.go_no2prediction_space, go_no2lm_results, self.accuracy_cut_off, self.percentage)
		self.lm_results_output(sys.stdout, go_no2lm_results)
		if self.needcommit:
			self.submit(curs, self.lm_table, go_no2lm_results)
			curs.execute("end")

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:f:t:s:m:l:a:p:j:b:nocru", ["help", "hostname=", \
			"dbname=", "schema=", "netmine_fname=", "table=", "splat_table=", "mcl_table=",  "lm_table=", \
			"accuracy_cut_off=", "percentage=", "judger_type=", "bit_string=", "no_redundancy", \
			"logistic", "commit", "report", "debug"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	netmine_fname = None
	table = None
	splat_table = None
	mcl_table = None
	lm_table = None
	accuracy_cut_off = 0.5
	percentage = 0.3
	judger_type = 0
	bit_string = '11111'
	no_redundancy = 0
	logistic = 0
	commit = 0
	report = 0
	debug = 0
	
	min_data_points = 5
	valid_space = 20
	recurrence_gap_size = 2
	connectivity_gap_size = 2

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
		elif opt in ("-f", "--netmine_fname"):
			netmine_fname = arg
		elif opt in ("-t", "--table"):
			table = arg
		elif opt in ("-s", "--splat_table"):
			splat_table = arg
		elif opt in ("-m", "--mcl_table"):
			mcl_table = arg
		elif opt in ("-l", "--lm_table"):
			lm_table = arg
		elif opt in ("-a", "--accuracy_cut_off"):
			accuracy_cut_off = float(arg)
		elif opt in ("-p", "--percentage"):
			percentage = float(arg)
		elif opt in ("-j", "--judger_type"):
			judger_type = int(arg)
		elif opt in ("-b", "--bit_string"):
			bit_string = arg
		elif opt in ("-n", "--no_redundancy"):
			no_redundancy = 1
		elif opt in ("-o", "--logistic"):
			logistic = 1
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-r", "--report"):
			report = 1
		elif opt in ("-u", "--debug"):
			debug = 1

	if netmine_fname:
		table = 'p_gene_%s_e5'%netmine_fname
		splat_table = 'splat_%s'%netmine_fname
		mcl_table = 'mcl_%s'%netmine_fname
		lm_table = 'lm_%s_e5'%netmine_fname
		
	if schema and table and splat_table and mcl_table:
		instance = p_gene_lm(hostname, dbname, schema, netmine_fname, table, splat_table, \
			mcl_table, lm_table, accuracy_cut_off, percentage, judger_type, bit_string, min_data_points, \
			no_redundancy, logistic, commit, report, debug, valid_space, recurrence_gap_size, connectivity_gap_size)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
