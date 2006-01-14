#!/usr/bin/env python
"""
Usage: DrawRpartValidationCurves.py -i -o

Option:
	-i ...,	inputfname
	-o ...,	outputfname prefix
	-p ...,	parameter list: cp,loss_ratio,prior(0.01,-1,0.5, default)
		for loss_ratio, use a(2nd in loss_matrix)/b(3rd) to denote a ratio
	-b,	enable debug flag
	-r,	enable report flag
	-h,	Display the usage infomation.
	
Examples:
	DrawRpartValidationCurves.py -i MpiRpartValidation.py_output
	-o /tmp/yuhuang/test -p -1,1/10,0.5

Description:
	Program to draw acc2cp and no_of_genes2cp plots.
	In parameter_list, -1 means variant, only one parameter could be -1.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path += [os.path.expanduser('~/lib64/python')]
else:   #32bit
	sys.path += [os.path.expanduser('~/lib/python')]
import sys, os, getopt, csv, math
from rpy import r

class DrawRpartValidationCurves:
	def __init__(self, inputfname=None, outputfname=None, parameter_list=[0.01, -1, 0.5], debug=0, report=0):
		self.inputfname = inputfname
		self.outputfname = outputfname
		self.parameter_list = parameter_list
		self.debug = int(debug)
		self.report = int(report)
		
		self.parameter_index2label = {0:'cp(log)',
			1:'loss_ratio',
			2:'prior'}
	
	def plot(self, outputfname, fix_index_ls, parameter_list, var_index, variant_ls, parameter_index2label, y_axis_ls, y_label):
		outputfname = '%s_%s2%s.png'%(outputfname, y_label, parameter_index2label[var_index])
		sys.stderr.write('Plotting %s'%outputfname)
		r.png(outputfname)
		r.plot(variant_ls, y_axis_ls, main='%s vs %s (%s=%s, %s=%s)'%(y_label, parameter_index2label[var_index],\
			parameter_index2label[fix_index_ls[0]], parameter_list[fix_index_ls[0]], parameter_index2label[fix_index_ls[1]],\
			parameter_list[fix_index_ls[1]]), xlab=parameter_index2label[var_index], ylab=y_label)
		r.dev_off()
		sys.stderr.write('Done.\n')
		
	
	def is_row_wanted(self, row, fix_index_ls, parameter_list):
		for index in fix_index_ls:
			if row[index]!=parameter_list[index]:
				return False
		return True
	
	def get_data_structure(self, inputfname, parameter_list):
		sys.stderr.write("Getting data structures...")
		fix_index_ls = []
		var_index = None
		for i in range(len(parameter_list)):
			if parameter_list[i] != -1:
				fix_index_ls.append(i)
			else:
				var_index = i
		
		if not var_index and len(fix_index_ls)!=len(parameter_list)-1:
			sys.stderr.write("parameter_list: %s error.\n"%parameter_list)
			sys.stderr.write("I got fix_index_ls: %s and var_index %s.\n"%(fix_index_ls, var_index))
			sys.exit(3)
		
		reader = csv.reader(open(inputfname, 'r'), delimiter='\t')
		reader.next()	#skip the header(1st line)
		variant_ls = []	#cp or loss_ratio or prior,	the x-axis in the plot
		data_dict = {'testing':[[], [], []],	#1st [] is for accuracy, 2nd [] is for no_of_predictions, 3rd [] is for no_of_genes
			'training':[[], [], []],
			'unknown_valid': [[], [], []],
			'unknown':[[], [], []],
			'known':[[], [], []]}
		for row in reader:
			#rpart_cp        loss_matrix     prior_prob      type    accuracy_avg    accuracy_std    
			#no_of_predictions_avg   no_of_predictions_std   no_of_genes_avg no_of_genes_std
			row[0] = float(row[0])
			loss_matrix = row[1]
			loss_matrix = loss_matrix[1:-1].split(',')
			loss_matrix = map(float, loss_matrix)
			loss_ratio = loss_matrix[1]/loss_matrix[2]
			row[1] = loss_ratio
			if row[2]:	#prior_prob might be None
				row[2] = float(row[2])
			else:
				row[2] = 0.0
			type = row[3]
			if self.is_row_wanted(row, fix_index_ls, parameter_list):
				row[0] = math.log(row[0], 10)	#log the cp
				if type == 'testing':	#to avoid redundant variant values
					variant_ls.append(row[var_index])
				#accuracy
				data_dict[type][0].append(float(row[4]))
				#no_of_predictions
				data_dict[type][1].append(float(row[6]))
				#position of no_of_genes depends on type
				if type=='known' or type=='unknown':
					data_dict[type][2].append(float(row[7]))
				else:
					data_dict[type][2].append(float(row[8]))
		sys.stderr.write('Done.\n')
		return fix_index_ls, var_index, variant_ls, data_dict
	
	def run(self):
		"""
		11-28-05
		11-30-05
			add training_accuracy and known_accuracy
		"""
		fix_index_ls, var_index, variant_ls, data_dict = self.get_data_structure(self.inputfname, self.parameter_list)
		self.plot(self.outputfname, fix_index_ls, self.parameter_list, var_index, variant_ls, self.parameter_index2label, \
			data_dict['testing'][0], 'testing_accuracy')
		self.plot(self.outputfname, fix_index_ls, self.parameter_list, var_index, variant_ls, self.parameter_index2label, \
			data_dict['training'][0], 'training_accuracy')	#11-30-05
		self.plot(self.outputfname, fix_index_ls, self.parameter_list, var_index, variant_ls, self.parameter_index2label, \
			data_dict['known'][2], 'known_no_genes')
		self.plot(self.outputfname, fix_index_ls, self.parameter_list, var_index, variant_ls, self.parameter_index2label, \
			data_dict['known'][0], 'known_accuracy')	#11-30-05
		self.plot(self.outputfname, fix_index_ls, self.parameter_list, var_index, variant_ls, self.parameter_index2label, \
			data_dict['unknown'][2], 'unknown_no_genes')
		

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hi:o:p:br", ["help"])
	except:
		print __doc__
		sys.exit(2)
	
	inputfname = None
	outputfname = None
	parameter_list = [0.01,-1,0.5]
	debug = 0
	report = 0
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-i"):
			inputfname = arg
		elif opt in ("-o"):
			outputfname = arg
		elif opt in ("-p"):
			parameter_list = arg
			parameter_list = parameter_list.split(',')
			if parameter_list[1]!='-1':	#11-28-05 parse the a/b into float
				numerator, denominator = parameter_list[1].split('/')
				parameter_list[1] = float(numerator)/float(denominator)
			parameter_list = map(float, parameter_list)
		elif opt in ("-b"):
			debug = 1
		elif opt in ("-r"):
			report = 1
	if inputfname and outputfname and len(parameter_list)==3:
		instance = DrawRpartValidationCurves(inputfname, outputfname, parameter_list,  \
			debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
