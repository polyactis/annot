#!/usr/bin/env python
'''
Usage:	p_value_cor.py [OPTIONS] OUTPUT_TABLE_FILE

Option:
	OUTPUT_TABLE_FILE is a file to store the p_value correlation table.
	-l ... --df_lower=...,	Lower degree_of_freedom, 5(default)
	-u ... --df_upper=...,	Upper degree_of_freedom, 100(default)
	-c ... --column=...,	specify which column to plot via R, 4(default)

Examples:
	p_value_cor.py p_value_cor.output
	p_value_cor.py -l 5 -u 200 -c 5 p_value_cor.output

Description:
	We assume that correlation value between two gene expression 
	vectors is distributed as t-distribution. The degree_of_freedom
	is the number of measurements.

	Program is to generate correlation values under different p_values
	and different degree_of_freedom. 

	Output will dumped to stdout. File, 'p_value_cor.pdf' is generated in
	the current directory, which is a plot, correlation v.s df.

'''

import sys, getopt
from rpy import r
from numarray import *

class p_value_cor:
	def __init__(self, df_lower, df_upper, column, output):
		self.p_value_list = [0.25, 0.1, 0.05, 0.025, 0.01, 0.005, 0.001, 0.0005, 0.0001]
		self.df_lower = int(df_lower)
		self.df_upper = int(df_upper)
		#index is 1 less than column
		self.column = int(column-1)
		self.result_array = []
		self.of = open(output, 'w')
	
	def run(self):
		for df in range(self.df_lower, self.df_upper+1):
			cor_list = []
			for p_value in self.p_value_list:
				t=r.qt(p_value,df,lower_tail=r.FALSE)
				cor = r.sqrt(t*t/(t*t+df))
				cor_list.append(cor)
			self.result_array.append(cor_list)
		#output the table and plot a sample
		self.output()
		
	def output(self):
		self.matrix = array(self.result_array)
		p_value_list = map(str, self.p_value_list)
		self.of.write('p_value\t%s\n'%'\t'.join(p_value_list))
		df = self.df_lower
		for cor_list in self.result_array:
			cor_list = map(str, cor_list)	#string can be 'join'ed. easy to output
			self.of.write('%d\t%s\n'%(df, '\t'.join(cor_list)))
			df += 1
		r.pdf('p_value_cor.pdf')
		#select a column to plot
		cor_list = self.matrix[:,self.column]
		p_value_label = self.p_value_list[self.column]
		df_list = range(self.df_lower, self.df_upper+1)
		r.plot(df_list, cor_list, type='o', pch='*', xlab='df', ylab='correlation', main='p_value: %s'%p_value_label)
		r.dev_off()


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hl:u:c:", ["help", "df_lower=", "df_upper=", "column="])
	except:
		print __doc__
		sys.exit(2)
	
	df_lower = 5
	df_upper = 100
	column = 4
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-l", "--df_lower"):
			df_lower = int(arg)
		elif opt in ("-u", "--df_upper"):
			df_upper = int(arg)
		elif opt in ("-c", "--column"):
			column = int(arg)

	if len(args) == 1:
		instance = p_value_cor(df_lower, df_upper, column, args[0])
		instance.run()
	else:
		print __doc__
		sys.exit(2)
