#!/usr/bin/env python
'''
We assume that correlation value between two gene expression vectors is
distributed as t-distribution. The degree_of_freedom is the number of measurements.

Program is to generate correlation values under different p_values and different
degree_of_freedom. 

Usage: p_value_cor.py

Output will dumped to stdout. File, 'p_value_cor.png' is generated in the current
directory, which is a plot, correlation v.s df.

'''

import sys
from rpy import r
from numarray import *

class p_value_cor:
	def __init__(self):
		self.p_value_list = [0.25, 0.1, 0.05, 0.025, 0.01, 0.005, 0.001, 0.0005, 0.0001]
		self.df_lower = 5
		self.df_upper = 100
		self.result_array = []
		
	def run(self):
		for df in range(self.df_lower, self.df_upper+1):
			cor_list = []
			for p_value in self.p_value_list:
				t=r.qt(p_value,df,lower_tail=r.FALSE)
				cor = r.sqrt(t*t/(t*t+df))
				cor_list.append(cor)
			self.result_array.append(cor_list)
	
	def output(self, of=sys.stdout):
		self.matrix = array(self.result_array)
		p_value_list = map(str, self.p_value_list)
		of.write('p_value\t%s\n'%'\t'.join(p_value_list))
		df = self.df_lower
		for cor_list in self.result_array:
			cor_list = map(str, cor_list)
			of.write('%d\t%s\n'%(df, '\t'.join(cor_list)))
			df += 1
		r.pdf('p_value_cor.pdf')
		cor_list = self.matrix[:,3]
		p_value_label = self.p_value_list[3]
		df_list = range(self.df_lower, self.df_upper+1)
		r.plot(df_list, cor_list, type='o', pch='*', xlab='df', ylab='correlation', main='p_value: %s'%p_value_label)
		r.dev_off()
		
if __name__ == '__main__':
	instance = p_value_cor()
	instance.run()
	instance.output()
