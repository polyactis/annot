#!/usr/bin/env python
"""
argv[1] is the input file containing stats from graph.cc
argv[2] is output file to store the percentage data
three plots will be generated as pdf files(started with per_xxx).
"""

import sys, os, csv, getopt
from rpy import r

class edge_percentage_plot:
	def __init__(self, infname, outfname):
		self.csv_reader = csv.reader(file(infname), delimiter='\t')
		self.csv_writer = csv.writer(open(outfname,'w'), delimiter = '\t')
		self.dataset_no = []
		self.per_05 = []
		self.per_025 = []
		self.per_01 = []
		
	def run(self):
		for ls in self.csv_reader:
			ls[1:] = map(float, ls[1:])
			self.dataset_no.append(int(ls[0]))
			self.per_05.append(ls[2]/ls[1])
			self.per_025.append(ls[3]/ls[1])
			self.per_01.append(ls[4]/ls[1])
		for i in range(len(self.dataset_no)):
			ls = [self.dataset_no[i], self.per_05[i], self.per_025[i], self.per_01[i]]
			self.csv_writer.writerow(ls)
		self.plot()
		
	def plot(self):
		r.pdf("per_p_value05.pdf")
		r.plot(self.dataset_no, self.per_05, type='o', pch='*', xlab='dataset no.',\
			ylab='percentage', main='p_value: 0.05. #edges compared with correlation cut_off 0.6')
		r.dev_off()
		r.pdf("per_p_value025.pdf")
		r.plot(self.dataset_no, self.per_025, type='o', pch='*', xlab='dataset no.',\
			ylab='percentage', main='p_value: 0.025. #edges compared with correlation cut_off 0.6')
		r.dev_off()
		r.pdf("per_p_value01.pdf")
		r.plot(self.dataset_no, self.per_01, type='o', pch='*', xlab='dataset no.',\
			ylab='percentage', main='p_value: 0.01. #edges compared with correlation cut_off 0.6')
		r.dev_off()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	instance = edge_percentage_plot(sys.argv[1], sys.argv[2])
	instance.run()
