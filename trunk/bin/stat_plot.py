#!/usr/bin/env python
"""
Usage: python stat_plot.py [options]

Options:
  -d ..., --dbname=...   use specified database name
  -h, --help              show this help

Examples:

"""

import sys, os, cStringIO, psycopg, getopt, pickle
from gene_stat import gene_stat
from rpy import r

class stat_plot:
	def __init__(self, dbname):
		self.dbname = dbname
		self.connectivity_list = [0.8,0.9]
		self.limit_list = [5,10,15,20]
		self.p_value_cut_off_list = [0.1, 0.01, 0.001, 0.0005, 0.0001, 0.00001]
		self.stat_data = []
		self.no_of_plots = 0
		
	def run(self):
		for connectivity in self.connectivity_list:
			for limit in self.limit_list:
				plot_data = []
				for p_value_cut_off in self.p_value_cut_off_list:
					instance = gene_stat(self.dbname, p_value_cut_off, limit, connectivity)
					instance.dstruc_loadin()
					instance.run()
					self.stat_data.append([connectivity, limit, p_value_cut_off, \
					 instance.tp, instance.tn, instance.fp, instance.fn])
					plot_data.append([instance.tp, instance.tn, instance.fp, instance.fn])
					del instance
				self.plot(connectivity, limit, plot_data)
				
	def plot(self, connectivity, limit, plot_data):
		r.png('%d.png'%self.no_of_plots)
		sensitivity_list = []
		false_positive_ratio_list = []
		for entry in plot_data:
			sensitivity_list.append(entry[0]/(entry[0]+entry[3]))
			false_positive_ratio_list.append(entry[2]/(entry[1]+entry[2]))
		
		r.plot(sensitivity_list, false_positive_ratio_list, type='o',pch='*',xlab='sensitivity', \
			ylab='false_positive_ratio', main='Connectivity: %f  Limit: %d'%(connectivity,limit))
		r.dev_off()
		self.no_of_plots += 1
	
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hd:", ["help", "dbname="])
	except:
		print __doc__
		sys.exit(2)
		
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-d", "--dbname"):
			instance = stat_plot(arg)
			instance.run()
			pickle_fname = os.path.join(os.path.expanduser('~'),'pickle/stat_data')
			pickle.dump(instance.stat_data,open(pickle_fname,'w'))
