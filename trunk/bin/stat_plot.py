#!/usr/bin/env python
"""
Usage: stat_plot.py -f FIXED -v VARYING -t TAG [options] OFNAME

Options:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...		database to connect, graphdb(default)
	-f ..., --fixed=...	fixed parameter values, like urs:0.25,8,8
	-v ..., --var=...	varying parameter, one of unrs
	-t ..., --tag=...	just tag
	-x ..., --xlim=...	the range for the x axis, 0.0-1000(default)
	-y ..., --ylim=...	the range for the y axis, 0.0-1.0(default)
	-c, --based_on_clusters	default is based_on_functions
	-l, --l1	L1(sibling) is counted as true positive.
	-h, --help		show this help

Examples:
	stat_plot.py -t sc_known

Description:
	unrs: u denotes unknown_cut_off, n denotes connectivity_cut_off
		r denotes recurrence_cut_off, s denotes cluster_size_cut_off
		
"""

import sys, os, psycopg, getopt
from rpy import r

class option_attr:
	def __init__(self, label=None, value=None):
		self.label = label
		self.value = value

class stat_plot:
	def __init__(self, hostname, dbname, fixed, var, tag, xlim, ylim, based_on_clusters, l1, ofname):
		self.conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to graph")
		(self.fixed_label, self.fixed_value) = fixed.split(':')
		self.fixed_value_list = self.fixed_value.split(',')
		self.var = var
		self.tag = tag
		self.x_range = xlim.split('-')
		self.x_range = map(float, self.x_range)
		self.y_range = ylim.split('-')
		self.y_range = map(float, self.y_range)
		self.based_on_clusters = int(based_on_clusters)
		self.l1 = int(l1)
		self.ofname = ofname
		self.option_label_dict = {'u': 'unknown_cut_off',
			'n':'connectivity_cut_off',
			'r':'recurrence_cut_off',
			's':'cluster_size_cut_off'}
		self.option_num_dict = {}
		self.stat_data = []
		self.no_of_curves = 0
	
	def option_parsing(self):
		for i in range(len(self.fixed_label)):
			label = self.option_label_dict[self.fixed_label[i]]
			self.option_num_dict[i] = option_attr(label, self.fixed_value_list[i])
		self.option_num_dict[3] = option_attr(self.option_label_dict[self.var])

	def run(self):
		self.option_parsing()
		self.curs.execute("select distinct %s, %s, %s, %s, tag from\
			stat_plot_data where %s=%s and %s=%s and %s=%s and tag='%s' order by %s \
			"%(self.option_num_dict[0].label, self.option_num_dict[1].label, self.option_num_dict[2].label,\
			self.option_num_dict[3].label, self.option_num_dict[0].label, self.option_num_dict[0].value, \
			self.option_num_dict[1].label, self.option_num_dict[1].value, self.option_num_dict[2].label, \
			self.option_num_dict[2].value, self.tag, self.option_num_dict[3].label))
		rows = self.curs.fetchall()
		r.png('%s'%self.ofname)
		for row in rows:
			self.plot(row)
		r.dev_off()
		
	def plot(self, row):
		self.no_of_curves += 1
		x_list = []
		y_list = []
		self.curs.execute("select tp, tp_m, tp1, tp1_m, tn, fp, fp_m, fn from stat_plot_data where\
			%s=%s and %s=%s and %s=%s and %s=%s and tag='%s' order by p_value_cut_off"%\
			(self.option_num_dict[0].label, row[0], self.option_num_dict[1].label, row[1], self.option_num_dict[2].label,\
			row[2], self.option_num_dict[3].label, row[3], row[4]))
		plot_data = self.curs.fetchall()
		for entry in plot_data:
			tn = entry[4]
			fn = entry[7]
			if self.based_on_clusters:
				#using the tp_m, tp1_m and fp_m
				tp = entry[1]
				tp1 = entry[3]
				fp = entry[6]
			else:
				#using the tp, tp1, fp
				tp = entry[0]
				tp1 = entry[2]
				fp = entry[5]
			if self.l1:
				#tp1 is counted as true positive
				tp += tp1
			else:
				#tp1 is counted as false positive
				fp += tp1
			x_list.append(tp)
			y_list.append(float(tp)/(tp+fp))

		if self.no_of_curves==1:
			r.plot(x_list, y_list, type='o',pch='*',xlab='true positive',xlim=self.x_range,ylim=self.y_range, \
			ylab='ratio', main='%s'%(self.option_num_dict[3].label), col=self.no_of_curves)
		else:
			r.lines(x_list, y_list, type='o',pch='*',col=self.no_of_curves)


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_option_list = ["help", "hostname=", "dbname=", "fixed=", "var=", "tag=", "xlim=", "ylim=", "based_on_clusters", "l1"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:f:v:t:x:y:cl", long_option_list)
	except:
		print __doc__
		sys.exit(2)
		
	hostname = 'zhoudb'
	dbname = 'graphdb'
	fixed = None
	var = None
	tag = None
	xlim = '0.0-1000'
	ylim = '0.0-1.0'
	based_on_clusters = 0
	l1 = 0
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-z", "--hostname"):
			hostname = arg
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-f", "--fixed"):
			fixed = arg
		elif opt in ("-v", "--var"):
			var = arg
		elif opt in ("-t", "--tag"):
			tag = arg
		elif opt in ("-x", "--xlim"):
			xlim = arg
		elif opt in ("-y", "--ylim"):
			ylim = arg
		elif opt in ("-c", "--based_on_clusters"):
			based_on_clusters = 1
		elif opt in ("-l", "--l1"):
			l1 = 1
	
	if fixed and var and tag and len(args)==1:
		instance = stat_plot(hostname, dbname, fixed, var, tag, xlim, ylim, based_on_clusters, l1, args[0])
		instance.run()
	else:
		print __doc__
		sys.exit(2)
