#!/usr/bin/env python
"""
Usage: stat_plot.py -f FIXED -v VARYING -t TAG [options] OFNAME

Options:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...		database to connect, graphdb(default)
	-f ..., --fixed=...	fixed parameter values, like ursp:0.25,8,8,0.001
	-v ..., --var=...	varying parameter, one of unrsp
	-t ..., --tag=...	just tag
	-x ..., --xlim=...	the range for the x axis, 0.0-1000(default)
	-y ..., --ylim=...	the range for the y axis, 0.0-1.0(default)
	-p ..., --type=...	the type of running, 0(defaut), 1, 2
	-c, --based_on_clusters	default is based_on_functions
	-l, --l1	L1(sibling) is counted as true positive.
	-h, --help		show this help

Examples:
	stat_plot.py -t sc_known -f ursp:0.25,8,8,0.001 -v n connectivity.jpg

Description:
	unrsp: u denotes unknown_cut_off, n denotes connectivity_cut_off
		r denotes recurrence_cut_off, s denotes cluster_size_cut_off
		p denotes p_value_cut_off
	type of running:
		0:	two varying parameters, plotting
		1:	one varying parameter, plotting
		2:	table_output, preparation for R
"""

import sys, os, psycopg, getopt, csv
from rpy import r

class option_attr:
	def __init__(self, label=None, value=None):
		self.label = label
		self.value = value

class stat_plot:
	def __init__(self, hostname, dbname, fixed, var, tag, xlim, ylim, based_on_clusters, type, l1, ofname):
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
		self.type = int(type)
		self.l1 = int(l1)
		self.ofname = ofname
		self.option_label_dict = {'u': 'unknown_cut_off',
			'n':'connectivity_cut_off',
			'r':'recurrence_cut_off',
			's':'cluster_size_cut_off',
			'p':'p_value_cut_off'}
		self.option_num_dict = {}
		self.stat_data = []
		self.no_of_curves = 0
		# a list of varying values, for legend
		self.varying_list = []
	
	def option_parsing(self):
		#in option_num_dict, 0,1,2,3 are the fixed parameters. 4 is the varying parameter.
		for i in range(len(self.fixed_label)):
			label = self.option_label_dict[self.fixed_label[i]]
			self.option_num_dict[i] = option_attr(label, self.fixed_value_list[i])
		self.option_num_dict[4] = option_attr(self.option_label_dict[self.var])

	def run(self):
		if self.type == 0:
			self.option_parsing()
			self.plot()
		elif self.type == 1:
			self.option_parsing()
			self.single_plot()
		elif self.type == 2:
			self.table_output()
	
	def table_output(self):
		'''
		Output a three-column table for 3D plotting via R's wireframe from lattice.
		The three columns are recurrence_cut_off, connectivity_cut_off, accuracy.
		'''
		writer = csv.writer(open(self.ofname,'w'), delimiter='\t')
		self.curs.execute("select  tp, tp_m, tp1, tp1_m, tn, fp, fp_m, fn, recurrence_cut_off, \
			connectivity_cut_off from stat_plot_data where tag='%s' order by recurrence_cut_off,\
			connectivity_cut_off"%(self.tag))
		plot_data = self.curs.fetchall()
		for entry in plot_data:
			tn = entry[4]
			fn = entry[7]
			recurrence_cut_off = entry[8]
			connectivity_cut_off = entry[9]
			
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
			if (tp+fp) != 0:
				accuracy = float(tp)/(tp+fp)
			else:
				#skip those unavailable data
				continue
			writer.writerow([recurrence_cut_off, connectivity_cut_off, accuracy])
		del writer		
		
	def single_plot(self):
		#this function deals with 4 fixed parameters and 1 varying parameter
		r.png('%s'%self.ofname)
		x_list = []
		y_list = []
		self.curs.execute("select  tp, tp_m, tp1, tp1_m, tn, fp, fp_m, fn from\
			stat_plot_data where %s=%s and %s=%s and %s=%s and %s=%s and tag='%s' order by %s \
			"%(self.option_num_dict[0].label, self.option_num_dict[0].value, \
			self.option_num_dict[1].label, self.option_num_dict[1].value, self.option_num_dict[2].label, \
			self.option_num_dict[2].value, self.option_num_dict[3].label, self.option_num_dict[3].value, \
			self.tag, self.option_num_dict[4].label))
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
		
		r.plot(x_list, y_list, type='o',pch='*',xlab='consistent predictions',xlim=self.x_range,ylim=self.y_range, \
				ylab='percentage', main='%s'%(self.option_num_dict[4].label), col=1)
		r.dev_off()
		
	def plot(self):
		#this function deals with 3 fixed parameters and 1 varying parameter
		self.curs.execute("select distinct %s, %s, %s, %s, tag from\
			stat_plot_data where %s=%s and %s=%s and %s=%s and tag='%s' order by %s \
			"%(self.option_num_dict[0].label, self.option_num_dict[1].label, self.option_num_dict[2].label,\
			self.option_num_dict[3].label, self.option_num_dict[0].label, self.option_num_dict[0].value, \
			self.option_num_dict[1].label, self.option_num_dict[1].value, self.option_num_dict[2].label, \
			self.option_num_dict[2].value, self.tag, self.option_num_dict[3].label))
		rows = self.curs.fetchall()
		r.png('%s'%self.ofname)
		for row in rows:
			#position 0,1,2 are fixed values, 3 is varying value, 4 is the tag value.
			self._plot(row)
		#add the legend
		r.legend(self.x_range[1], self.y_range[1], self.varying_list, col=range(1, self.no_of_curves+1), lty=1, pch='*', xjust=1)
		r.dev_off()
		
	def _plot(self, row):
		self.no_of_curves += 1
		#position 3 in row is the varying value
		self.varying_list.append(row[3])
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
			r.plot(x_list, y_list, type='o',pch='*',xlab='consistent predictions',xlim=self.x_range,ylim=self.y_range, \
			ylab='percentage', main='%s'%(self.option_num_dict[3].label), col=self.no_of_curves)
		else:
			r.lines(x_list, y_list, type='o',pch='*',col=self.no_of_curves)


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_option_list = ["help", "hostname=", "dbname=", "fixed=", "var=", "tag=", "xlim=", "ylim=", "based_on_clusters", "type=", "l1"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:f:v:t:x:y:p:cl", long_option_list)
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
	type = 0
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
		elif opt in ("-p", "--type"):
			type = int(arg)
		elif opt in ("-l", "--l1"):
			l1 = 1
	
	if fixed and var and tag and len(args)==1:
		instance = stat_plot(hostname, dbname, fixed, var, tag, xlim, ylim, based_on_clusters, type, l1, args[0])
		instance.run()
	else:
		print __doc__
		sys.exit(2)
