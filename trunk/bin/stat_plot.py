#!/usr/bin/env python
"""
Usage: stat_plot.py -c CONNECTIVITY -l LIMIT -t TAG [options]

Options:
  -d ..., --dbname=...		database to connect, graphdb(default)
  -c ..., --connectivity=...	connectivity_cut_off
  -l ..., --limit=...	max no. of clusters related to one gene
  -t ..., --tag=...	just tag
  -x ..., --xlim=...	the range for the x axis, 0.0-1.0(default)
  -y ..., --ylim=...	the range for the y axis, 0.0-1.0(default)
  -o ..., --output=...	output file name
  -h, --help		show this help

Examples:
  stat_plot.py -c 0.8 -l 10 -t sc_known
  stat_plot.py -c 0.8-0.9 -l 10 -t sc_known
  stat_plot.py -c 0.8-0.9 -l 5-10 -t sc_known
  stat_plot.py -x 0.2-0.8 -y 0.3-0.6 -l 0 -t sc_shu_known
"""

import sys, os, cStringIO, psycopg, getopt
from rpy import r

class stat_plot:
	def __init__(self, dbname, c_range, l_range, tag, xlim, ylim, ofname):
		self.conn = psycopg.connect('dbname=%s'%dbname)
		self.curs = self.conn.cursor()
		self.curs.execute("set search_path to graph")
		self.x_range = xlim.split('-')
		self.x_range = map(float, self.x_range)
		self.y_range = ylim.split('-')
		self.y_range = map(float, self.y_range)
		self.c_range = c_range.split('-')
		if len(self.c_range) == 1:
			self.c_range.append(self.c_range[0])
		self.l_range = l_range.split('-')
		if len(self.l_range) == 1:
			self.l_range.append(self.l_range[0])
		self.tag = tag
		self.ofname = ofname
		self.stat_data = []
		self.no_of_curves = 0
		
	def run(self):
		self.curs.execute("select distinct connectivity, limit_of_cluster, tag from\
			stat_plot_data where connectivity >= %f and\
			connectivity <= %f and limit_of_cluster >= %d and limit_of_cluster <= %d and\
			tag='%s' order by connectivity, limit_of_cluster\
			"%(float(self.c_range[0]), float(self.c_range[1]),\
			int(self.l_range[0]), int(self.l_range[1]), self.tag))
		rows = self.curs.fetchall()
		r.pdf('%s.pdf'%self.ofname)
		for connectivity,limit,tag in rows:
			self.plot(connectivity,limit,tag)
		r.dev_off()
		
	def plot(self, connectivity, limit, tag):
		self.no_of_curves += 1
		sensitivity_list = []
		positive_predictive_value_list = []
		self.curs.execute("select tp,tn,fp,fn from stat_plot_data where\
			connectivity=%f and limit_of_cluster=%d and tag='%s' order by p_value_cut_off"%\
			(connectivity, limit, tag))
		plot_data = self.curs.fetchall()
		for entry in plot_data:
			sensitivity_list.append(float(entry[0])/(entry[0]+entry[3]))
			positive_predictive_value_list.append(float(entry[0])/(entry[0]+entry[2]))

		if self.no_of_curves==1:
			r.plot(sensitivity_list, positive_predictive_value_list, type='o',pch='*',xlab='sensitivity',xlim=self.x_range,ylim=self.y_range, \
			ylab='positive_predictive_value', main='Connectivity: %s  Limit: %s'%(repr(self.c_range),repr(self.l_range)), col=self.no_of_curves)
		else:
			r.lines(sensitivity_list, positive_predictive_value_list, type='o',pch='*',col=self.no_of_curves)


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_option_list = ["help", "dbname=", "connectivity=", "limit=", "tag=", "xlim=", "ylim=", "output="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hd:c:l:t:x:y:o:", long_option_list)
	except:
		print __doc__
		sys.exit(2)
	
	dbname = 'graphdb'
	connectivity = None
	limit = None
	tag = None
	xlim = '0.0-1.0'
	ylim = '0.0-1.0'
	output =''
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-c", "--connectivity"):
			connectivity = arg
		elif opt in ("-l", "--limit"):
			limit = arg
		elif opt in ("-t", "--tag"):
			tag = arg
		elif opt in ("-x", "--xlim"):
			xlim = arg
		elif opt in ("-y", "--ylim"):
			ylim = arg
		elif opt in ("-o", "--output"):
			output = arg
	
	if connectivity and limit and tag:
		if output=='':
			output = '%s_%s_%s'%(tag,connectivity,limit)
		instance = stat_plot(dbname, connectivity, limit, tag, xlim, ylim, output)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
