#!/usr/bin/env python
"""
Usage: triplet_plot.py  [OPTION] TRIPLET_STAT_FILE PLOT_FNAME

Option:
	-i ..., --init_recurrence=...	1(default)
	-x ..., --xlim=...	the range for the x axis, 0.0-20(default)
	-y ..., --ylim=...	the range for the y axis, 0.0-1.0(default)
	-t, --transcription	plot the homogeneity curve of transcription(function, default)
	-h, --help              show this help
	
Examples:
	triplet_plot.py  -i 2 2-14 function_homogeneity.png

Description:
	A program to draw homogeneity curves based on a matrix from TRIPLET_STAT_FILE.
	Columns of odd number are function_homogeneity ratios.
	Columns of even number are transcription_homogeneity ratios.
"""

import sys, os, getopt
from rpy import *

class triplet_plot:
	def __init__(self, init_recurrence, xlim, ylim, transcription, infname, ofname):
		self.infname = infname
		self.ofname = ofname
		self.init_recurrence = init_recurrence
		self.x_range = xlim.split('-')
		self.x_range = map(float, self.x_range)
		self.y_range = ylim.split('-')
		self.y_range = map(float, self.y_range)
		self.transcription = int(transcription)
		self.recurrence_list = []
		self.avg_list = []
		self.var_list = []
		
	def dstruc_loadin(self):
		data = with_mode(0, r.read_table)(self.infname)
		self.matrix = r.as_matrix(data)
		#column 0,2,4... is functional homogeneity
		#column 1,3,5... is transcriptional homogeneity
		
	def stat(self):
		for i in range(len(self.matrix[0,:])):
			if (i%2) == self.transcription:
				#construct three lists
				self.recurrence_list.append(self.init_recurrence+(i/2))
				self.avg_list.append(average(self.matrix[:,i]))
				self.var_list.append(r.sd(self.matrix[:,i]))
	
	def plot(self):
		self.conf_interval = r.qnorm(0.95)
		r.png('%s'%self.ofname)
		#first draw the main line
		r.plot(self.recurrence_list, self.avg_list, type='o',pch='*',xlab='recurrence',xlim=self.x_range,ylim=self.y_range, \
			ylab='homogeneity')
		for i in range(len(self.recurrence_list)):
			#draw a segment to show confidence interval
			y1 = self.avg_list[i] - (self.conf_interval*self.var_list[i]/sqrt(100))
			y2 = self.avg_list[i] + (self.conf_interval*self.var_list[i]/sqrt(100))
			r.segments(self.recurrence_list[i], y1, self.recurrence_list[i], y2)
		#don't forget to close the device.
		r.dev_off()

	def run(self):
		self.dstruc_loadin()
		self.stat()
		self.plot()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_option_list = ["help", "init_recurrence=", "xlim=", "ylim=", "transcription"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hi:x:y:t", long_option_list)
	except:
		print __doc__
		sys.exit(2)
	
	init_recurrence = 1
	xlim = '0.0-20'
	ylim = '0.0-1.0'
	transcription = 0
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-i", "--init_recurrence"):
			init_recurrence = int(arg)
		elif opt in ("-x", "--xlim"):
			xlim = arg
		elif opt in ("-y", "--ylim"):
			ylim = arg
		elif opt in ("-t", "--transcription"):
			transcription = 1
	
	if len(args)==2:
		instance = triplet_plot(init_recurrence, xlim, ylim, transcription, args[0], args[1])
		instance.run()
	else:
		print __doc__
		sys.exit(2)
