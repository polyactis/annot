#!/usr/bin/env python
"""
03-08-05:
	Something wrong after I modified copath's output filename. Need to filter some overheads out.
Overhead is like this:

***********************************************************
cluster: 34
***********************************************************
(1 35 1)(1 51 1)(1 124 1)(1 225 1)(6 215 1)(6 225 1)
(35 51 1)(124 215 1)(215 225 1)

Usage:
	copath_output_filter.py input_file output_file

"""

import sys, os, getopt

if len(sys.argv)!=3:
	print __doc__
	sys.exit(1)

inf = open(sys.argv[1], 'r')
outf = open(sys.argv[2], 'w')

for line in inf:
	first_character = line[0]
	if first_character=='c' or first_character == '(' or first_character == '*':
		continue
	else:
		outf.write(line)
