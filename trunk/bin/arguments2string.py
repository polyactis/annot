#!/usr/bin/env python
"""
Usage: arguments2string.py ARGUMENTS

	
Examples:
	arguments2string.py -s 100 -m -t 0.58 -u 0.2 -i das
		--prints s100mt58u2idas
	
Description:
	Return a string form of the arguments.

"""

import sys, os

def arguments2string(arguments_list):
	"""
	10-05-05
	"""
	string = ''
	for argument in arguments_list:
		if argument[0] == '-':
			if argument[1] == '-':	#long option
				string += argument[2:]
			else:
				string += argument[1:]
		else:
			if argument[:2] == '0.':
					argument = argument[2:]
			string += argument
	return string

if __name__ == '__main__':
	print arguments2string(sys.argv[1:])
