#!/usr/bin/env python
"""
Usage: arguments2string.py ARGUMENTS

	
Examples:
	arguments2string.py -s 100 -m -t 0.58 -u 0,2 -i da-s
		--prints s100mt0_58u0_2ida_s
	
Description:
	Return a string form of the arguments.
	Warning: Front - can only be used in option,
	not argument. Otherwise, that argument
	is regarded as argument. i.e -c -1
	-1 is regarded as option. Prints c1.

"""

import sys, os

def arguments2string(arguments_list):
	"""
	10-05-05
	"""
	string = ''
	for argument in arguments_list:
		if argument[0] == '-':	#this is option
			if argument[1] == '-':	#long option
				string += argument[2:]
			else:	#short option
				string += argument[1:]
		else:	#this is argument
			argument = argument.replace('-','_')	#replace '-' with '_'
			argument = argument.replace('.','_')	#replace '.' with '_'
			string += argument.replace(',','_')	#replace ',' with '_'
			
	return string

if __name__ == '__main__':
	if len(sys.argv)==1:
		print __doc__
		sys.exit(0)
	else:
		print arguments2string(sys.argv[1:])
