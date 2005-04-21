#!/usr/bin/env python
"""
04-12-05
Usage:
	test_biclustering.py INPUTFILE HSCORE_CUT_OFF

Decription:
	A program to test the biclustering module.
	INPUTFILE is tab-delimited format. First column is the label.
	'NA' denotes Not avaible.
	
"""
import csv,sys
import biclustering
from Numeric import array
import random

if len(sys.argv)!=3:
	print __doc__
	sys.exit(0)

filename=sys.argv[1]
hscore_cut_off = float(sys.argv[2])
inf = open(filename,'r')
reader = csv.reader(inf,delimiter='\t')
ls = []
for row in reader:
	one_row = []
	for data in row[1:]:	#first column is gene name
		if data == 'NA':
			data=random.randint(-800,800)
		elif data=='':	#some datasets' last column is None
			continue
		else:
			data = float(data)
		one_row.append(data)
	ls.append(one_row)

instance = biclustering.biclustering(hscore_cut_off,5,6,100)
biclustering.set_module_and_type('Numeric', 'ArrayType')
instance.data_read_in(array(ls))
ls1 = instance.return_matrix_data()

#output the matrix read in by the class
writer = csv.writer(open('/tmp/test_biclustering.out', 'w'), delimiter='\t')
for row in ls1:
	writer.writerow(row)
del writer

print "the first row is "
print ls1[0]
print "the penultimate row is "
print ls1[-2]
print "the final row is "
print ls1[-1]
print "total number of rows is "
print len(ls)
result = instance.getbicluster()
while result:
	print "Score is %s"%result.score
	print "%s rows and their indices are %s"%(len(result.row_index_list), repr(result.row_index_list))
	print "%s columns and their indices are %s"%(len(result.column_index_list), repr(result.column_index_list))
	print "consensus is %s"%repr(result.consensus_list)
	y = raw_input("Continue:(Y/n)?")
	if y=="n":
		break
	result = instance.getbicluster()
