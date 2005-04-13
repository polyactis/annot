#!/usr/bin/env python
"""
04-12-05
Usage:
	test_biclustering.py INPUTFILE

Decription:
	A program to test the biclustering module.
	INPUTFILE is tab-delimited format. First column is the label.
	'NA' denotes Not avaible.
	
"""
import csv,sys
import biclustering

if len(sys.argv)!=2:
	print __doc__
	sys.exit(0)

instance = biclustering.biclustering(50,5,6,100)
filename=sys.argv[1]
inf = open(filename,'r')
reader = csv.reader(inf,delimiter='\t')
ls = []
for row in reader:
	if row[-1]:	#some datasets' last column is None
		ls.append(row[1:])
	else:
		ls.append(row[1:-1])
ls[0]
ls[432]
len(ls)
instance.data_read_in(ls)
ls1 = instance.return_matrix_data()
result = instance.getbicluster()
result = instance.getbicluster()
result = instance.getbicluster()
print result
