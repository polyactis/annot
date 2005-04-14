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
instance.data_read_in(ls)
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
	y = raw_input("Continue:(Y/n)?")
	if y=="n":
		break
	result = instance.getbicluster()
