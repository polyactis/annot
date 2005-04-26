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

"""
04-26-05
	procedures to check output examples of MpiBiclustering.py
#########check the biclustering score (Cheng2000), first through python

from graph import biclustering
from codense.common import get_edge_vector_by_id, db_connect
conn,curs = db_connect("zhoudb","graphdb","sc_54_6661")
import Numeric
import MLab
edge_id_list = [3450, 4360, 4412, 4602, 4691]
col_index_list =   [0, 4, 30, 34, 35, 38]

cor_2d_list, sig_2d_list = get_edge_vector_by_id(curs, edge_id_list)
cor_2d_array = Numeric.array(cor_2d_list)



seed_array = Numeric.take(cor_2d_array, col_index_list,1)
colMean = MLab.mean(seed_array)
rowMean = MLab.mean(seed_array,1)
mean = sum(sum(seed_array))/(seed_array.shape[0]*seed_array.shape[1])
score_array = MLab.transpose(seed_array-colMean)-rowMean +mean

Hscore = 0
for i in range(seed_array.shape[0]):
	for j in range(seed_array.shape[1]):
		r = seed_array[i,j]-rowMean[i]-colMean[j]+mean
		r = r*r
		Hscore += r

print Hscore/(seed_array.shape[0]*seed_array.shape[1])

##########check the extended candidate edges

cand_id_list =  [757, 760, 1110, 1651, 3278, 3641, 3730, 4310]
cor_2d_list, sig_2d_list = get_edge_vector_by_id(curs, edge_id_list)
cor_2d_array = Numeric.array(cor_2d_list)
cand_array = Numeric.take(cor_2d_array, col_index_list,1)

from rpy import r

r.dist(r.rbind(can_array[0,:],colMean))
r.cor(can_array[0,:],colMean)

##########check the c++ module

instance = biclustering.biclustering(0.01,5,6,100)

biclustering.set_module_and_type('Numeric', 'ArrayType')
instance.data_read_in(seed_array)
ls1 = instance.return_matrix_data()
result = instance.getbicluster()

"""
