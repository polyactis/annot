#!/usr/bin/env python
"""
functions or classes common to all programs

"""

import csv

def index_plus_one(i):
	'''
	This small function is used in map() to increase the indices by 1.
	'''
	return int(i)+1
	
def divided_by_1000(i):
	return '%1.3f'%(float(i)/float(1000))


def get_haiyan_no2gene_no(table_file):
	reader = csv.reader(file(table_file), delimiter='\t')
	haiyan_no2gene_no={}
	no = 0
	for row in reader:
		gene_no = int(row[1])
		haiyan_no2gene_no[no] = gene_no
		no += 1
	del reader
	return haiyan_no2gene_no


def get_haiyan_no2gene_id(table_file):
	reader = csv.reader(file(table_file), delimiter='\t')
	haiyan_no2gene_id={}
	no = 0
	for row in reader:
		gene_id = row[0]
		haiyan_no2gene_id[no] = gene_id
		no += 1
	del reader
	return haiyan_no2gene_id


def dict_map(dict, ls):
	new_list = []
	for item in ls:
		new_list.append(dict[item])
	return new_list
