#!/usr/bin/env python
"""
functions or classes common to all programs

"""

import csv, sys
import psycopg

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
		new_list.append(dict.get(item))
	return new_list

def db_connect(hostname, dbname, schema):
	"""
	02-28-05
		establish database connection, return (conn, curs).
		copied from CrackSplat.py
	"""
	conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
	curs = conn.cursor()
	curs.execute("set search_path to %s"%schema)
	return (conn, curs)

def get_gene_no2direct_go(curs, schema=None, table='association'):
	"""
	03-02-05
		return gene_no2direct_go
	"""
	sys.stderr.write("Getting gene_no2direct_go...")
	gene_no2direct_go = {}
	if schema:
		curs.execute("set search_path to %s"%schema)
	curs.execute("select ge.gene_no, a.go_id from graph.association a, gene ge\
			where ge.gene_id=a.gene_id")
	rows = curs.fetchall()
	for row in rows:
		if row[0] in gene_no2direct_go:
			gene_no2direct_go[row[0]].append(row[1])
		else:
			gene_no2direct_go[row[0]] = [row[1]]
	sys.stderr.write("Done\n")
	return gene_no2direct_go
	

def get_gene_no2gene_id(curs, schema=None, table='gene'):
	"""
	03-02-05
		return gene_no2gene_id
	"""
	sys.stderr.write("Getting gene_no2gene_id...")
	gene_no2gene_id = {}
	if schema:
		curs.execute("set search_path to %s"%schema)
	curs.execute("select gene_no, gene_id from %s"%table)
	rows = curs.fetchall()
	for row in rows:
		gene_no2gene_id[row[0]] = row[1]
	sys.stderr.write("Done\n")
	return gene_no2gene_id

def get_go_no2go_id(curs, schema=None, table='go'):
	"""
	03-02-05
		return go_no2go_id
	"""
	sys.stderr.write("Getting go_no2go_id...")
	go_no2go_id = {}
	if schema:
		curs.execute("set search_path to %s"%schema)
	curs.execute("select go_no, go_id, depth from %s"%table)
	rows = curs.fetchall()
	for row in rows:
		go_no2go_id[row[0]] = row[1]
	sys.stderr.write("Done\n")
	return go_no2go_id

def get_go_no2name(curs, schema=None, table='go'):
	"""
	03-03-05
		return go_no2name
	"""
	sys.stderr.write("Getting go_no2name...")
	go_no2name = {}
	if schema:
		curs.execute("set search_path to %s"%schema)
	curs.execute("select go_no, name from %s"%table)
	rows = curs.fetchall()
	for row in rows:
		go_no2name[row[0]] = row[1]
	sys.stderr.write("Done\n")
	return go_no2name
