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
	sys.stderr.write("Getting haiyan_no2gene_no...")
	reader = csv.reader(file(table_file), delimiter='\t')
	haiyan_no2gene_no={}
	no = 0
	for row in reader:
		gene_no = int(row[1])
		haiyan_no2gene_no[no] = gene_no
		no += 1
	del reader
	sys.stderr.write("Done\n")
	return haiyan_no2gene_no


def get_haiyan_no2gene_id(table_file):
	sys.stderr.write("Getting haiyan_no2gene_id...")
	reader = csv.reader(file(table_file), delimiter='\t')
	haiyan_no2gene_id={}
	no = 0
	for row in reader:
		gene_id = row[0]
		haiyan_no2gene_id[no] = gene_id
		no += 1
	del reader
	sys.stderr.write("Done\n")
	return haiyan_no2gene_id


def dict_map(dict, ls):
	new_list = []
	for item in ls:
		value = dict.get(item)
		if value:
			new_list.append(value)
	return new_list

def db_connect(hostname, dbname, schema=None):
	"""
	02-28-05
		establish database connection, return (conn, curs).
		copied from CrackSplat.py
	03-08-05
		parameter schema is optional
	"""
	conn = psycopg.connect('host=%s dbname=%s'%(hostname, dbname))
	curs = conn.cursor()
	if schema:
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


def get_gene_id2gene_no(curs, schema=None, table='gene'):
	"""
	03-02-05
		return gene_id2gene_no
	"""
	sys.stderr.write("Getting gene_id2gene_no...")
	gene_id2gene_no = {}
	if schema:
		curs.execute("set search_path to %s"%schema)
	curs.execute("select gene_no, gene_id from %s"%table)
	rows = curs.fetchall()
	for row in rows:
		gene_id2gene_no[row[1]] = row[0]
	sys.stderr.write("Done\n")
	return gene_id2gene_no


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

def get_go_no2term_id(curs, schema=None, term_table='go.term'):
	"""
	03-04-05
		get the go_no2term_id dictionary
	"""
	sys.stderr.write("Getting go_no2term_id...")
	if schema:
		curs.execute("set search_path to %s"%schema)
	go_no2term_id = {}
	curs.execute("select g.go_no, t.name, t.id from go g, %s t where g.go_id=t.acc"%(term_table))
	rows = curs.fetchall()
	for row in rows:
		go_no2term_id[row[0]] = row[2]
	sys.stderr.write("Done\n")
	return go_no2term_id

def get_gene_no2go_no(curs, schema=None, gene_table='gene'):
	"""
	03-09-05
		get the gene_no2go_no
	"""
	sys.stderr.write("Getting gene_no2go_no...")
	if schema:
		curs.execute("set search_path to %s"%schema)
	
	gene_no2go_no = {}
	
	curs.execute("select gene_no,go_functions from gene")
	rows = curs.fetchall()
	for row in rows:
		gene_no2go_no[row[0]] = []
		go_functions_list = row[1][1:-1].split(',')
		for go_no in go_functions_list:
			gene_no2go_no[row[0]].append(int(go_no))
	sys.stderr.write("Done\n")
	return gene_no2go_no
