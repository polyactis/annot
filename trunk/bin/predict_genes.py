#!/usr/bin/env python
"""
Usage: predict_genes.py [OPTIONS] GENE_FILE TABLE_FILE

Option:
	GENE_FILE is the file containing all the yeast gene orfnames.
	TABLE_FILE is the output file.
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database, sc_yh60_splat(default)
	-t ..., --table=...	cluster_stat_sup_3_2(default)
	-m ..., --mcl_table=...	mcl_result_sup_3_2(default), corresponding to above table.
	-p ..., --p_value_cut_off=...	p_value_cut_off, 0.001 (default)
	-y ..., --recurrence_cut_off=...	6 (default), minimum recurrences
	-n ..., --connectivity_cut_off=...	0.8 (default), minimum connectivity of a mcl cluster
	-e ..., --email=...	email address to receive the result
	-r, --report	report the progress(a number)
	-h, --help              show this help

Examples:
	predict_genes.py genes.txt tab.csv

Description:
	A wrapper of gene_stat from gene_stat_plot.py for predict_genes.cgi.
	The program will output prediction table for the genes in the GENE_FILE.

"""

import sys, os, getopt, csv, smtplib
from gene_stat_plot import gene_stat
from sets import Set
from email.MIMEText import MIMEText

class predict_genes:
	def __init__(self, gene_infname='tmp/in', tab_ofname='tmp/out', hostname='zhoudb',\
		dbname='graphdb', schema='sc_yh60_splat', table='cluster_stat_sup_3_2', \
		mcl_table='mcl_result_sup_3_2', p_value_cut_off=0.001, \
		connectivity_cut_off=0.8, recurrence_cut_off=6, email='', report=0):
		self.gene_infname = gene_infname
		self.tab_ofname = tab_ofname
		self.tab_of = csv.writer(open(self.tab_ofname, 'w'), delimiter='\t')
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.table = table
		self.mcl_table = mcl_table
		self.p_value_cut_off = float(p_value_cut_off)
		self.recurrence_cut_off = int(recurrence_cut_off)
		self.connectivity_cut_off = float(connectivity_cut_off)
		self.email = email
		self.report = int(report)
		
		self.gene_set = Set()


	def dstruc_loadin(self):
		sys.stderr.write("Loading Data STructure...")
		reader = csv.reader(file(self.gene_infname), delimiter='\t')
		for row in reader:
			#the gene names must be capitalized.
			self.gene_set.add(row[0].upper())
		sys.stderr.write("Done\n")
	
	def run(self):
		self.gs_instance = gene_stat(self.hostname, self.dbname, self.schema, self.table, self.mcl_table, \
			self.p_value_cut_off, 1.0, self.connectivity_cut_off, self.recurrence_cut_off, 1000,\
			1, 1, self.report)
		self.gs_instance.dstruc_loadin()
		self.gs_instance.core()
		self.gs_instance.final()
		self.gs_instance.stat_output()
		self.gs_instance.go_no_accuracy()
		self.output()
		if self.email != '':
			self.send_email()

	def output(self):
		header = ['gene_id', 'function_known', 'function_predicted', 'is_correct', 'average p_value', \
		'expected accuracy', '#supporting clusters', 'cluster_context']
		self.tab_of.writerow(header)
		for gene_no in self.gs_instance.gene_prediction_dict:
			gene_id = self.gs_instance.gene_no2gene_id[gene_no]
			if gene_id in self.gene_set:
				#remove the predicted gene from gene_set
				self.gene_set.remove(gene_id)
				unit = self.gs_instance.gene_prediction_dict[gene_no].p_functions_struc_dict
				self.gs_instance._table_output(self.tab_of, gene_no, unit)
		
		for gene_id in self.gene_set:
			#the remaining unpredicted genes
			self.tab_of.writerow([gene_id, 'No Prediction'])
			
		#delete the csv writer and close the file handler.
		del self.tab_of
		
	def send_email(self):
		smtp_server = smtplib.SMTP('localhost')
		fromaddr = 'yuhuang@usc.edu'
		toaddr = self.email
		
		fp = open(self.tab_ofname, 'rb')
		msg = MIMEText(fp.read())
		fp.close()
		
		msg['Subject'] = 'Gene Prediction Results'
		msg['From'] = fromaddr
		msg['To'] = toaddr
		smtp_server.sendmail(fromaddr, toaddr, msg.as_string())
		smtp_server.close()

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema=", "table=", "mcl_table=", "p_value_cut_off=",\
		"connectivity_cut_off=", "recurrence_cut_off=", "email=", "report"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:m:p:n:y:e:r", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = 'sc_yh60_splat'
	table = 'cluster_stat_sup_3_2'
	mcl_table = 'mcl_result_sup_3_2'
	p_value_cut_off = 0.001
	connectivity_cut_off = 0.8
	recurrence_cut_off = 6
	email = ''
	report = 0
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-z", "--hostname"):
			hostname = arg
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-k", "--schema"):
			schema = arg
		elif opt in ("-t", "--table"):
			table = arg
		elif opt in ("-m", "--mcl_table"):
			mcl_table = arg
		elif opt in ("-p", "--p_value_cut_off"):
			p_value_cut_off = float(arg)
		elif opt in ("-n", "--connectivity_cut_off"):
			connectivity_cut_off = float(arg)
		elif opt in ("-y", "--recurrence_cut_off"):
			recurrence_cut_off = int(arg)
		elif opt in ("-e", "--email"):
			email = arg
		elif opt in ("-r", "--report"):
			report = 1
			
	if len(args)==2:
		instance = predict_genes(args[0], args[1], hostname, dbname, schema, table, mcl_table, \
			p_value_cut_off, connectivity_cut_off, recurrence_cut_off, email, report)
		instance.dstruc_loadin()
		instance.run()
	else:
		print __doc__
		sys.exit(2)
