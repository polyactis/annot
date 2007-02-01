#!/usr/bin/env python
"""
Usage: ConfirmWrongPredictionsWithUpdatedGOinfo.py -k SCHEMA -p xx -g xx -o xx
	-o xxx [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-a ...,	raw association table, 'gene.raw_association_2006_12_19'(default)
	-x ...,	term table, 'go200612.term'(default)
	-y ...,	term2term table, 'go200612.term2term'(default)
	-p ...,	p_gene table
	-g ...,	gene_p table
	-o ...,	output file
	-m ...,	organism, 'Homo sapiens'(default)
	-t ...,	prediction type, 0(default, wrong), -1(unknown), 1(correct)
	-b,	debug version.
	-h,	Display the usage infomation.
	
Examples:
	ConfirmWrongPredictionsWithUpdatedGOinfo.py -k hs_fim_65 -p p_gene_hs_fim_65_n2s175_m5x65s4l5_ft2_e5
	-g gene_p_hs_fim_65_n2s175_m5x65s4l5_ft2_e5_000001a60 -o /tmp/p_gene_hs_fim_65_n2s175_m5x65s4l5_ft2_e5.wrong.pred.confirm
	
	
Description:
	Program to confirm the wrong predictions with updated GO info.
	Check /research/annot/log-2006-10 for more info.
"""


import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import sys, os, getopt, csv, math
from codense.common import db_connect, pg_1d_array2python_ls
from sets import Set
from kjbuckets import *
from graphlib import Graph

class ConfirmWrongPredictionsWithUpdatedGOinfo:
	"""
	2006-12-19
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, \
		raw_association_table=None, term_table=None, term2term_table=None, p_gene_table=None,\
		gene_p_table=None, output_fname=None, organism='Homo sapiens', prediction_type=0, debug=0):
		"""
		2007-01-31 add prediction_type
		"""
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.raw_association_table = raw_association_table
		self.term_table = term_table
		self.term2term_table = term2term_table
		self.p_gene_table = p_gene_table
		self.gene_p_table = gene_p_table
		self.output_fname = output_fname
		self.organism = organism
		self.prediction_type = int(prediction_type)
		self.debug = int(debug)
		
		self.branch = 2	#default is biological_process, hidden parameter
		
		self.root_dict = {1:'GO:0003674',
			2:'GO:0008150',
			3:'GO:0005575'}
		self.root = self.root_dict[self.branch]
		self.branch_name_dict = {1:'molecular_function',
			2:'biological_process',
			3:'cellular_component'}
		self.branch_unknown_acc_dict = {1:'GO:0005554',
			2:'GO:0000004',
			3:'GO:0008372'}
			
	def dstruc_loadin_from_db(self, curs, raw_association_table, term_table, term2term_table, organism, \
		branch_name_dict, branch_unknown_acc_dict, branch):
		"""
		"""
		sys.stderr.write("Loading Data STructure...")
		direct_gene_id2go_id_set = {}
		#not obsolete, not unknown, biological_process, organism
		curs.execute("select a.go_id, a.gene_id from %s a, %s t\
			where a.organism='%s' and t.acc=a.go_id and\
			t.term_type='%s' and t.acc!='%s' and t.is_obsolete=0"%(raw_association_table, term_table,\
			organism, branch_name_dict[branch], branch_unknown_acc_dict[branch]))
		
		rows = curs.fetchall()
		for row in rows:
		#setup the go_id2gene_id_dict structure
			go_id = row[0]
			gene_id = int(row[1])	#2006-12-19 turn it into integer
			if gene_id not in direct_gene_id2go_id_set:
				direct_gene_id2go_id_set[gene_id] = Set()
			direct_gene_id2go_id_set[gene_id].add(go_id)
		
		#get the non-obsolete biological_process GO DAG
		curs.execute("select t2t.term1_id, t2t.term2_id, t1.acc, t2.acc from \
			%s t2t, %s t1, %s t2 where t2t.term1_id=t1.id and \
			t2t.term2_id=t2.id and t1.is_obsolete=0 and t2.is_obsolete=0 and \
			t1.term_type='%s' and t2.term_type='%s' "%(term2term_table, term_table, term_table, \
			branch_name_dict[branch], branch_name_dict[branch]))
		
		go_graph = kjGraph()
		complement_go_graph = kjGraph()	#used to check parents of a node
		
		rows = curs.fetchall()
		for row in rows:
		#setup the go_graph structure
			go_graph.add((row[2], row[3]))
			complement_go_graph.add((row[3], row[2]))
		
		#setup go_id2go_name
		go_id2go_name = {}
		curs.execute("select acc, name from %s"%term_table)
		rows = curs.fetchall()
		for row in rows:
			go_id2go_name[row[0]] = row[1]
		
		sys.stderr.write("Done\n")
		return direct_gene_id2go_id_set, go_graph, complement_go_graph, go_id2go_name
	
	def get_wrong_predictions(self, curs, p_gene_table, gene_p_table, prediction_type=0):
		"""
		2006-12-19
		2007-01-31
			add prediction_type
		"""
		sys.stderr.write("Getting wrong predictions ...")
		curs.execute("select p.gene_no, go.go_id from %s p, %s g, go where p.p_gene_id=g.p_gene_id and \
			p.go_no=go.go_no and p.is_correct=%s"%(p_gene_table, gene_p_table, prediction_type))
		rows = curs.fetchall()
		
		gene_id2go_id_set = {}
		for row in rows:
			gene_no, go_id = row
			if gene_no not in gene_id2go_id_set:
				gene_id2go_id_set[gene_no] = Set()
			gene_id2go_id_set[gene_no].add(go_id)
		no_of_total_wrong_predictions = len(rows)
		sys.stderr.write("Done\n")
		return gene_id2go_id_set, no_of_total_wrong_predictions
	
	def checking_wrong_predictions(self, go_graph, pred_gene_id2go_id_set, direct_gene_id2go_id_set, \
		no_of_total_wrong_predictions, go_id2go_name, output_fname):
		"""
		2006-12-20
			fix a bug (direct_go_id==go_id)
		"""
		sys.stderr.write("Checking wrong predictions ...")
		writer = csv.writer(open(output_fname, 'w'), delimiter='\t')
		no_of_correct_wrong_predictions = 0
		for gene_id, go_id_set in pred_gene_id2go_id_set.iteritems():
			if gene_id in direct_gene_id2go_id_set:
				direct_go_id_set = direct_gene_id2go_id_set[gene_id]
				for go_id in go_id_set:
					evidence_direct_go_id_list_for_this_prediction = []
					for direct_go_id in direct_go_id_set:
						if direct_go_id in go_graph.reachable(go_id).items() or direct_go_id==go_id:	#2006-12-20
							evidence_direct_go_id_list_for_this_prediction.append(direct_go_id)
					if len(evidence_direct_go_id_list_for_this_prediction)>0:
						no_of_correct_wrong_predictions += 1
						writer.writerow([gene_id, '%s(%s)'%(go_id, go_id2go_name[go_id]), evidence_direct_go_id_list_for_this_prediction])
					#else:
						#writer.writerow([gene_id, '%s'%(go_id), list(direct_go_id_set)])
		del writer
		sys.stderr.write("Done\n")
		return no_of_correct_wrong_predictions
		
	def run(self):
		"""
		2006-12-19
		2007-01-31
			get_wrong_predictions()
		"""
		conn, curs = db_connect(self.hostname, self.dbname, self.schema)
		direct_gene_id2go_id_set, go_graph, complement_go_graph, go_id2go_name = self.dstruc_loadin_from_db(curs,\
			self.raw_association_table, self.term_table, self.term2term_table, \
			self.organism, self.branch_name_dict, self.branch_unknown_acc_dict, self.branch)
		pred_gene_id2go_id_set, no_of_total_wrong_predictions = self.get_wrong_predictions(curs, self.p_gene_table, self.gene_p_table, self.prediction_type)
		no_of_correct_wrong_predictions = self.checking_wrong_predictions(go_graph, pred_gene_id2go_id_set, \
			direct_gene_id2go_id_set, no_of_total_wrong_predictions, go_id2go_name, self.output_fname)
		accuracy = no_of_correct_wrong_predictions/float(no_of_total_wrong_predictions)
		print "%s/%s=%s wrong predictions are correct."%(no_of_correct_wrong_predictions, no_of_total_wrong_predictions, accuracy)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:a:x:y:p:g:o:m:t:b", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	raw_association_table = 'gene.raw_association_2006_12_19'
	term_table = 'go200612.term'
	term2term_table = 'go200612.term2term'
	p_gene_table = None
	gene_p_table = None
	output_fname = None
	organism = 'Homo sapiens'
	prediction_type = 0
	debug = 0
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
		elif opt in ("-a",):
			raw_association_table = arg
		elif opt in ("-x",):
			term_table = arg
		elif opt in ("-y",):
			term2term_table = arg
		elif opt in ("-p",):
			p_gene_table = arg
		elif opt in ("-g",):
			gene_p_table = arg
		elif opt in ("-o",):
			output_fname = arg
		elif opt in ("-m",):
			organism = arg
		elif opt in ("-t",):
			prediction_type = int(arg)
		elif opt in ("-b",):
			debug = 1
	if schema and p_gene_table and gene_p_table and output_fname:
		instance = ConfirmWrongPredictionsWithUpdatedGOinfo(hostname, dbname, schema,  \
			raw_association_table, term_table, term2term_table, p_gene_table, gene_p_table, \
			output_fname, organism, prediction_type, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)