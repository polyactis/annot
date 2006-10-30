#!/usr/bin/env python
"""
Usage: IntActXML2db.py [OPTIONS] -i INPUTFILE_LIST

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database, protein_interaction(default)
	-i ...,	input file list, ',' separated, or just one
	-a ...,	acc_file, gene2accession from NCBI Gene
	-x ...,	tax_id list, ',' delimited, just for acc_file
	-t ...,	interaction table, 'intact_interaction'(default)
	-e ...,	experiment table, 'intact_expt'(default)
	-c,	commit this database transaction
	-b,	enable debug
	-r,	enable report

Examples:
	IntActXML2db.py -i yeast_uetz-2000-2.xml  -c
	
Description:
	Dump IntAct interaction xml files into database.
	Before inserting an experiment, the program will check if the database 
		already has it or not.
	
"""

import sys, getopt, os
sys.path += [os.path.join(os.path.expanduser('~/script/annot/bin'))]
sys.path += [os.path.join(os.path.expanduser('..'))]
sys.path += [os.path.join(os.path.expanduser('../../../microarray/bin'))]
sys.path += [os.path.expanduser('~/script/microarray/bin')]
from MdbId2GeneId import MdbId2GeneId
from codense.common import db_connect, get_gene_symbol2gene_id
import cElementTree as ElementTree
from sets import Set

class interaction_attribute:
	def __init__(self):
		self.uniprot_id_array = []
		self.gene_id_array = []
		self.interaction_type_id = None
		self.expt_id_array = []
		self.tax_id = None
		self.is_cross_species = 0
		self.intact_id = None

class expt_attribute:
	def __init__(self):
		self.expt_id = None
		self.short_label = None
		self.full_name = None
		self.pubmed_id = None
		

class IntActXML2db:
	"""
	12-29-05
	"""
	def __init__(self,hostname='zhoudb', dbname='graphdb', schema=None,  input_fname_list=[],\
		acc_file=None, tax_id_list=[], interaction_table='intact_interaction', expt_table='intact_expt', \
		commit=0, debug=0, report=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_fname_list = input_fname_list
		self.acc_file = acc_file
		self.tax_id_list = tax_id_list
		self.interaction_table = interaction_table
		self.expt_table = expt_table
		self.commit = int(commit)
		self.debug = int(debug)
		self.report = int(report)
		self.gene_symbol2gene_id = {}
	
	def submit_expt_table(self, curs, expt_attrib, expt_table):
		"""
		12-29-05
			WATCH: a little different format to solve the '-within-text problem
			first check if this expt_id is already present or not
		"""
		curs.execute("select expt_id from %s where expt_id='%s'"%(expt_table, expt_attrib.expt_id))
		rows = curs.fetchall()
		if len(rows)==0:	#not in table yet
			curs.execute("insert into " + expt_table +"(expt_id, short_label, full_name, pubmed_id) \
				values(%s, %s, %s, %s)", (expt_attrib.expt_id, \
				expt_attrib.short_label, expt_attrib.full_name, expt_attrib.pubmed_id))
	
	def submit_interaction_table(self, curs, interaction_attrib, interaction_table):
		curs.execute("insert into %s(uniprot_id_array, interaction_type_id, expt_id_array, \
			tax_id, is_cross_species, intact_id, gene_id_array) values(ARRAY%s, '%s', ARRAY%s, %s, %s, '%s', ARRAY%s)"%(interaction_table,\
			repr(interaction_attrib.uniprot_id_array), interaction_attrib.interaction_type_id,\
			repr(interaction_attrib.expt_id_array), interaction_attrib.tax_id, interaction_attrib.is_cross_species,\
			interaction_attrib.intact_id, repr(interaction_attrib.gene_id_array)))
	
	def parse_intact_xml_file(self, curs, input_fname, expt_table, interaction_table, acc_tax_id2gene_id_list):
		"""
		12-28-05
			xmlns="net:sf:psidev:mi" causes a namespace header to be added for each element
		12-29-05
			add acc_tax_id2gene_id_list
		"""		
		sys.stderr.write("Parsing %s...\n"%input_fname)
		namespace = '{net:sf:psidev:mi}'
		for event, elem in ElementTree.iterparse(input_fname):
			if elem.tag == '%sentry'%namespace:
				interactor_intact_id2uniprot_id_tax_id = {}	#later used in parsing interactionList
				for sub_elem in elem:
					if sub_elem.tag == '%sexperimentList'%namespace:
						for expt_desc_elem in sub_elem:
							expt_attrib = expt_attribute()
							expt_attrib.expt_id = expt_desc_elem.get("id")
							expt_attrib.short_label = expt_desc_elem.findtext('%snames/%sshortLabel'%(namespace, namespace))
							expt_attrib.full_name = expt_desc_elem.findtext("%snames/%sfullName"%(namespace, namespace))
							pubmed_ref_elem = expt_desc_elem.find("%sbibref/%sxref/%sprimaryRef"%(namespace, namespace, namespace))
							if pubmed_ref_elem.get("db")=="pubmed":
								expt_attrib.pubmed_id = int(pubmed_ref_elem.get("id"))
							else:
								expt_attrib.pubmed_id = -1
							self.submit_expt_table(curs, expt_attrib, expt_table)
							expt_desc_elem.clear()	#release memory
					if sub_elem.tag == '%sinteractorList'%namespace:
						for interactor_elem in sub_elem:
							interactor_intact_id = interactor_elem.get("id")
							uniprot_id_elem = interactor_elem.find('%sxref/%sprimaryRef'%(namespace, namespace))
							uniprot_id = uniprot_id_elem.get("id")
							tax_id_elem = interactor_elem.find("%sorganism"%namespace)
							if not tax_id_elem:
								interactor_intact_id2uniprot_id_tax_id[interactor_intact_id] = (uniprot_id.upper(), None, None)
							else:
								tax_id = int(tax_id_elem.get("ncbiTaxId"))
								#12-29-05
								key = (uniprot_id.upper(), tax_id)
								gene_id_list = acc_tax_id2gene_id_list.get(key)
								if gene_id_list and len(gene_id_list)==1:
									interactor_intact_id2uniprot_id_tax_id[interactor_intact_id] = (uniprot_id.upper(), int(gene_id_list[0]), tax_id)
								else:
									#sys.stderr.write("\t Warning: %s gets entrez gene_id_list: %s\n"%(uniprot_id, gene_id_list))
									interactor_intact_id2uniprot_id_tax_id[interactor_intact_id] = (uniprot_id.upper(), None, tax_id)	#12-29-05 use None
									#10-25-06 mrinal
									if not self.gene_symbol2gene_id.has_key(str(tax_id)):
										print "Getting gene symbol mappings for tax_id",tax_id
										self.gene_symbol2gene_id[str(tax_id)] = get_gene_symbol2gene_id(curs, tax_id)
									gs2gid = self.gene_symbol2gene_id[str(tax_id)]
	
									names_elem = interactor_elem.find('%snames'%namespace)
									gene_names = Set()
									gene_ids = Set()
									for alias_elem in names_elem:
										if alias_elem.tag == '%salias'%namespace:
											alias_type = alias_elem.get("type")
											if alias_type =="gene name":
												gene_names.add(alias_elem.text)
											elif alias_type =="gene name synonym":
												gene_names.add(alias_elem.text)
									for gene_name in gene_names:
										if gs2gid.has_key(gene_name):
											gene_ids.add(gs2gid[gene_name])
									if len(gene_ids)==1:
										interactor_intact_id2uniprot_id_tax_id[interactor_intact_id] = (uniprot_id.upper(), int(gene_ids.pop()), tax_id)
									else:
										sys.stderr.write("\t Warning: Couldn't find gene ids for interactor id %s\n"%(interactor_intact_id))
							interactor_elem.clear()	#release memory
					if sub_elem.tag == "%sinteractionList"%namespace:
						skipped=0
						wrote=0
						for interaction_elem in sub_elem:
							interaction_attrib = interaction_attribute()
							interaction_attrib.expt_id_array = [expt_ref_elem.text \
								for expt_ref_elem in interaction_elem.find("%sexperimentList"%namespace)]
							interaction_attrib.interaction_type_id = \
								interaction_elem.find('%sinteractionType/%sxref/%sprimaryRef'%(namespace, namespace, namespace)).get("id")
							interaction_attrib.intact_id = interaction_elem.find("%sxref/%sprimaryRef"%(namespace, namespace)).get("id")
							for prot_part_elem in interaction_elem.find("%sparticipantList"%namespace):
								prot_intact_id = prot_part_elem.find("%sinteractorRef"%namespace).text
								#12-29-05
								uniprot_id, gene_id, tax_id = interactor_intact_id2uniprot_id_tax_id[prot_intact_id]
								interaction_attrib.uniprot_id_array.append(uniprot_id)
								
								if interaction_attrib.tax_id and tax_id!=interaction_attrib.tax_id:
									sys.stderr.write("\t Warning: interaction %s has >1 tax_id: %s, %s(ignored).\n"%\
										(interaction_attrib.intact_id, interaction_attrib.tax_id, tax_id))
									interaction_attrib.is_cross_species = 1	#interaction not just within one species
								else:
									interaction_attrib.tax_id = tax_id

								#10-25-06 (mrinal)
								if gene_id:
									interaction_attrib.gene_id_array.append(gene_id)
								else:
									sys.stderr.write("\t Warning: prot_intact_id %s doesn't have proper NCBI gene id.\n"%prot_intact_id)
									interaction_attrib.is_cross_species = 1 # tag it as bad
								
							if interaction_attrib.uniprot_id_array and interaction_attrib.gene_id_array:	#12-29-05 not empty
								self.submit_interaction_table(curs, interaction_attrib, interaction_table)
								wrote+=1
							else:
								skipped+=1
							interaction_elem.clear()	#release memory
					
					sub_elem.clear()	#release the sub_elem
		sys.stderr.write("Done.\n")
		print "wrote:",wrote,"skipped:",skipped
	
	def run(self):
		conn, curs = db_connect(self.hostname, self.dbname, self.schema)
		MdbId2GeneId_instance = MdbId2GeneId()
		acc_tax_id2gene_id_list = MdbId2GeneId_instance.setup_acc2gene_id(self.acc_file, Set(self.tax_id_list))
		for input_fname in self.input_fname_list:
			self.parse_intact_xml_file(curs, input_fname, self.expt_table, self.interaction_table, acc_tax_id2gene_id_list)
		if self.commit:
			curs.execute("end")


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:a:x:t:e:cbr", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = 'protein_interaction'
	#schema = 'mrinal_pi'
	input_fname_list = []
	acc_file = None
	tax_id_list = []
	interaction_table = 'intact_interaction'
	expt_table = 'intact_expt'
	commit = 0
	debug = 0
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
		elif opt in ("-i",):
			input_fname_list = arg.split(',')
		elif opt in ("-a",):
			acc_file = arg
		elif opt in ("-x",):
			tax_id_list = arg.split(',')
			tax_id_list = map(int, tax_id_list)
		elif opt in ("-t",):
			interaction_table = arg
		elif opt in ("-e",):
			expt_table = arg
		elif opt in ("-c",):
			commit = 1
		elif opt in ("-b",):
			debug = 1
		elif opt in ("-r",):
			report = 1
	if input_fname_list and schema and acc_file and tax_id_list:
		instance = IntActXML2db(hostname, dbname, schema, input_fname_list, acc_file, \
			tax_id_list, interaction_table, expt_table, commit, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
