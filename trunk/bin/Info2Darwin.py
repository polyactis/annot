#!/usr/bin/env python
"""
Usage: Info2Darwin.py -o OUTPUTFILE [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database, graph(default)
	-t ...,	tax_id, 9606(default, human)
	-o ...,	output filename
	-l ...,	running bit, 1111111(default, all on)
	-b,	debug version
	-h, --help              show this help
	
Examples:
	#only family_size info
	Info2Darwin.py -o -l 0000001 /tmp/yuhuang/graphdb_info.darwin

Description:
	Convert following info of an organism from database to darwin format.
	
	NOTE:
		If the gene_symbol corresponding to a gene_id is not found in my local database, 
			i juse use gene_id as gene_symbol.
		
	1. GeneOntology biological_process direct association.
		go_bp:=[
			[gene_symbol, [bp1, bp2, ... ] ],
			...
			[]]:
	2. GeneOntology cellular_component direct association.
		go_cc:=[
			[gene_symbol, [cc1, cc2, ... ] ],
			...
			[]]:
	3. alternative splicing info(EBI's asd)
		as:=[
			[gene_symbol, no_of_events],
			...
			[]]:
	4. differential promoter(from Kim2005, PI: Ren, Bing)
		dp:=[
			[gene_symbol, no_of_promoters],
			...
			[]]:
	5. gene age(from NCBI Homologene)
		gene_age:=[
			[gene_symbol, [common_ancestor_depth:tax_id:short_organism, ..., ...] ],
			...
			[]]:
			
			common_ancestor_depth is based on the taxonomy tree.
			tax_id is searchable from NCBI taxonomy database.
			short_organism is the initial of the scientific name of the organism.
	6. tissue info(from LSAT)
		gene_tissue:=[
			[gene_symbol, [tissue1, tissue2, ... ] ],
			...
			[]]:
	7. gene_family_size(from ENSEMBL)
		gene_family_size:=[
			[gene_symbol, family_size],
			...
			[]]:
"""

import sys, os, getopt, re
sys.path += [os.path.expanduser('~/script/annot/bin')]
from codense.common import db_connect, get_gene_id2gene_symbol, \
	get_gene_id2no_of_events, get_org_from_tax_id, dict_map, get_gene_id2family_size
from sets import Set
if sys.version_info[:2] < (2, 3):       #python2.2 or lower needs some extra
		from python2_3 import *

def get_tg_tax_id2ca_depth_tax_id_short_org(curs, src_tax_id, tax_id_relationship_table = 'homologene.tax_id_relationship'):
	"""
	12-28-05
		for a tg_tax_id, get a string form for it, specifically for jasmine's darwin format
	"""
	sys.stderr.write("Getting tg_tax_id2ca_depth_tax_id_short_org ...")
	tg_tax_id2ca_depth_tax_id_short_org = {}
	curs.execute("SELECT tg_tax_id, common_ancestor_depth, taxonomy.tax_id2org(tg_tax_id) from %s \
		where src_tax_id=%s"%(tax_id_relationship_table, src_tax_id))
	rows = curs.fetchall()
	for row in rows:
		tg_tax_id, common_ancestor_depth, organism = row
		short_org = organism.split()
		short_org = '%s%s'%(short_org[0][0], short_org[1][0])
		short_org = short_org.lower()
		tg_tax_id2ca_depth_tax_id_short_org[tg_tax_id] = '%s:%s:%s'%(common_ancestor_depth, tg_tax_id, short_org)
	sys.stderr.write("Done.\n")
	return tg_tax_id2ca_depth_tax_id_short_org

def get_gene_id2go_term(curs, table='graph.raw_association', term_type='biological_process', organism='Homo sapiens'):
	"""
	12-28-05
		get gene_id2go_term (direct association) for a certain branch
	"""
	sys.stderr.write("Getting gene_id2go_term for %s..."%term_type)
	curs.execute("select r.gene_id, r.go_id, t.name from %s r, go.term t where  t.acc=r.go_id and \
		t.is_obsolete=0 and t.term_type='%s' and r.organism='%s'"%(table, term_type, organism))
	gene_id2go_term = {}
	rows = curs.fetchall()
	for row in rows:
		gene_id, go_id, go_name = row
		gene_id = int(gene_id)	#convert it to integer
		if gene_id not in gene_id2go_term:
			gene_id2go_term[gene_id] = []
		gene_id2go_term[gene_id].append('%s %s'%(go_id, go_name))
	
	sys.stderr.write("Done.\n")
	return gene_id2go_term

def get_gene_id2ortholog_tax_id_set(curs, src_tax_id, homologene_table='homologene.homologene'):
	sys.stderr.write("Getting gene_id2ortholog_tax_id_set ...")
	gene_id2ortholog_tax_id_set = {}
	curs.execute("SELECT hid, tax_id, gene_id from %s order by hid"%homologene_table)
	rows = curs.fetchall()
	prev_hid = None
	prev_src_gene_id_ls = []
	ortholog_tax_id_set = Set()
	for row in rows:
		hid, tax_id, gene_id = row
		if hid!=prev_hid:	#not the 1st row, prev_src_gene_id_ls not empty
			for gene_id in prev_src_gene_id_ls:
				gene_id2ortholog_tax_id_set[gene_id] = ortholog_tax_id_set
			prev_hid = hid
			prev_src_gene_id_ls = []
			ortholog_tax_id_set = Set()
		if tax_id == src_tax_id:
			prev_src_gene_id_ls.append(gene_id)
		ortholog_tax_id_set.add(tax_id)
	#don't forget the last block
	for gene_id in prev_src_gene_id_ls:
		gene_id2ortholog_tax_id_set[gene_id] = ortholog_tax_id_set
	sys.stderr.write("Done.\n")
	return gene_id2ortholog_tax_id_set

def get_gene_id2tissue_list(curs, tax_id, ensembl_id2tissue_table='graph.ensembl_id2tissue', \
	ensembl_id2gene_id_table='graph.ensembl_id2gene_id', schema=None):
	sys.stderr.write("Getting gene_id2tissue_list...")
	gene_id2tissue_list = {}
	curs.execute("SELECT e2.gene_id, e1.tissue from %s e1, %s e2 where e1.ensembl_id=e2.ensembl_id \
		and e2.tax_id=%s"%(ensembl_id2tissue_table, ensembl_id2gene_id_table, tax_id))
	rows = curs.fetchall()
	for row in rows:
		gene_id, tissue = row
		gene_id = int(gene_id)
		if gene_id not in gene_id2tissue_list:
			gene_id2tissue_list[gene_id] = []
		gene_id2tissue_list[gene_id].append(tissue)
	sys.stderr.write("Done.\n")
	return gene_id2tissue_list

class Info2Darwin:
	def __init__(self, hostname='zhoudb', dbname='mdb', schema='protein_interaction', \
		tax_id =9606, output_fname=None, running_bit='1111111', commit=0, debug=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.tax_id = int(tax_id)
		self.output_fname = output_fname
		self.running_bit = running_bit
		self.commit = int(commit)
		self.debug = int(debug)
	
	def dict2darwin(self, dc, dc_name, key_map, outf):
		"""
		12-28-05
			output dictionary into darwin format
		"""
		sys.stderr.write("Outputting %s..."%dc_name)
		outf.write("%s:=[\n"%dc_name)
		for key, value in dc.iteritems():
			if key in key_map:	#if key's not in key_map, just use plain key
				key = key_map[key]
			else:
				key = str(key)	#convert it to string form
			row = [key, value]
			outf.write(repr(row))
			outf.write(',\n')
		outf.write("[]]:\n")
		sys.stderr.write("Done.\n")
	
	def run(self):
		"""
		12-28-05
		"""
		conn, curs = db_connect(self.hostname, self.dbname, self.schema)
		organism = get_org_from_tax_id(curs, self.tax_id)
		#get the key_map
		gene_id2symbol = get_gene_id2gene_symbol(curs, self.tax_id)
		#open output here
		outf = open(self.output_fname, 'w')
		
		if len(self.running_bit)>=1 and self.running_bit[0] =='1':
			gene_id2go_bp_term = get_gene_id2go_term(curs, term_type='biological_process', organism=organism)
			self.dict2darwin(gene_id2go_bp_term, 'go_bp', gene_id2symbol, outf)
		if len(self.running_bit)>=2 and self.running_bit[1] =='1':
			gene_id2go_cc_term = get_gene_id2go_term(curs, term_type='cellular_component', organism=organism)
			self.dict2darwin(gene_id2go_cc_term, 'go_cc', gene_id2symbol, outf)
		if len(self.running_bit)>=3 and self.running_bit[2] =='1':
			gene_id2no_of_events = get_gene_id2no_of_events(curs, self.tax_id, ensembl2no_of_events_table='graph.ensembl2no_of_events')
			self.dict2darwin(gene_id2no_of_events, 'as', gene_id2symbol, outf)
		if len(self.running_bit)>=4 and self.running_bit[3] =='1':
			gene_id2no_of_promoters = get_gene_id2no_of_events(curs, self.tax_id, ensembl2no_of_events_table='graph.ensembl_id2no_of_promoters')
			self.dict2darwin(gene_id2no_of_promoters, 'dp', gene_id2symbol, outf)
		if len(self.running_bit)>=5 and self.running_bit[4] =='1':
			tg_tax_id2ca_depth_tax_id_short_org = get_tg_tax_id2ca_depth_tax_id_short_org(curs, self.tax_id)
			gene_id2ortholog_tax_id_set = get_gene_id2ortholog_tax_id_set(curs, self.tax_id, homologene_table='homologene.homologene')
			#convert gene_id2ortholog_tax_id_set to gene_id2ca_depth_tax_id_short_org_list
			gene_id2ca_depth_tax_id_short_org_list = {}
			for gene_id, ortholog_tax_id_set in gene_id2ortholog_tax_id_set.iteritems():
				ca_depth_tax_id_short_org_list = dict_map(tg_tax_id2ca_depth_tax_id_short_org, list(ortholog_tax_id_set))
				ca_depth_tax_id_short_org_list.sort()
				gene_id2ca_depth_tax_id_short_org_list[gene_id] = ca_depth_tax_id_short_org_list
			self.dict2darwin(gene_id2ca_depth_tax_id_short_org_list, 'gene_age', gene_id2symbol, outf)
		if len(self.running_bit)>=6 and self.running_bit[5] =='1':
			gene_id2tissue_list = get_gene_id2tissue_list(curs, self.tax_id)
			self.dict2darwin(gene_id2tissue_list, 'gene_tissue', gene_id2symbol, outf)
		if len(self.running_bit)>=7 and self.running_bit[6] =='1':
			gene_id2family_size = get_gene_id2family_size(curs, self.tax_id)
			self.dict2darwin(gene_id2family_size, 'gene_family_size', gene_id2symbol, outf)
		#close output
		outf.close()
	
	
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:o:l:b", ["help", "hostname=", \
			"dbname=", "schema="])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = 'graph'
	tax_id = 9606
	output_fname = None
	running_bit = '1111111'
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
		elif opt in ("-t",):
			tax_id = int(arg)
		elif opt in ("-o",):
			output_fname = arg
		elif opt in ("-l",):
			running_bit = arg
		elif opt in ("-b",):
			debug = 1
	if schema and output_fname and tax_id:
		instance = Info2Darwin(hostname, dbname, schema, tax_id, output_fname, running_bit, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
