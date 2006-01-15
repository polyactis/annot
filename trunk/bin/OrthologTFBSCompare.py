#!/usr/bin/env python
"""
Usage: OrthologTFBSCompare.py -p output_prefix [OPTIONS]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database(IGNORE)
	-i ...,	ReformatTRANSFACOutput_fname list, separated by ','
	-p ...,	output prefix
	-y ...,	type, 1(gene_id2mt_no_set from transfac.binding_site, default),
		2(from ReformatTRANSFACOutput_fname)
	-b	enable debugging, no debug by default
	-r	report the progress(a number)
	-h, --help              show this help

Examples:
	OrthologTFBSCompare.py -p /tmp/yuhuang/hs_mm_tfbs_homo -r

Description:
	Program to draw histogram of tfbs_similarity_ls between two groups
		of homologenes.
"""
import os, sys, getopt, random, csv, fileinput
sys.path += [os.path.join(os.path.expanduser('~/script/annot/bin'))]
from codense.common import db_connect
from sets import Set
from rpy import r

class OrthologTFBSCompare:
	"""
	01-05-06
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, input_fname_list=None, \
		output_prefix=None, type=1, debug=0, report=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_fname_list = input_fname_list
		self.output_prefix = output_prefix
		self.type = int(type)
		self.debug = int(debug)
		self.report = int(report)
	
	def get_prom_id2gene_id(self, curs, prom_seq_table='transfac.prom_seq'):
		"""
		01-06-06
			promoter id to Entrez gene id
			only upstream(prom_type_id=1)
		"""
		sys.stderr.write("Getting prom_id2gene_id...\n")
		prom_id2gene_id = {}
		curs.execute("DECLARE crs CURSOR FOR select id, prom_acc from %s where prom_type_id=1"%prom_seq_table)
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		counter = 0
		while rows:
			for row in rows:
				prom_id, gene_id = row
				gene_id = int(gene_id)
				prom_id2gene_id[prom_id] = gene_id
				counter += 1
			if self.report:
				sys.stderr.write("%s%s"%('\x08'*20, counter))
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		curs.execute("close crs")
		sys.stderr.write("Done.\n")
		return prom_id2gene_id
	
	def get_mt_id2no(self, curs, matrix_table='transfac.matrix'):
		"""
		01-06-06
		"""
		sys.stderr.write("Getting mt_id2no...\n")
		curs.execute("select mt_id, id from %s"%matrix_table)
		mt_id2no = {}
		rows = curs.fetchall()
		for row in rows:
			mt_id, no = row
			mt_id2no[mt_id] = no
		sys.stderr.write("Done.\n")
		return mt_id2no
	
	def get_gene_id2mt_no_set_from_file(self, ReformatTRANSFACOutput_fname_list, mt_id2no, prom_id2gene_id):
		"""
		01-06-06
		"""
		sys.stderr.write("Getting gene_id2mt_no_set from ReformatTRANSFACOutput_fname...\n")
		aggregated_inf = fileinput.input(ReformatTRANSFACOutput_fname_list)
		reader = csv.reader(aggregated_inf, delimiter='\t')
		gene_id2mt_no_set = {}
		counter = 0
		inner_counter = 0
		for row in reader:
			mt_id, prom_id, strand, bs_disp_start, bs_disp_end, core_similarity_score, matrix_similarity_score, sequence, order = row
			prom_id = int(prom_id)
			if prom_id in prom_id2gene_id:
				gene_id = prom_id2gene_id[prom_id]
				mt_no = mt_id2no[mt_id]
				if gene_id not in gene_id2mt_no_set:
					gene_id2mt_no_set[gene_id] = Set()
				gene_id2mt_no_set[gene_id].add(mt_no)
				inner_counter += 1
			counter += 1
			if counter%100000==0:
				sys.stderr.write("%s%s/%s"%('\x08'*20, counter, inner_counter))
		sys.stderr.write("Done.\n")
		return gene_id2mt_no_set
	
	def get_gene_id2mt_no_set(self, curs, binding_site_table='transfac.binding_site', \
		prom_seq_table='transfac.prom_seq', matrix_table='transfac.matrix'):
		sys.stderr.write("Getting gene_id2mt_no_set...\n")
		curs.execute("DECLARE crs CURSOR  FOR select p.prom_acc, m.id from %s b, %s p, %s m \
			where p.id=b.prom_id  and b.mt_id=m.mt_id and p.prom_type_id=1"%(binding_site_table, \
			prom_seq_table, matrix_table))
			#only upstream
		gene_id2mt_no_set = {}
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		counter = 0
		while rows:
			for row in rows:
				gene_id, mt_no = row
				gene_id = int(gene_id)
				if gene_id not in gene_id2mt_no_set:
					gene_id2mt_no_set[gene_id] = Set()
				gene_id2mt_no_set[gene_id].add(mt_no)
				counter += 1
			if self.report:
				sys.stderr.write("%s%s"%('\x08'*20, counter))
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		curs.execute("close crs")
		sys.stderr.write("Done.\n")
		return gene_id2mt_no_set
	
	def get_homo_genes_for_2_tax_ids(self, homo_group, tax_id1=9606, tax_id2=10090):
		"""
		01-06-06, just hs and mm
		"""
		tax_id1_gene_id_list = []
		tax_id2_gene_id_list = []
		for tax_id, gene_id in homo_group:
			if tax_id==tax_id1:
				tax_id1_gene_id_list.append(gene_id)
			elif tax_id==tax_id2:
				tax_id2_gene_id_list.append(gene_id)
		return tax_id1_gene_id_list, tax_id2_gene_id_list
	
	def get_ortholog_pairs(self, curs, bad_ortholog_fname, homologene_table='homologene.homologene'):
		"""
		01-06-06, just hs and mm
		01-06-06, fix a small bug, not crucial for histograms. replace hid with prev_hid when
			reporting non 1to1 ortholog pairs
		01-11-06
			add bad_ortholog_fname
		"""		
		sys.stderr.write("Getting ortholog pairs...\n")
		curs.execute("select hid, tax_id, gene_id from %s order by hid"%homologene_table)
		rows = curs.fetchall()
		bad_ortholog_f = open(bad_ortholog_fname, 'w')
		prev_hid = None
		hs_gene_id_list = []
		mm_gene_id_list = []
		homo_group = []
		tax_id1 = 9606
		tax_id2 = 10090
		for row in rows:
			hid, tax_id, gene_id = row
			if prev_hid == None:
				prev_hid = hid
			if hid!=prev_hid:
				tax_id1_gene_id_list, tax_id2_gene_id_list = \
					self.get_homo_genes_for_2_tax_ids(homo_group, tax_id1, tax_id2)
				if tax_id1_gene_id_list and tax_id2_gene_id_list:
					hs_gene_id_list.append(tax_id1_gene_id_list[0])	#only the 1st one
					mm_gene_id_list.append(tax_id2_gene_id_list[0])
				if len(tax_id1_gene_id_list)!=1 or len(tax_id2_gene_id_list)!=1:
					bad_ortholog_f.write("hid %s has %s tax_id=%s genes, %s tax_id=%s genes.\n"%(prev_hid, \
						len(tax_id1_gene_id_list), tax_id1, len(tax_id2_gene_id_list), tax_id2))	#here it is prev_hid
				homo_group = []
				prev_hid = hid
			homo_group.append((tax_id, gene_id))
		bad_ortholog_f.close()
		sys.stderr.write("Done.\n")
		return hs_gene_id_list, mm_gene_id_list
	
	def calculate_tfbs_similarity(self, mt_no_set1, mt_no_set2):
		return float(len(mt_no_set1&mt_no_set2))/len(mt_no_set1|mt_no_set2)
		
	def get_tfbs_similarity_ls(self, hs_gene_id_list, mm_gene_id_list, gene_id2mt_no_set):
		sys.stderr.write("Calculating tfbs_similarity.\n")
		tfbs_similarity_ls = []
		for i in range(len(hs_gene_id_list)):
			hs_gene_id = hs_gene_id_list[i]
			mm_gene_id = mm_gene_id_list[i]
			if hs_gene_id not in gene_id2mt_no_set or mm_gene_id not in gene_id2mt_no_set:
				continue
			tfbs_similarity_ls.append(self.calculate_tfbs_similarity(gene_id2mt_no_set[hs_gene_id], \
				gene_id2mt_no_set[mm_gene_id]))
		sys.stderr.write("Done.\n")
		return tfbs_similarity_ls
	
	def draw_tfbs_similarity_ls_histogram(self, tfbs_similarity_ls, output_fname):
		sys.stderr.write("Drawing histogram for tfbs_similarity_ls...")
		if len(tfbs_similarity_ls)>10:
			r.png('%s'%output_fname)
			r.hist(tfbs_similarity_ls, main='histogram',xlab='tfbs_similarity',ylab='freq')
			r.dev_off()
			sys.stderr.write("Done.\n")
		else:
			sys.stderr.write("too short: %s, aborted\n"%tfbs_similarity_ls)
	
	def run(self):
		"""
		01-06-06
		01-14-06
			add kolmogorov-smirnov test
		
			--db_connect()
			--get_gene_id2mt_no_set()
			--get_ortholog_pairs()
				--get_homo_genes_for_2_tax_ids()
			--get_tfbs_similarity_ls()
				--calculate_tfbs_similarity()
			--draw_tfbs_similarity_ls_histogram()
		"""
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		if self.type==1:
			gene_id2mt_no_set = self.get_gene_id2mt_no_set(curs)
		elif self.type==2:
			if self.input_fname_list==[]:
				sys.stderr.write("Please specify -i, type=2 needs that.\n")
				sys.exit(2)
			prom_id2gene_id = self.get_prom_id2gene_id(curs)
			mt_id2no = self.get_mt_id2no(curs)
			gene_id2mt_no_set = self.get_gene_id2mt_no_set_from_file(self.input_fname_list, mt_id2no, prom_id2gene_id)
		#get ortholog pairs
		bad_ortholog_fname = '%s.bad_ortholog'%self.output_prefix
		hs_gene_id_list, mm_gene_id_list = self.get_ortholog_pairs(curs, bad_ortholog_fname)
		
		tfbs_similarity_ls = self.get_tfbs_similarity_ls(hs_gene_id_list, mm_gene_id_list, gene_id2mt_no_set)
		ortholog_pair_hist_fname = '%s_ortholog_pair_hist.png'%self.output_prefix
		self.draw_tfbs_similarity_ls_histogram(tfbs_similarity_ls, ortholog_pair_hist_fname)
		
		#random
		random.shuffle(mm_gene_id_list)
		random_tfbs_similarity_ls = self.get_tfbs_similarity_ls(hs_gene_id_list, mm_gene_id_list, gene_id2mt_no_set)
		random_pair_hist_fname = '%s_random_pair_hist.png'%self.output_prefix
		self.draw_tfbs_similarity_ls_histogram(random_tfbs_similarity_ls, random_pair_hist_fname)
		
		#do a kolmogorov-smirnov test
		ks_result = r.ks_test(tfbs_similarity_ls,random_tfbs_similarity_ls)
		print "kolmogorov-smirnov test p-value:", ks_result['p.value']

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["help", "hostname=", "dbname=", "schema="]
	opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:p:y:br", long_options_list)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = None
	input_fname = None
	output_prefix = None
	type = 1
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
			input_fname = arg.split(',')
		elif opt in ("-p",):
			output_prefix = arg
		elif opt in ("-y",):
			type = int(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-r",):
			report = 1
		
	if output_prefix:
		instance = OrthologTFBSCompare(hostname, dbname, schema, input_fname, output_prefix, \
			type, debug, report)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
