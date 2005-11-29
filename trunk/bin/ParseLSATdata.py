#!/usr/bin/env python
"""
Usage: ParseLSATdata.py abstract_ens040505_forweb.txt abstract_tissue.txt abstract_isoforminfo.txt

Description:
	Program to parse data downloaded from LSAT.
	Output ensembl_id2tissue and ensembl_id2no_of_isoforms.
"""

import sys, os, csv, re
sys.path += [os.path.expanduser('~/script/annot/bin')]
from sets import Set

"""
11-29-05
	for abstract_ens040505_forweb.txt
"""
def get_ensembl_id2pmid(inputfname):
	ensembl_id2pmid_set = {}
	reader = csv.reader(open(inputfname, 'r'), delimiter='|')
	for row in reader:
		pmid, ensembl_id, organism = row
		pmid = int(pmid)
		dot_index = ensembl_id.find('.')
		if dot_index!=-1:
			ensembl_id = ensembl_id[:dot_index]
		if ensembl_id not in ensembl_id2pmid_set:
			ensembl_id2pmid_set[ensembl_id] = Set()
		ensembl_id2pmid_set[ensembl_id].add(pmid)
	del reader
	return ensembl_id2pmid_set

"""
11-29-05
	for abstract_tissue.txt
"""
def get_pmid2tissue_set(inputfname):
	pmid2tissue_set = {}
	reader = csv.reader(open(inputfname, 'r'), delimiter='|')
	for row in reader:
		pmid, tissue = row
		pmid = int(pmid)
		tissue = tissue.strip()	#remove beginning and ending spaces
		tissue = tissue.lower()	#lower case
		if pmid not in pmid2tissue_set:
			pmid2tissue_set[pmid] = Set()
		pmid2tissue_set[pmid].add(tissue)
	del reader
	return pmid2tissue_set


"""
11-29-05
	for abstract_isoforminfo.txt
"""
def get_pmid2no_of_isoforms_ls(inputfname):
	pmid2no_of_isoforms_ls = {}
	literal2no = {'one':1,
		'two':2,
		'three':3,
		'four':4,
		'five':5,
		'six':6,
		'seven':7,
		'eight':8,
		'first':1,
		'second':2,
		'third':3,
		'fourth':4}
	reader = csv.reader(open(inputfname, 'r'), delimiter='|')
	for row in reader:
		pmid, no_of_isoforms = row
		pmid = int(pmid)
		no_of_isoforms = no_of_isoforms.strip()	#remove beginning and ending spaces
		no_of_isoforms  = no_of_isoforms.lower()	#lower case
		if no_of_isoforms not in literal2no:
			print no_of_isoforms, "not in literal2no"
			break
		if pmid not in pmid2no_of_isoforms_ls:
			pmid2no_of_isoforms_ls[pmid] = []
		pmid2no_of_isoforms_ls[pmid].append(literal2no[no_of_isoforms])
	del reader
	return pmid2no_of_isoforms_ls

"""
11-29-05 convert pmid2no_of_isoforms_ls above into pmid2avg_no_of_isoforms
"""
def get_pmid2avg_no_of_isoforms(pmid2no_of_isoforms_ls):
	pmid2avg_no_of_isoforms = {}
	for pmid, no_of_isoforms_ls in pmid2no_of_isoforms_ls.iteritems():
		avg_no_of_isoforms = sum(no_of_isoforms_ls)/float(len(no_of_isoforms_ls))
		pmid2avg_no_of_isoforms[pmid] = avg_no_of_isoforms
	return pmid2avg_no_of_isoforms

"""
11-29-05
	for any file the 1st column is pmid(almost all)
"""
def get_pmid_set(inputfname):
	pmid_set = Set()
	reader = csv.reader(open(inputfname, 'r'), delimiter='|')
	for row in reader:
		pmid = int(row[0])
		pmid_set.add(pmid)
	del reader
	return pmid_set
	
if __name__ == '__main__':
	
	import sys,os
	sys.path += [os.path.expanduser('~/script/annot/bin')]
	if len(sys.argv)==1:
		print __doc__
		#print "\t parsed_alternative_polyadenylation.txt parsed_alternative_promoter.txt parsed_alternative_splicing.txt"
		sys.exit(0)
	ensembl_id2pmid_set_infname, pmid2tissue_set_infname, pmid2no_of_isoforms_ls_infname = sys.argv[1:]
		#pmid_set_ap_infname, pmid_set_dp_infname, pmid_set_as_infname	=  sys.argv[1:]
	
	ensembl_id2pmid_set = get_ensembl_id2pmid(ensembl_id2pmid_set_infname)
	pmid2tissue_set = get_pmid2tissue_set(pmid2tissue_set_infname)
	pmid2no_of_isoforms_ls = get_pmid2no_of_isoforms_ls(pmid2no_of_isoforms_ls_infname)
	pmid2avg_no_of_isoforms = get_pmid2avg_no_of_isoforms(pmid2no_of_isoforms_ls)
	"""
	pmid_set_ap = get_pmid_set(pmid_set_ap_infname)
	pmid_set_dp = get_pmid_set(pmid_set_dp_infname)
	pmid_set_as = get_pmid_set(pmid_set_as_infname)
	"""
	#output the ensembl_id2tissue
	ensembl_id2tissue_outf = open('ensembl_id2tissue','w')
	for ensembl_id, pmid_set in ensembl_id2pmid_set.iteritems():
		for pmid in pmid_set:
			if pmid in pmid2tissue_set:
				for tissue in pmid2tissue_set[pmid]:
					ensembl_id2tissue_outf.write('%s\t%s\n'%(ensembl_id, tissue))
	ensembl_id2tissue_outf.close()
	
	ensembl_id2no_of_isoforms = {}
	for ensembl_id, pmid_set in ensembl_id2pmid_set.iteritems():
		no_of_isoforms_ls = []
		for pmid in pmid_set:
			if pmid in pmid2avg_no_of_isoforms:
				no_of_isoforms_ls.append(pmid2avg_no_of_isoforms[pmid])
		if len(no_of_isoforms_ls)>0:
			ensembl_id2no_of_isoforms[ensembl_id] = sum(no_of_isoforms_ls)/float(len(no_of_isoforms_ls))
		else:
			ensembl_id2no_of_isoforms[ensembl_id] = 2
	
	ensembl_id2no_of_isoforms_outf = open('ensembl_id2no_of_isoforms', 'w')
	for ensembl_id, no_of_isoforms in ensembl_id2no_of_isoforms.iteritems():
		ensembl_id2no_of_isoforms_outf.write('%s\t%s\n'%(ensembl_id, no_of_isoforms))
	del ensembl_id2no_of_isoforms_outf
