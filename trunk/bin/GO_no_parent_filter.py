#!/usr/bin/env python
"""
Usage: GO_no_parent_filter.py [OPTION] DATADIR Outputfile

Option:
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database, go(default)
	-t ..., --type=...	which branch of GO, 0, 1(default) or 2(IGNORE it)
	-h, --help              show this help

Examples:
	GO_no_parent_filter.py /tmp/cfgs probe2go
	GO_no_parent_filter.py -d GO /tmp/cfgs/ probe2go
	
Description:
	Each file in DATADIR is Paul pavlidis' chip annotation files.
	THis program filters out those parent annotations.
"""

import sys, os, psycopg, getopt, csv
from GO_graphxml import GO_graphxml

class GO_no_parent_filter(GO_graphxml):
	def __init__(self, dbname, schema, type, output):
		GO_graphxml.__init__(self, dbname, schema, type, output)
		
	def batch_filter(self, dir):
		#loadin three data structures
		self.dstruc_loadin()
		files = os.listdir(dir)
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
		of = open(self.ofname, 'w')
		
		for f in files:
			sys.stderr.write("%d/%d:\t%s\n"%(files.index(f)+1,len(files),f))
			src_file = os.path.join(dir, f)
			reader = csv.reader(file(src_file), delimiter='\t')
			#throw away the first line
			reader.next()
			for row in reader:
					probe_id = row[0]
					#go term list is on the fourth column
					if row[3] == '':
						continue
					acc_list = row[3].split('|')
					node_list = []
					#transform the go term to termid,
					for acc in acc_list:
						#may not be obsolete, but still not exist
						if acc in self.acc2id_dict:
							node_list.append(self.acc2id_dict[acc])
					go_subgraph = self.go_graph.subgraph_from_node_list(node_list)
					for node in node_list:
						#terminal node in this subgraph is the real GO annotation, also not obsolete
						if self.go_graph.out_degree(node) == 0 and self.termid_dict[node].is_obsolete == 0:
							of.write('%s\t%s\n'%(probe_id, self.termid_dict[node].acc))


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)

	long_options_list = ["help", "dbname=", "schema=", "type="]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hd:k:t:", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	dbname = 'graphdb'
	schema = 'go'
	type = 1
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			print __doc__
			sys.exit(2)
		elif opt in ("-d", "--dbname"):
			dbname = arg
		elif opt in ("-k", "--schema"):
			schema = arg
		elif opt in ("-t", "--type"):
			type = int(arg)

	if len(args) == 2:
		'''
		args[0] is the DATADIR
		args[1] is the Outputfile
		'''
		instance = GO_no_parent_filter(dbname, schema, type, args[1])
		instance.batch_filter(args[0])
	else:
		print __doc__
		sys.exit(2)
