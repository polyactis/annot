#!/usr/bin/env python
"""
Usage: GO_no_parent_filter.py [OPTION] DATADIR Outputfile

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database, go(default)
	-t ..., --type=...	which branch of GO, 0, 1(default) or 2(IGNORE it)
	-v, --visualize	visualize the subgraph GraphDot
	-h, --help              show this help

Examples:
	GO_no_parent_filter.py /tmp/cfgs probe2go
	GO_no_parent_filter.py -d GO /tmp/cfgs/ probe2go
	
Description:
	Each file in DATADIR is Paul pavlidis' chip annotation files.
	THis program filters out those parent annotations.
"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import sys, os, psycopg, getopt, csv
from GO_graphxml import GO_graphxml
from graphlib import GraphDot

class GO_no_parent_filter(GO_graphxml):
	def __init__(self, hostname, dbname, schema, type, dir, output, visualize):
		GO_graphxml.__init__(self, hostname, dbname, schema, type, output)
		self.dir = dir
		self.visualize = visualize
	
	def batch_filter(self):
		#loadin three data structures
		self.dstruc_loadin()
		files = os.listdir(self.dir)
		sys.stderr.write("\tTotally, %d files to be processed.\n"%len(files))
		if self.visualize == 0:
			of = open(self.ofname, 'w')
		
		for f in files:
			sys.stderr.write("%d/%d:\t%s\n"%(files.index(f)+1,len(files),f))
			src_file = os.path.join(self.dir, f)
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
				#get the subgraph constituted by these nodes
				go_subgraph = self.go_graph.subgraph_from_node_list(node_list)
				#visualize it
				if self.visualize==1:
					dot = GraphDot.Dot(go_subgraph)
					dot.display()
					yes = raw_input("Continue? (y/n):\t")
					if yes == 'y':
						continue
					else:
						sys.exit(2)
				#get the terminal nodes
				else:
					for node in node_list:
						#terminal node in this subgraph is the real GO annotation, output the third column to indicate obsolete or not
						if go_subgraph.out_degree(node) == 0:
							of.write('%s\t%s\t%d\t%d\n'%(probe_id, self.termid_dict[node].acc, node, self.termid_dict[node].is_obsolete))


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)

	long_options_list = ["help", "hostname=", "dbname=", "schema=", "type=", "visualize"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:v", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = 'go'
	type = 1
	visualize = 0
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
		elif opt in ("-t", "--type"):
			type = int(arg)
		elif opt in ("-v", "--visualize"):
			visualize = 1

	if len(args) == 2:
		'''
		args[0] is the DATADIR
		args[1] is the Outputfile
		'''
		instance = GO_no_parent_filter(hostname, dbname, schema, type, args[0], args[1], visualize)
		instance.batch_filter()
	else:
		print __doc__
		sys.exit(2)
