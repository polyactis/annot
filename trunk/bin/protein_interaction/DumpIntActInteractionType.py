#!/usr/bin/env python
"""
Usage: DumpIntActInteractionType.py -i INPUTFILE [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database, protein_interaction(default)
	-i ...,	input filename
	-t ...,	table, intact_interaction_type(default)
	-c,	commit
	-b,	debug version
	-h, --help              show this help
	
Examples:
	DumpIntActInteractionType.py -i CvInteractionType.def -c

Description:
	Dump CvInteractionType.def from IntAct into database.
	
"""

import sys, os, getopt, re
sys.path += [os.path.expanduser('~/script/annot/bin')]
sys.path += [os.path.expanduser('..')]
from codense.common import db_connect
if sys.version_info[:2] < (2, 3):       #python2.2 or lower needs some extra
		from python2_3 import *

class interaction_type_attribute:
	def __init__(self):
		self.is_a = []
		self.tag2variable = {
			'name': None,
			'id': None,
			'def': None,
			'exact_synonym': None,
			'is_obsolete': 0,
			}

class DumpIntActInteractionType:
	def __init__(self, hostname='zhoudb', dbname='mdb', schema='protein_interaction', \
		input_fname=None, table='intact_interaction_type', commit=0, debug=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.input_fname = input_fname
		self.table = table
		self.commit = int(commit)
		self.debug = int(debug)
		self.regex = re.compile('"(.*)" \[.*\]')
		self.isaregex = re.compile(' ! .*')
	
	def submit2table(self, curs, table, interaction_type_attrib):
		"""
		12-29-05
			solve the '-within-text problem
			change acc to interaction_type_id
		"""
		curs.execute("insert into " + table +" (interaction_type_id, name, definition, is_obsolete, exact_synonym) \
			values(%s, %s, %s, %s, %s)", (interaction_type_attrib.tag2variable['id'],\
			interaction_type_attrib.tag2variable['name'], interaction_type_attrib.tag2variable['def'], \
			interaction_type_attrib.tag2variable['is_obsolete'], \
			interaction_type_attrib.tag2variable['exact_synonym']))
		for isa in interaction_type_attrib.is_a:
			curs.execute("insert into " + table+"_rel (interaction_type_id, is_a_interaction_type_id) \
				values(%s, %s)", (interaction_type_attrib.tag2variable['id'],isa))
	
	def parse_block(self, block):
		"""
		term: acetylation: acetylation reaction
		id: MI:0192
		definition: Reaction, that can affect K,C,A,D,E,Q,G,I,K,M,P,S,T,Y,V residues.
		definition_reference: RESID:AA0055 , AA0056 , AA0041 , AA0042 , AA0043 , AA0044
		, AA0045 , AA0046 , AA0047 , AA0048 , AA0049 , AA0050 , AA0051 , AA0052 , AA0053 , AA0054


		"""
		interaction_type_attrib = interaction_type_attribute()
		for line in block:
			colon_position = line.find(':')
			tag = line[:colon_position]
			value = line[colon_position+2:-1]	#+2 to skip ': '
			if tag == 'is_obsolete':
				if value.find('true') == 0:
					interaction_type_attrib.tag2variable[tag] = 1
			elif tag == 'is_a':
				value = self.isaregex.sub('',value)
				interaction_type_attrib.is_a.append(value)
			else:
				#value = value.replace("'", '"')	#replace ' with "
				value = self.regex.sub('\\1',value)
				#print "%s = %s"%(tag,value)
				interaction_type_attrib.tag2variable[tag] = value
		if self.debug:
			print interaction_type_attrib.tag2variable['id']
			print interaction_type_attrib.tag2variable['name']
			print interaction_type_attrib.tag2variable['def']
			print interaction_type_attrib.tag2variable['is_obsolete']
			for isa in interaction_type_attrib.is_a:
				print "is a %s"%(isa)
			raw_input("continue?(Y/n)")
		return interaction_type_attrib
	
	def parse(self, curs, input_fname, table):
		sys.stderr.write("Parsing...")
		inf = open(input_fname, 'r')
		block = []
		state = 0
		for line in inf:
			#print line
			if state==0:
				# searching for [Term]
				if line=='[Term]\n':
					state=1
			elif state==1:
				# Accepting text until newline
				if line=='\n':
					interaction_type_attrib = self.parse_block(block)
					self.submit2table(curs, table, interaction_type_attrib)
					block = []
					state=0
				else:
					block.append(line)
		del inf
		sys.stderr.write("done.\n")
	
	def run(self):
		"""
		12-28-05
		"""
		(conn, curs) =  db_connect(self.hostname, self.dbname, self.schema)
		self.parse(curs, self.input_fname, self.table)
		if self.commit:
			curs.execute("end")


if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	
	long_options_list = ["hostname=", "dbname=", "schema=", "help"]
	try:
		opts, args = getopt.getopt(sys.argv[1:], "z:d:k:i:t:cbh", long_options_list)
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = 'protein_interaction'
	#schema = 'mrinal_pi'
	input_fname = ''
	table = 'intact_interaction_type'
	commit = 0
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
		elif opt in ("-i",):
			input_fname = arg
		elif opt in ("-t",):
			table = arg
		elif opt in ("-b",):
			debug = 1
		elif opt in ("-c",):
			commit = 1

	if hostname and dbname and schema and input_fname:
		instance = DumpIntActInteractionType(hostname, dbname, schema, input_fname,\
			table, commit, debug)
		instance.run()
	else:
		print __doc__
		sys.exit(2)
