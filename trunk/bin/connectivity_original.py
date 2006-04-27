#!/usr/bin/env python
"""
Usage:	connectivity_original.py -k SCHEMA -t TABLE [OPTION]

Option:
	-z ..., --hostname=...	the hostname, zhoudb(default)
	-d ..., --dbname=...	the database name, graphdb(default)
	-k ..., --schema=...	which schema in the database
	-t ..., --table=...	a 'patttern_xxx' -like table, get the vertex_set from this table
	-e ...,	'edge_sig_vector' filename
	-w ..., --min_weight=...	the minimum weight of an edge. 6(default)(IGNORE)
	-s ..., --dict_threshold=...	the threshold of the internal dictionary, 10,000,000(default) 26.5% of 6G (IGNORE)
	-b, --debug	enable debugging, no debug by default
	-c, --commit	commit the database transaction (IGNORE)
	-r, --report	report the progress(a number)
	-h, --help	show this help
	
Examples:
	connectivity_original.py -k sc_fim_50  -t pattern_table -e sc_fim_50_3.sig_vector

Description:
	This program computes the original connectivity of a cluster, which is the
	average connectivity of the vertex set of the cluster in each occurrent dataset.
	It'll pop out scatterplot and histograms in the end. (04-18-06)

"""
import sys, os, math
bit_number = math.log(sys.maxint)/math.log(2)
if bit_number>40:       #64bit
	sys.path.insert(0, os.path.expanduser('~/lib64/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script64/annot/bin')))
else:   #32bit
	sys.path.insert(0, os.path.expanduser('~/lib/python'))
	sys.path.insert(0, os.path.join(os.path.expanduser('~/script/annot/bin')))
import sys, os, psycopg, re, getopt, csv
from codense.common import db_connect
#04-18-06
from MpiClusterBsStat import MpiClusterBsStat
import MLab

class connectivity_original:
	"""
	04-01-05
	
	--run
		--db_connect
		--alter_table
		--data_fetch
			--_connectivity_original
				--get_edge
		--update_connectivity_original
	"""
	def __init__(self, hostname='zhoudb', dbname='graphdb', schema=None, \
		table=None, edge_table='edge_cor_vector', min_weight=6, \
		dict_threshold=10000000, debug=0, report=0, needcommit=0):
		self.hostname = hostname
		self.dbname = dbname
		self.schema = schema
		self.table = table
		self.edge_table = edge_table
		self.min_weight = int(min_weight)
		self.dict_threshold = int(dict_threshold)
		self.debug = int(debug)
		self.report = int(report)
		self.needcommit = int(needcommit)
	
		#a dictionary limited by dict_threshold, key is (vertex1,vertex2), ascending order, value is its weight.
		self.edge_dict = {}
		#mcl_id v.s. connectivity_original mapping.
		self.mcl_id2connectivity = {}
		#counter
		self.no_of_records = 0
	
	def fill_edge2encodedOccurrence(self, sig_vector_fname, total_vertex_set=None):
		"""
		04-18-06
			not realy encoded
		"""
		sys.stderr.write("Getting edge2encodedOccurrence...\n")
		from MpiFromDatasetSignatureToPattern import encodeOccurrenceBv
		edge2encodedOccurrence = {}
		reader = csv.reader(open(sig_vector_fname), delimiter='\t')
		no_of_datasets = 0
		counter = 0
		for row in reader:
			edge = row[:2]
			edge = map(int, edge)
			if edge[0] not in total_vertex_set or edge[1] not in total_vertex_set:
				continue
			edge.sort()	#04-06-06 in ascending order
			sig_vector = row[2:]
			sig_vector = map(int, sig_vector)
			if no_of_datasets==0:
				no_of_datasets = len(sig_vector)
			edge2encodedOccurrence[tuple(edge)] = sig_vector
		sys.stderr.write("Done.\n")
		del reader
		return edge2encodedOccurrence, no_of_datasets
	
	def data_fetch(self, curs, table, edge2encodedOccurrence, min_weight=6):
		"""
		04-01-05
			_connectivity_original
		04-18-06
			use edge2encodedOccurrence to replace edge_table
			mcl_id => id
			add recurrence_array
		"""
		sys.stderr.write("Computing cluster's density in its source graph...\n")
		curs.execute("DECLARE crs CURSOR FOR select id, vertex_set, recurrence_array, connectivity \
			from %s"%(table))
		curs.execute("fetch 5000 from crs")
		rows = curs.fetchall()
		connectivity_list = []
		original_avg_connectivity_list = []
		while rows:
			for row in rows:
				connectivity = row[3]
				connectivity_list.append(connectivity)
				original_avg_connectivity = self._connectivity_original(curs, row, edge2encodedOccurrence, min_weight)
				original_avg_connectivity_list.append(original_avg_connectivity)
				self.no_of_records+=1
			if self.report:
				sys.stderr.write('%s%s'%('\x08'*20, self.no_of_records))
			
			curs.execute("fetch 5000 from crs")
			rows = curs.fetchall()
		curs.execute("close crs")
		sys.stderr.write("Done\n")
		return connectivity_list, original_avg_connectivity_list
	
	def _connectivity_original(self, curs, row, edge2encodedOccurrence, min_weight):
		"""
		04-01-05
			input: mcl_id and vertex_set(in row)
			output: update its connectivity_original
		04-18-06
			use edge2encodedOccurrence to replace edge_table
			add recurrence_array
			calculate the connectivity in each occurrent dataset, and  then get the average as the original connectivity
			min_weight not used.
		"""
		mcl_id = row[0]
		if self.debug:
			print "mcl_id: %s"%mcl_id
		vertex_set = row[1][1:-1].split(',')
		vertex_set = map(int, vertex_set)
		#sort it first.
		vertex_set.sort()
		recurrence_array = row[2][1:-1].split(',')
		recurrence_array = map(float, recurrence_array)
		
		occurrence_vector = []	#index starts from 0
		for i in range(len(recurrence_array)):
			if recurrence_array[i] == 1:
				occurrence_vector.append(i)
		
		connectivity_vector = [0.0]*len(occurrence_vector)
		no_of_nodes = len(vertex_set)
		no_of_edges = 0
		for i in range(no_of_nodes):
			for j in range(i+1, no_of_nodes):
				if vertex_set[i] < vertex_set[j]:
					edge = (vertex_set[i], vertex_set[j])
				else:
					edge = (vertex_set[j], vertex_set[i])
				if edge in edge2encodedOccurrence:
					sig_vector = edge2encodedOccurrence[edge]
					if self.debug:
						print "edge", edge, "sig_vector", sig_vector
					for k in range(len(occurrence_vector)):
						if sig_vector[occurrence_vector[k]]==1:
							connectivity_vector[k] += 1
							if self.debug:
								print "edge (%s, %s) in dataset: %s"%(vertex_set[i], vertex_set[j], occurrence_vector[k])
		common_factor = no_of_nodes*(no_of_nodes-1)/2.0
		common_factor_func  = lambda x: x/common_factor
		connectivity_vector = map(common_factor_func,  connectivity_vector)
		connectivity = MLab.average(connectivity_vector)
		if self.debug:
			print "connectivity_vector", connectivity_vector
			print "connectivity: %s"%connectivity
			raw_input("WAIT:\t")
		return connectivity
	
	def get_edge(self, curs, edge, edge_table='edge_cor_vector', min_weight=6):
		"""
		04-01-05
			input: two vertices of an edge
			output: return None if no such edge, its weight if it exists.
		04-18-06
			defunct
		"""
		if edge in self.edge_dict:
			#the internal dictionary has it
			weight = self.edge_dict[edge]
			if weight>=min_weight:
				if self.debug:
					print "edge %s in internal dictionary with weight %s"%(repr(edge), weight)
				return weight
			else:
				if self.debug:
					print "edge %s in internal dictionary but with weight %s<%s"%(repr(edge), weight,min_weight)
				return None
		else:		
			#check database
			weight = self._get_edge(curs, edge, edge_table)
			if weight:
				if len(self.edge_dict)>=self.dict_threshold:
					#over the threshold, throw away the first item.
					self.edge_dict.popitem()
				self.edge_dict[edge] = weight
				
				if weight>=min_weight:
					#bigger than the min_weight
					if self.debug:
						print "edge %s in database with weight %s"%(repr(edge), weight)
					return weight
				else:
					if self.debug:
						print "edge %s in database but with weight %s<%s"%(repr(edge), weight, min_weight)
					return None
			else:
				if self.debug:
					print "edge %s not found"%(repr(edge))
				return None
	
	def _get_edge(self, curs, edge, edge_table='edge_cor_vector'):
		"""
		04-03-05
			split from get_edge(), to be class independent
		04-18-06
			defunct
		"""
		curs.execute("select sig_vector from %s where edge_name='{%s,%s}'"%(edge_table, edge[0], edge[1]))
		rows = curs.fetchall()
		if len(rows)>0:
			#edge is present
			edge_data = rows[0][0][1:-1]
			edge_data = edge_data.split(',')
			edge_data = map(int, edge_data)
			weight = sum(edge_data)
		else:
			weight = None
		return weight
	
	def alter_table(self, curs, table):
		"""
		04-01-05
			add a column, connectivity_original to the table
		04-18-06
			defunct
		"""
		sys.stderr.write("Adding connectivity_original to %s..."%table)
		curs.execute("alter table %s add connectivity_original float"%table)
		sys.stderr.write("Done.\n")
		
	def update_connectivity_original(self, curs, table, mcl_id2connectivity):
		"""
		04-01-05
			update connectivity_original in table
		04-18-06
			defunct
		"""
		sys.stderr.write("Database updating %s..."%table)
		for mcl_id,connectivity in mcl_id2connectivity.iteritems():
			curs.execute("update %s set connectivity_original=%s where mcl_id=%d"%\
			(table, connectivity, mcl_id))
		sys.stderr.write("Done.\n")
	
	def run(self):
		"""
		04-01-05
		04-18-06
			overhaulled
		"""
		(conn, curs) = db_connect(self.hostname, self.dbname, self.schema)
		MpiClusterBsStat_instance = MpiClusterBsStat()
		total_vertex_set = MpiClusterBsStat_instance.return_total_vertex_set(curs, self.table)
		edge2encodedOccurrence, no_of_datasets = self.fill_edge2encodedOccurrence(\
			self.edge_table, total_vertex_set)
		"""
		if self.needcommit:
			self.alter_table(curs, self.table)
		"""
		connectivity_list, original_avg_connectivity_list = self.data_fetch(curs, self.table, edge2encodedOccurrence, self.min_weight)
		
		import pylab
		print "scatter plot of connectivity vs. original_avg_connectivity"
		pylab.scatter(connectivity_list, original_avg_connectivity_list)
		pylab.show()
		print "histogram of connectivity_list"
		pylab.hist(connectivity_list, 200)
		pylab.show()
		print "histogram of original_avg_connectivity_list"
		pylab.hist(original_avg_connectivity_list, 200)
		pylab.show()
		"""
		import rpy
		rpy.r.postscript('scatter.ps')
		rpy.r.plot(connectivity_list, original_avg_connectivity_list, xlab='connectivity', ylab='original_avg_connectivity', main='scatter plot')
		rpy.r.dev_off()
		"""
		"""
		if self.needcommit:
			self.update_connectivity_original(curs, self.table, self.mcl_id2connectivity)
			curs.execute("end")
		"""
		
if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
		
	try:
		opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:t:e:w:s:brc", \
			["help", "hostname=", "dbname=", "schema=", "table=", "edge_table=", \
			"min_weight=", "dict_threshold=", "debug", "report", "commit"])
	except:
		print __doc__
		sys.exit(2)
	
	hostname = 'zhoudb'
	dbname = 'graphdb'
	schema = ''
	table = ''
	edge_table = ''
	min_weight = 6
	dict_threshold = 10000000
	debug = 0
	commit = 0
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
		elif opt in ("-e", "--edge_table"):
			edge_table = arg
		elif opt in ("-w", "--min_weight"):
			min_weight = int(arg)
		elif opt in ("-s", "--dict_threshold"):
			dict_threshold = int(arg)
		elif opt in ("-b", "--debug"):
			debug = 1
		elif opt in ("-c", "--commit"):
			commit = 1
		elif opt in ("-r", "--report"):
			report = 1


	if schema and table and edge_table:
		instance = connectivity_original(hostname, dbname, schema, table, \
			edge_table, min_weight, dict_threshold, debug, report, commit)
		instance.run()

	else:
		print __doc__
		sys.exit(2)
