#!/usr/bin/env python
"""
Usage: DrawPatternExpr.py

Option:
        -z ..., --hostname=...        the hostname, zhoudb(default)
        -d ..., --dbname=...        the database name, graphdb(default)
        -k ..., --schema=...        which schema in the database
        -i ..., --inputfile=...        the input dir containing datasets
        -o ..., --outputfile=...        where to store graphs
        -c ..., --clusterbstable=...        cluster_bs table name
        -p ..., --patterntable=...        pattern table name
	-P id1,id2,..., --patterns=id1,id2,...	pattern ids
        
Examples:
        DrawPatternExpr.py -k hs_fim_65 -i ~/datasets/hs_fim_65/hs_dataset -o ~/graphs -P 97,98
        
Description:
        For the given patterns, draws graphs containing the expression levels
        If pattern list is not specified, draws graphs for all patterns
"""

RECURRENCE_THRESHOLD=0.8
TAXONOMY_ID=9606

import sys, os, getopt, csv
from codense.common import *
import matplotlib
from matplotlib.font_manager import *
from pylab import *

class DataSet:
	def __init__(self, filename):
		self.filename=filename
		self.expression_levels={}
		self.loaddata()
        
	def percentile(self, min, max, value):
		if (value=='NA'):
			return value
		else:
                	max = float(max)-float(min)
	                value = float(value)-float(min)
	                return (value/max)*100.0
	
	def get_expr_levels(self, geneID):
		return self.expression_levels[str(geneID)]

	def contains_gene(self, geneID):
		return self.expression_levels.has_key(str(geneID))
	
	def mymax(self, a, b):
		if a=='NA':
			return b
		elif b=='NA':
			return a
		else:
			return str(max(float(a),float(b)))
 
	def mymin(self, a, b):
		if a=='NA':
			return b
		elif b=='NA':
			return a
		else:
			return str(min(float(a),float(b)))
	
	def loaddata(self):
                reader = csv.reader(file(self.filename), delimiter='\t')
		counter = 0
                for line in reader:
			if counter==0:
				columnMax = (len(line)-1)*['NA']
				columnMin = (len(line)-1)*['NA']
				counter = 1
			else:
				columnMax = map(self.mymax, line[1:], columnMax)
        	                columnMin = map(self.mymin, line[1:], columnMin)
			self.expression_levels[line[0]]=line[1:]
                del reader
		for geneID,values in self.expression_levels.iteritems():
			self.expression_levels[geneID]=map(self.percentile, columnMin, columnMax, values)

class PatternExpressionGraphGenerator:
        def __init__(self, hostname=None, dbname=None, schema=None, indir=None, \
                outdir=None, clustertable=None, patterntable=None, patterns=None):
                self.hostname=hostname
                self.dbname=dbname
                self.schema=schema
                self.indir=indir
                self.outdir=outdir
                self.clustertable=schema+"."+clustertable
                self.patterntable=schema+"."+patterntable
		self.datasetcache={}
                self.patterns=patterns
                
	def removeNAs(self, values):
		index = 0
		indices=[]
		vals=[]
		for v in values:
			if v!='NA':
				indices.append(index)
				vals.append(v)
			index += 1
		return (indices,vals)

	def plotGeneSet(self, set, dataset, myTitle=None, xLabel=None, yLabel=None):
                for id in set:
			if (dataset.contains_gene(id)):
				(indices,values)=self.removeNAs(dataset.get_expr_levels(id))
        	                plot(indices, values, marker='o', markersize=5, label=self.gene_id_to_gene_symbol[id])
                legend(loc='best',prop=FontProperties(size=8))
		if not yLabel==None:
			ylabel(yLabel,fontsize=8)
		if not xLabel==None:
			xlabel(xLabel,fontsize=8)
		if not myTitle==None:
			title(myTitle,fontsize=10)
		locs,labels=xticks()
		setp(labels,fontsize=8)
		locs,labels=yticks()
		setp(labels,fontsize=8)
		#axis('tight')

		

        def generateGraphForDataSet(self, dataSetIndex, vertexSet, patternID, bsNoList):
		if not self.datasetcache.has_key(dataSetIndex):
			self.datasetcache[dataSetIndex] = DataSet(self.indir+"%s"%(dataSetIndex))
		dataset=self.datasetcache[dataSetIndex]
		subplot(211)
		self.plotGeneSet(bsNoList, dataset, yLabel='Expression level (percentile)', myTitle='Expression levels of bs_no_list of pattern ID %s in dataset %s'%(patternID,dataSetIndex))

		subplot(212)
		self.plotGeneSet(vertexSet, dataset, xLabel='Condition', yLabel='Expression level (percentile)',myTitle='Expression levels of genes of pattern ID %s in dataset %s'%(patternID,dataSetIndex))
                savefig('%s/pattern%s_dataset%s'%(self.outdir,patternID,dataSetIndex))
                clf()
        
        def generateGraphForPattern(self, patternID, vertexSet, recurrenceArray, bsNoList):
                index=0
                print "pattern id: %s"%(patternID),
		print "bslist = [",
		for bs in bsNoList:
			print "%s"%(bs),
                print "]"
                for rec in recurrenceArray:
                       index+=1
                       if (rec>RECURRENCE_THRESHOLD):
                               print "\tGenerating graph for dataset: %s (recurrence = %s)"%(index,rec)
                               self.generateGraphForDataSet(str(index), vertexSet, patternID, bsNoList)                              
                #for vertex in vertexSet:
                #        print "vertex: %s"%(vertex)
                
        def run(self):
                matplotlib.use('GD')
                self.conn, self.curs = db_connect(self.hostname, self.dbname)
		self.gene_id_to_gene_symbol = get_gene_id2gene_symbol(self.curs, TAXONOMY_ID)
                whereClause = ' '
                if not self.patterns==None:
                        whereClause = ' and pat.id in ('+self.patterns+') '
                #self.curs.execute("DECLARE crs CURSOR FOR select id,vertex_set,recurrence_array from %s"%(self.patterntable))
                self.curs.execute("DECLARE crs CURSOR FOR select pat.id,pat.vertex_set,pat.recurrence_array,clus.bs_no_list from %s as pat, %s as clus where pat.id=clus.mcl_id"%(self.patterntable,self.clustertable)+whereClause+"order by pat.id")
                self.curs.execute("fetch 1000 from crs")
                rows = self.curs.fetchall()
                counter = 0
		newRow=True
                while rows:
                        for row in rows:
                		patternID, vertexSet, recurrenceArray, bsNoList = row
				if newRow:
					newRow=False
					curPatternID=patternID
					curVertexSet=map(int,vertexSet[1:-1].split(','))
					curRecurrenceArray=map(float,recurrenceArray[1:-1].split(','))
					curBsNoList=map(int,bsNoList[1:-1].split(','))
				else:
					if patternID==curPatternID:
						#We're still processing the same id
						for newid in map(int,bsNoList[1:-1].split(',')):
							if (curBsNoList.count(newid)==0):
								curBsNoList.append(newid)
					else:
						newRow=True
                                		self.generateGraphForPattern(curPatternID, curVertexSet, curRecurrenceArray, curBsNoList)
                                counter +=1
                        self.curs.execute("fetch 1000 from crs")
                        rows = self.curs.fetchall()            
		self.generateGraphForPattern(curPatternID, curVertexSet, curRecurrenceArray, curBsNoList)

if __name__ == '__main__':
        if len(sys.argv) == 1:
                print __doc__
                sys.exit(2)        
        try:
                opts, args = getopt.getopt(sys.argv[1:], "hz:d:k:i:o:c:p:P:", ["help", "hostname=", \
                        "dbname=", "schema=", "inputdir=", "outputdir=", "clusterbstable=", "patterntable=", "patterns="])
        except:
                print __doc__
                sys.exit(2)
        
        hostname = 'zhoudb'
        dbname = 'graphdb'
        clustertable = 'cluster_bs_hs_fim_65_m5x65s4l5e0_1p001expt2'
        patterntable = 'pattern_hs_fim_65_m5x65s4l5'
        schema = None
        indir = None
        outdir = None
	patterns = None
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
                elif opt in ("-i", "--inputfile"):
                        indir = arg
                elif opt in ("-o", "--outputfile"):
                        outdir = arg
                elif opt in ("-c", "--clusterbstable"):
                        clustertable = arg
                elif opt in ("-p", "--patterntable"):
                        patterntable = arg
                elif opt in ("-P", "--patterns"):
                        patterns = arg
        if indir and outdir and schema:
                instance = PatternExpressionGraphGenerator(hostname, dbname, schema, indir, outdir, clustertable, patterntable, patterns)
                instance.run()
                
        else:
                print __doc__
                sys.exit(2)
