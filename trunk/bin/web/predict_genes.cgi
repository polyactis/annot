#!/usr/bin/env python

import cgi, os, time, httplib, cStringIO
import cgitb

cgitb.enable()
print "Content-Type: text/html\n\n"
form=cgi.FieldStorage()

genes=form['gene'].value
fileitem=form['genefile']
p_value = float(form['p_value'].value)
recurrence = int(form['recurrence'].value)
connectivity = float(form['connectivity'].value)
email = form['email'].value

#the 'genes' area is blank, get it from 'genefile'
if not genes:
	genes=fileitem.file.read()
#'genefile' is empty
if not genes:
	print '<h1 align=center>Where are the genes?</h1>'
else:
	if not os.path.isdir('tmp'):
				os.makedirs('tmp')
	infname = 'tmp/gene_%f'%(time.time())
	outfname = 'tmp/tab_%f'%(time.time())
	#write the genes into the input file
	inf = open(infname, 'w')
	inf.write(genes)
	inf.close()
	print "<pre>"
	print "Output file is at http://hto-pc44.usc.edu/%s.\n"%(outfname)
	print "It usually takes more than half an hour to finish.\n"
	print "If you've given an email address, the result will be emailed to your mailbox.\n"
	print "Thanks.\n"
	wl = ['predict_genes.py', '-p', repr(p_value), '-n', repr(connectivity), '-y', repr(recurrence), '-e', email, infname, outfname]
	pid = os.spawnvp(os.P_NOWAIT, '/home/yh/script/annot/bin/predict_genes.py', wl)
	print pid

	print "</pre>"
