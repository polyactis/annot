#!/usr/bin/env python

import cgi, os, time, httplib, cStringIO
import cgitb
from GO_common_ancestor import GO_common_ancestor
from GO_graphxml import GO_graphxml

cgitb.enable()
print "Content-Type: text/html\n\n"
form=cgi.FieldStorage()

if form.has_key('go_acc_list'):
	go_acc_list=form['go_acc_list'].value
	if not os.path.isdir('tmp'):
		os.makedirs('tmp')
	inf = open('tmp/go_acc_list', 'w')
	inf.write(go_acc_list)
	inf.close()
	instance = GO_common_ancestor('zhoudb', 'graphdb', 'go', 1, '/var/www/repositary/hypergraph/go.xml', 'sc', 'tmp/go_acc_list')
	instance.run()
	
elif form.has_key('go_acc'):
	go_acc = form['go_acc'].value
	backward = 0
	if form.has_key('backward'):
		backward = form['backward'].value
		backward = 1
	instance = GO_graphxml('zhoudb', 'graphdb', 'go', 1, '/var/www/repositary/hypergraph/go.xml', 'sc')
	instance.cgi_run(go_acc, backward)
