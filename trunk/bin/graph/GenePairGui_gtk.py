#!/usr/bin/env python

import os, sys, pygtk, Image
pygtk.require('2.0')

import gtk, gtk.glade
from GenePair import GenePair
from scipy import fromimage
"""
02-19-05 
	use fromimage to convert a PIL's Image to a numeric array and 
	pass the array to gtk.gdk.pixbuf_new_from_array.
"""
def foreach_cb(model, path, iter, pathlist):
	pathlist.append(path)	


class GenePairGui:
	def __init__(self):
		xml = gtk.glade.XML('genepairgui.glade')
		xml.signal_autoconnect(self)
		self.window1 = xml.get_widget("window1")
		self.window1.connect("destroy", self.destroy)
		#to get gene_ids input by user
		self.textview1 = xml.get_widget("textview1")
		#to draw the plot
		self.textview2 = xml.get_widget("textview2")
		#to get gene_ids input by user
		self.treeview1 = xml.get_widget("treeview1")
		# create the TreeViewColumn to display the data
		self.tvcolumn = gtk.TreeViewColumn('Gene Id')
		# add tvcolumn to treeview
		self.treeview1.append_column(self.tvcolumn)
		# create a CellRendererText to render the data
		self.cell = gtk.CellRendererText()
		# add the cell to the tvcolumn and allow it to expand
		self.tvcolumn.pack_start(self.cell, True)
		
		# set the cell "text" attribute to column 0 - retrieve text
		# from that column in liststore
		self.tvcolumn.add_attribute(self.cell, 'text', 0)
		
		#setting the selection mode
		self.treeselection = self.treeview1.get_selection()
		
		#the backend module, non-initilized
		self.GenePair_instance = None
		
	def file_ok_sel(self, w):
		filename = self.filew.get_filename()
		self.filew.destroy()
		self.window1.set_title(filename)
		self.liststore_gene_fill(filename)
	
	def liststore_gene_fill(self, filename):
		self.GenePair_instance = GenePair(filename)
		gene_ids = self.GenePair_instance.gene_id2expr_array.keys()
		gene_ids.sort()

		self.liststore = gtk.ListStore(str)
		# we'll add some data now -
		for gene_id in gene_ids:
			self.liststore.append([gene_id])
		
		# set the TreeView mode to be liststore
		self.treeview1.set_model(self.liststore)

		# make it searchable
		self.treeview1.set_search_column(0)
		# Allow sorting on the column
		self.tvcolumn.set_sort_column_id(0)
		# Allow drag and drop reordering of rows
		self.treeview1.set_reorderable(True)

		self.treeselection.set_mode(gtk.SELECTION_MULTIPLE)

	def on_button1_clicked(self, button):
		self.filew = gtk.FileSelection("File selection")
		#self.filew.connect("destroy", self.destroy)
		# Connect the ok_button to file_ok_sel method
		self.filew.ok_button.connect("clicked", self.file_ok_sel)
		# Connect the cancel_button to destroy the widget
		self.filew.cancel_button.connect("clicked", lambda w: self.filew.destroy())
		# Lets set the filename, as if this were a save dialog,
		# and we are giving a default filename
		#self.filew.set_filename("penguin.png")
		
		self.filew.show()
	
	def on_button2_clicked(self, button):
		'''
		gene selected
		'''
		pathlist = []
		gene_id_list = []
		self.treeselection.selected_foreach(foreach_cb, pathlist)
		if len(pathlist) >0:
			for i in range(len(pathlist)):
				gene_id_list.append(self.liststore[pathlist[i][0]][0])
			filename = self.GenePair_instance.gene_pair_analyze(gene_id_list)
			if filename:
				self.show_plot(filename)
		else:
			sys.stderr.write("Have you selected genes?\n")

	def on_button3_clicked(self, button):
		gtk.main_quit()
		
	def on_button4_clicked(self, button):
		'''
		input gene_ids done
		'''
		textbuffer = self.textview1.get_buffer()
		startiter, enditer = textbuffer.get_bounds()
		text  = textbuffer.get_text(startiter, enditer)
		gene_id_list = text.split("\n")
		if self.GenePair_instance:
			if len(gene_id_list) >0:
				filename = self.GenePair_instance.gene_pair_analyze(gene_id_list)
				if filename:
					self.show_plot(filename)
		else:
			sys.stderr.write("Have you selected a dataset?\n")

	def show_plot(self, filename):
		textbuffer = self.textview2.get_buffer()
		im = Image.open(filename)
		#r.postscript() gives an image with a strange direction, rotate it back
		im = im.rotate(-90)
		ar = fromimage(im)
		pixbuf = gtk.gdk.pixbuf_new_from_array(ar, gtk.gdk.COLORSPACE_RGB, 8)
		startiter = textbuffer.get_start_iter()
		textbuffer.insert_pixbuf(startiter, pixbuf)
		#cleanup the temporary files and others
		del im
		del ar
		os.remove(filename)
		
	def destroy(self, widget):
		gtk.main_quit()
	

test = GenePairGui()
gtk.main()
