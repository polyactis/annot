#!/usr/bin/env python

import os,sys
from Tkinter import *
from GenePair import GenePair

class GenePairGui:
	def __init__(self, master, dir):
		self.dir = dir
		frame = Frame(master)
		frame.pack()
		self.scrollbar_file = Scrollbar(frame, orient=VERTICAL)
		self.scrollbar_gene_pair = Scrollbar(frame, orient=VERTICAL)
		
		self.file_listbox = Listbox(frame, yscrollcommand=self.scrollbar_file.set)
		self.gene_pair_listbox = Listbox(frame, selectmode=EXTENDED, yscrollcommand=self.scrollbar_gene_pair.set)
		
		self.button_file = Button(frame, text="select file", command=self.file_selected)
		self.button_gene_pair = Button(frame, text="select gene pair", command=self.gene_pair_selected)
		self.button_quit = Button(frame, text="QUIT", fg="red", command=frame.quit)
		
		#geometry arrangement
		self.scrollbar_file.config(command=self.file_listbox.yview)
		self.scrollbar_gene_pair.config(command=self.gene_pair_listbox.yview)
		self.scrollbar_file.pack(fill=Y, expand=1)
		self.scrollbar_gene_pair.pack(fill=Y, expand=1)
		
		self.file_listbox.grid(row=0, column=0)
		self.scrollbar_file.grid(row=0,column=1)
		self.gene_pair_listbox.grid(row=0, column=2)
		self.scrollbar_gene_pair.grid(row=0, column=3)
		
		self.button_file.grid(row=1, column=0, columnspan=2)
		self.button_gene_pair.grid(row=1, column=2, columnspan=2)
		self.button_quit.grid(row=2, column=1, columnspan=2)
	
		#fill the list box
		self.file_listbox_fill()
		#GenePair handler
		self.GenePair_instance = None
		
	def file_listbox_fill(self):
		files = os.listdir(self.dir)
		for item in files:
			self.file_listbox.insert(END, item)
		
	def file_selected(self):
		file = self.file_listbox.curselection()
		if len(file) != 1:
			sys.stderr.write("No file selected.\n")
		else:
			file = map(int, file)
			fname = self.file_listbox.get(file[0])
			file_path = os.path.join(self.dir, fname)
			if self.GenePair_instance:
				#clean it up
				self.GenePair_instance = None
			else:
				self.GenePair_instance = GenePair(file_path)
				gene_ids = self.GenePair_instance.gene_id2expr_array.keys()
				gene_ids.sort()
				for gene_id in gene_ids:
					self.gene_pair_listbox.insert(END, gene_id)

	def gene_pair_selected(self):
		gene_pair = self.gene_pair_listbox.curselection()
		if len(gene_pair) !=2:
			sys.stderr.write("Did you select two genes?\n")
		else:
			gene_pair = map(int, gene_pair)
			gene_id1 = self.gene_pair_listbox.get(gene_pair[0])
			gene_id2 = self.gene_pair_listbox.get(gene_pair[1])
			print gene_id1
			print gene_id2
			self.GenePair_instance.gene_pair_analyze(gene_id1, gene_id2)

if __name__ == '__main__':
	if len(sys.argv) == 1:
		print __doc__
		sys.exit(2)
	else:
		root = Tk()
		app = GenePairGui(root, sys.argv[1])
		root.mainloop()
