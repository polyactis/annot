#!/usr/bin/env python
"""
Gui program to analyze clusters.
"""
import os, sys, pygtk, Image
pygtk.require('2.0')

import gtk, gtk.glade
from scipy import fromimage
from codense.common import db_connect
from codense.common import foreach_cb, create_columns, fill_treeview
from codense.common import get_gene_no2gene_id, get_gene_no2go_no
from cluster_info import cluster_info
from codense.common import system_call, graphDotOutput

class GuiAnalyzer:
	"""
	04-18-05
		The integrated program to analyze clusters
	"""
	def __init__(self):
		xml = gtk.glade.XML('guianalyzer.glade')
		xml.signal_autoconnect(self)
		self.window_center = xml.get_widget("window_center")
		self.window_center.connect("destroy", self.destroy)
		self.button_cluster_info = xml.get_widget("button_cluster_info")
		self.button_cluster_accuracy = xml.get_widget("button_cluster_accuracy")
		self.entry_hostname = xml.get_widget("entry_hostname")
		self.entry_dbname = xml.get_widget("entry_dbname")
		self.entry_schema = xml.get_widget("entry_schema")
		self.entry_no_of_datasets = xml.get_widget("entry_no_of_datasets")
		self.entry_splat_table = xml.get_widget("entry_splat_table")
		self.entry_mcl_table = xml.get_widget("entry_mcl_table")
		self.entry_cluster_stat = xml.get_widget("entry_cluster_stat")
		self.togglebutton_hide_show = xml.get_widget("togglebutton_hide_show")

		self.button_db_connect = xml.get_widget("button_db_connect")
		
		self.window_cluster_info1 = xml.get_widget("window_cluster_info1")
		self.window_cluster_info1.connect("delete_event", self.subwindow_hide)
		self.window_cluster_info1.show()
		self.textview_subgraph = xml.get_widget("textview_subgraph")
		self.treeview_dataset = xml.get_widget("treeview_dataset")
		self.button_dataset_plot = xml.get_widget("button_dataset_plot")
		self.treeview_go_association = xml.get_widget("treeview_go_association")
		self.button_go_plot = xml.get_widget("button_go_plot")
		self.entry_cluster_id = xml.get_widget("entry_cluster_id")
		self.button_search_cluster_id = xml.get_widget("button_search_cluster_id")
		
		self.window_cluster_info2 = xml.get_widget("window_cluster_info2")
		self.window_cluster_info2.connect("delete_event", self.subwindow_hide)
		self.window_cluster_info2.show()
		self.treeview_recurrence = xml.get_widget("treeview_recurrence")
		self.label_total_recurrence = xml.get_widget("label_total_recurrence")
		self.label_total_edges = xml.get_widget("label_total_edges")
		self.treeview_connectivity = xml.get_widget("treeview_connectivity")
		self.label_avg_connectivity = xml.get_widget("label_avg_connectivity")
		self.label_splat_connectivity = xml.get_widget("label_splat_connectivity")
		self.label_original_connectivity = xml.get_widget("label_original_connectivity")
		self.treeview_edge_correlation = xml.get_widget("treeview_edge_correlation")
		self.button_edge_cor_plot = xml.get_widget("button_edge_cor_plot")
		self.treeview_edge_significance = xml.get_widget("treeview_edge_significance")
		self.button_edge_sig_plot = xml.get_widget("button_edge_sig_plot")
	
		self.dialog_cluster_accuracy = xml.get_widget("dialog_cluster_accuracy")
		self.dialog_cluster_accuracy.connect("delete_event", self.subwindow_hide)
		self.dialog_cluster_accuracy.show()
		self.textview_cluster_id = xml.get_widget("textview_cluster_id")
		self.okbutton_cluster_accuracy = xml.get_widget("okbutton_cluster_accuracy")
		self.treeview_cluster_accuracy = xml.get_widget("treeview_cluster_accuracy")
		self.entry_cluster_accuracy_p_gene_table = xml.get_widget("entry_cluster_accuracy_p_gene_table")
		self.entry_cluster_accuracy_p_value_cut_off = xml.get_widget("entry_cluster_accuracy_p_value_cut_off")
		self.comboboxentry_plot_type = xml.get_widget("comboboxentry_plot_type")	#06-14-05	to offer selection of plot type.
		self.comboboxentry_plot_type.set_active(0)	#06-14-05 set default
		
		#borrowed class instances
		self.cluster_info_instance = cluster_info()
		
		self.cluster_info_need_init = 1	#flag whether the windows of cluster_info need to be initialized
		self.curs = None	#used to check in on_button_cluster_info_clicked()
		self.cluster = None	#used to check in on_button_go_plot_clicked()
		self.dataset_liststore_dict = {54:self.dataset_liststore_54,
			55:self.dataset_liststore_55,
			73:self.dataset_liststore_73,
			74:self.dataset_liststore_74,
			79:self.dataset_liststore_79,
			80:self.dataset_liststore_80,
			100:self.dataset_liststore_100,
			101:self.dataset_liststore_101}
		
		
	def dataset_liststore_54(self):
		"""
		04-18-05
			return a 54 column liststore
		"""
		return gtk.ListStore(\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str)	#54 entries	
	
	def dataset_liststore_55(self):
		"""
		04-19-05
			return a 55 column liststore
		"""
		return gtk.ListStore(\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str, str)	#55 entries
	
	
	def dataset_liststore_73(self):
		"""
		05-19-05
			return a 73 column liststore
		"""
		return gtk.ListStore(\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str)	
	
	
	def dataset_liststore_74(self):
		"""
		05-19-05
			return a 74 column liststore
		"""
		return gtk.ListStore(\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str)	
	
	def dataset_liststore_79(self):
		"""
		04-19-05
			return a 79 column liststore
		"""
		return gtk.ListStore(\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str)	#79 entries
	
	def dataset_liststore_80(self):
		"""
		04-19-05
			return a 80 column liststore
		"""
		return gtk.ListStore(\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str)	#80 entries
	
	def dataset_liststore_100(self):
		"""
		08-04-05
		"""
		return gtk.ListStore(\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str)	#100 entries

	def dataset_liststore_101(self):
		"""
		08-04-05
		"""
		return gtk.ListStore(\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str,str,str,str,str, str,str,str,str,str,\
			str)	#101 entries

	"""
	window_center
	"""
	def on_button_cluster_info_clicked(self, button_cluster_info, *args):
		if self.curs==None:
			print "db_connect first"
			return
		self.no_of_datasets = int(self.entry_no_of_datasets.get_text())
		self.treeview_init(self.no_of_datasets)
		self.gene_no2gene_id = get_gene_no2gene_id(self.curs)
		self.gene_no2go_no = get_gene_no2go_no(self.curs)

		self.window_cluster_info1.show()
		self.window_cluster_info2.show()
	
	def on_button_cluster_accuracy_clicked(self, button_cluster_accuracy, *args):
		self.dialog_cluster_accuracy.show()
	
	def on_button_db_connect_clicked(self, button_db_connect, *args):
		hostname = self.entry_hostname.get_text()
		dbname = self.entry_dbname.get_text()
		schema = self.entry_schema.get_text()
		self.conn, self.curs = db_connect(hostname, dbname, schema)
		
	def on_togglebutton_hide_show_toggled(self, widget, data=None):
		"""
		04-24-05
			a toggle button to hide or show subwindows
		"""
		if widget.get_active():
			self.window_cluster_info1.hide()
			self.window_cluster_info2.hide()
			self.dialog_cluster_accuracy.hide()
			widget.set_label("Subwindow Show")
		else:
			self.window_cluster_info1.show()
			self.window_cluster_info2.show()
			self.dialog_cluster_accuracy.show()
			widget.set_label("Subwindow Hide")
	"""
	cluster_info part
	"""
	
	def on_button_search_cluster_id_clicked(self, button_search_cluster_id, *args):	
		"""
		04-19-05
			add edge_name
		"""
		cluster_id = self.entry_cluster_id.get_text()
		splat_table = self.entry_splat_table.get_text()
		mcl_table = self.entry_mcl_table.get_text()
		self.cluster = self.cluster_info_instance.get_cluster_dstructure(self.curs, cluster_id, splat_table, mcl_table)
		
		if self.cluster==None:
			return

		#deal with go information
		self.liststore_go_association = gtk.ListStore(int, str, str, int, int, int, float)	#7 entries
		fill_treeview(self.treeview_go_association, self.liststore_go_association,\
			self.cluster.go_no2information.values())
		
		liststore_recurrence = self.dataset_liststore_dict[self.no_of_datasets]()
		fill_treeview(self.treeview_recurrence, liststore_recurrence, \
			[self.cluster.recurrence_array])
		
		self.label_total_recurrence.set_text("Total recurrence: %s"%(sum(self.cluster.recurrence_array)) )
		self.label_total_edges.set_text("Total no of edges: %s"%(len(self.cluster.edge_set)) )
		self.label_avg_connectivity.set_text("Average connectivity: %s"%(self.cluster.connectivity) )
		self.label_splat_connectivity.set_text("Splat connectivity: %s"%(self.cluster.splat_connectivity) )
		self.label_original_connectivity.set_text("Original connectivity: %s"%(self.cluster.connectivity_original) )
		
		
		liststore_edge_correlation = self.dataset_liststore_dict[self.no_of_datasets+1]()	#it's one column bigger than before
		name2edge_cor_2d_list = self.add_name2edge_2d_list(self.cluster.edge_set, \
			self.cluster.edge_cor_2d_list, self.gene_no2gene_id)
		fill_treeview(self.treeview_edge_correlation, liststore_edge_correlation,\
			name2edge_cor_2d_list)
		liststore_edge_significance = self.dataset_liststore_dict[self.no_of_datasets+1]()
		name2edge_sig_2d_list = self.add_name2edge_2d_list(self.cluster.edge_set,\
			self.cluster.edge_sig_2d_list, self.gene_no2gene_id)
		fill_treeview(self.treeview_edge_significance, liststore_edge_significance,\
			name2edge_sig_2d_list)
		
	def add_name2edge_2d_list(self, edge_set, edge_2d_list, label_dict):
		"""
		04-19-05
			add name ahead each list
		"""
		name2edge_2d_list = []
		for i in range(len(edge_set)):
			edge = edge_set[i]
			edge_name = '%s:%s'%(label_dict[edge[0]], label_dict[edge[1]])
			name2edge_2d_list.append([edge_name]+edge_2d_list[i])
		return name2edge_2d_list
	
	def treeview_init(self, no_of_datasets):
		"""
		04-18-05
			initialize all treeview column headers.
		"""
		go_association_label_list = ['go-no', 'go-id', 'name', 'depth', 'population size', 'local size', 'p-value']
		create_columns(self.treeview_go_association, go_association_label_list)
		
		dataset_label_list = map(str,range(1,no_of_datasets+1))
		create_columns(self.treeview_recurrence, dataset_label_list)
		create_columns(self.treeview_connectivity, dataset_label_list)
		#'edge' put first
		create_columns(self.treeview_edge_correlation, ['edge'] + dataset_label_list)
		create_columns(self.treeview_edge_significance, ['edge'] + dataset_label_list)

	
	def on_button_go_plot_clicked(self, widget):
		"""
		04-19-05
		
		06-14-05
			use graphDotOutput() and insert it into the textbuffer
		10-11-05
			add cluster_id and go_no for each graph
		"""
		if self.cluster == None:
			print "Cluster not loaded in"
			return
		
		#first construct the graph
		subgraph = self.cluster_info_instance.graph_from_node_edge_set(\
			self.cluster.vertex_set, self.cluster.edge_set)
			
		pathlist = []
		treeselection_go_association = self.treeview_go_association.get_selection()
		treeselection_go_association.selected_foreach(foreach_cb, pathlist)
		if len(pathlist) >0:
			for i in range(len(pathlist)):
				go_no = self.liststore_go_association[pathlist[i][0]][0]
				graphSrcFname = '/tmp/GuiAnalyzer.dot'
				graphFname = '/tmp/GuiAnalyzer.ps'
				graphSrcF = open(graphSrcFname, 'w')
				graphDotOutput(graphSrcF, subgraph, \
					self.gene_no2gene_id, self.gene_no2go_no, \
					function=int(go_no), weighted=0)
				graphSrcF.close()
				#06-14-05 get the plot_type_command from the comboboxentry_plot_type
				entry = self.comboboxentry_plot_type.child
				plot_type_command = entry.get_text()
				if plot_type_command=='neato' or plot_type_command=='twopi':
					plot_type_command += ' -Goverlap=false' 	#06-14-05	-Goverlap=false passed to neato or twopi
				commandline = '%s -Tps %s -o %s'%(plot_type_command, graphSrcFname, graphFname)
				
				system_call(commandline)
				
				#10-11-05 temporary adding this block to draw pictures in png format
				outputfname = '/tmp/%s_%s.png'%(self.cluster.cluster_id, go_no)
				commandline = '%s -Tpng %s -o %s'%(plot_type_command, graphSrcFname, outputfname)
				system_call(commandline)
				
				textbuffer = self.textview_subgraph.get_buffer()
				im = Image.open(graphFname)
				ar = fromimage(im)
				pixbuf = gtk.gdk.pixbuf_new_from_array(ar, gtk.gdk.COLORSPACE_RGB, 8)
				startiter = textbuffer.get_start_iter()
				textbuffer.insert_pixbuf(startiter, pixbuf)
				startiter = textbuffer.get_start_iter()
				textbuffer.insert(startiter, "cluster_id:%s go_no:%s"%(self.cluster.cluster_id, go_no))
				#cleanup the temporary files and others
				del im
				del ar
	
	"""
	cluster_accuracy part
	"""
	
	def on_cancelbutton_cluster_accuracy_clicked(self, widget):
		self.dialog_cluster_accuracy.hide()
		
	def subwindow_hide(self, widget, event, data=None):
		widget.hide()
		return True
	
	def hide_all_subwindows(self):
		"""
		04-20-05
			called after the class's initialization to hide subwindows.
		"""
		self.window_cluster_info1.hide()
		self.window_cluster_info2.hide()
		self.dialog_cluster_accuracy.hide()
	
	def destroy(self, widget):
		gtk.main_quit()
	

instance = GuiAnalyzer()
gtk.main()
