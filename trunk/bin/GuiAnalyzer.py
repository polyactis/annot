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
from cluster_info import cluster_info

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
		self.entry_p_gene_table = xml.get_widget("entry_p_gene_table")
		self.entry_gene_p_table = xml.get_widget("entry_gene_p_table")
		self.button_db_connect = xml.get_widget("button_db_connect")
		
		self.window_cluster_info1 = xml.get_widget("window_cluster_info1")
		self.window_cluster_info1.hide()
		self.textview_subgraph = xml.get_widget("textview_subgraph")
		self.treeview_dataset = xml.get_widget("treeview_dataset")
		self.button_dataset_plot = xml.get_widget("button_dataset_plot")
		self.treeview_go_association = xml.get_widget("treeview_go_association")
		self.button_go_plot = xml.get_widget("button_go_plot")
		self.entry_cluster_id = xml.get_widget("entry_cluster_id")
		self.button_search_cluster_id = xml.get_widget("button_search_cluster_id")
		
		self.window_cluster_info2 = xml.get_widget("window_cluster_info2")
		self.window_cluster_info2.hide()
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
	
		self.window_cluster_accuracy = xml.get_widget("window_cluster_accuracy")
		#self.window_cluster_accuracy.hide()
		self.textview_cluster_id = xml.get_widget("textview_cluster_id")
		self.button_cluster_id_input = xml.get_widget("button_cluster_id_input")
		self.treeview_cluster_accuracy = xml.get_widget("treeview_cluster_accuracy")

		#borrowed class instances
		self.cluster_info_instance = cluster_info()


	def yeast_dataset_liststore(self):
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
	
	def on_button_cluster_info_clicked(self, button_cluster_info, *args):
		self.treeview_init(self.no_of_datasets)
		self.window_cluster_info1.show()
		self.window_cluster_info2.show()
	
	def on_button_cluster_accuracy_clicked(self, button_cluster_accuracy, *args):
		self.window_cluster_accuracy.show()
	
	def on_button_db_connect_clicked(self, button_db_connect, *args):
		hostname = self.entry_hostname.get_text()
		dbname = self.entry_dbname.get_text()
		schema = self.entry_schema.get_text()
		self.no_of_datasets = int(self.entry_no_of_datasets.get_text())
		self.conn, self.curs = db_connect(hostname, dbname, schema)
	
	def on_button_search_cluster_id_clicked(self, button_search_cluster_id, *args):		
		cluster_id = self.entry_cluster_id.get_text()
		splat_table = self.entry_splat_table.get_text()
		mcl_table = self.entry_mcl_table.get_text()
		cluster = self.cluster_info_instance.get_cluster_dstructure(self.curs, cluster_id, splat_table, mcl_table)
		
		if cluster==None:
			return

		#deal with go information
		liststore_go_association = gtk.ListStore(int, str, str, int, int, int, float)	#7 entries
		fill_treeview(self.treeview_go_association, liststore_go_association,\
			cluster.go_no2information.values())
		
		liststore_recurrence = self.yeast_dataset_liststore()
		fill_treeview(self.treeview_recurrence, liststore_recurrence, \
			[cluster.recurrence_array])
		
		self.label_total_recurrence.set_text("Total recurrence: %f"%(sum(cluster.recurrence_array)) )
		self.label_total_edges.set_text("Total no of edges: %d"%(len(cluster.edge_set)) )
		self.label_avg_connectivity.set_text("Average connectivity: %f"%(cluster.connectivity) )
		self.label_splat_connectivity.set_text("Splat connectivity: %f"%(cluster.splat_connectivity) )
		self.label_original_connectivity.set_text("Original connectivity: %f"%(cluster.connectivity_original) )
		
		
		liststore_edge_correlation = self.yeast_dataset_liststore()
		fill_treeview(self.treeview_edge_correlation, liststore_edge_correlation,\
			cluster.edge_cor_2d_list)
		liststore_edge_significance = self.yeast_dataset_liststore()
		fill_treeview(self.treeview_edge_significance, liststore_edge_significance,\
			cluster.edge_sig_2d_list)
		
	
	def treeview_init(self, no_of_datasets):
		"""
		04-18-05
			initialize all treeview column headers.
		"""
		go_association_label_list = ['go-no', 'go-id', 'name', 'depth', 'population size', 'local size', 'p-value']
		create_columns(self.treeview_go_association, go_association_label_list)
		
		dataset_label_list = map(str,range(1,no_of_datasets+1))
		create_columns(self.treeview_edge_correlation, dataset_label_list)
		create_columns(self.treeview_edge_significance, dataset_label_list)
		
		create_columns(self.treeview_recurrence, dataset_label_list)
		create_columns(self.treeview_connectivity, dataset_label_list)
		
		
	
	def destroy(self, widget):
		gtk.main_quit()
	

instance = GuiAnalyzer()
gtk.main()
