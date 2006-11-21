#!/usr/bin/env python

go_node_size_cutoff = range(100, 251, 25)
no_of_categories_ls = [56, 49, 42, 40, 33, 33, 30]
accuracy_ls = [78.8, 77.1, 78.6, 77.3, 79.1, 76.4, 77.9]
no_of_unknowns_ls = [70, 90, 82, 124, 115, 147, 141]
no_of_knowns_ls = [413, 494, 581, 745, 785, 881, 899]

import matplotlib; matplotlib.use("Agg")
import pylab

def plot_func(x_ls, y_ls, xlab, ylab, output_fname):
	pylab.xlabel(xlab)
	pylab.ylabel(ylab)
	pylab.plot(x_ls, y_ls, marker='o')
	pylab.savefig(output_fname)
	pylab.clf()

plot_func(go_node_size_cutoff, no_of_categories_ls, 'GO node size cutoff','No of categories', 'go_node_size_cutoff_vs_no_of_categories.png' )
plot_func(go_node_size_cutoff, accuracy_ls, 'GO node size cutoff','prediction accuracy', 'go_node_size_cutoff_vs_accuracy.png' )
plot_func(go_node_size_cutoff, no_of_unknowns_ls, 'GO node size cutoff','No of unknowns', 'go_node_size_cutoff_vs_no_of_unknowns.png' )
plot_func(go_node_size_cutoff, no_of_knowns_ls, 'GO node size cutoff','No of knowns', 'go_node_size_cutoff_vs_no_of_knowns.png' )


