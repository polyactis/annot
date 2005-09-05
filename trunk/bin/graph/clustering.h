/*
*
*
*02-24-05
*	boost and boost::python conflict for the make_tuple function. use boost::python::make_tuple instead.*
*05-25-05
*	modify it to be a standalone program and do normalized cut in a way similar to netmine(copath).
*	cut the graph until its size is less than a threshold.
*/
#include <boost/config.hpp>
#include <iostream>                      // for std::cout
#include <utility>                       // for std::pair
#include <algorithm>		//for std::sort and std::copy	05-30-05
#include "areig.h"	//for AREig()	05-30-05
#include <boost/utility.hpp>             // for boost::tie
#include <boost/graph/subgraph.hpp>	//for boost::subgraph
#include <boost/graph/graph_traits.hpp>  // for boost::graph_traits
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/connected_components.hpp>	//for connected_components

//for parsing input file
#include <fstream>
#include <string>
#include <boost/tokenizer.hpp>	//for tokenizer, parse input file
#include <boost/tuple/tuple.hpp>

#include <gsl/gsl_math.h>		//for gsl_matrix, gsl_vector stuff
#include <gsl/gsl_eigen.h>		//for gsl_eigen stuff

#include <boost/python.hpp>		//for python module and dict, tuple, make_tuple

#include <getopt.h>	//05-25-05	to parse program options.

using namespace boost;

using namespace boost::python;

typedef subgraph<adjacency_list<vecS, vecS, undirectedS,
	property<vertex_name_t, int>, property<edge_index_t, int> > > Graph;
typedef graph_traits < Graph >::vertex_descriptor vertexDescriptor;
typedef property_map<Graph, vertex_name_t>::type vertexNamePropertyMap;
typedef property_map<Graph, boost::edge_index_t>::type EdgeIndexMap;

class clustering
{
	public:
		clustering();
		clustering(double connectivity, int which_eigen_vector, int cluster_size);
		clustering(std::string input_filename_option, std::string output_filename_option, int max_size_option,
			int cut_loop_num_option, double density_cutoff_option, int min_edge_weight_option, int matrix_format_option);	//05-25-05
		~clustering();
		void init_graph_from_dict(dict graph_dict, Graph &graph);
		void init_graph_from_file(std::string input_filename, Graph &graph, int min_edge_weight);	//05-25-05
		void init_graph_from_file_matrix(std::string input_filename, Graph &graph, int min_edge_weight);	//05-26-05
		dict graph2dict(Graph &subgraph, Graph &graph);
		gsl_matrix *graph2gsl_matrix(Graph &graph);
		gsl_vector *return_eigen_vector(gsl_matrix* graph_matrix, int which_eigen_vector);
		gsl_vector *return_eigen_vector_by_arpack(Graph &graph, int which_eigen_vector);
		double connectivity_of_graph(Graph &graph);
		std::vector<Graph> subgraph_components(Graph &subgraph, Graph &graph, std::vector<int> component, int no_of_components);
		void normalized_cut(std::ofstream &outf, Graph &subgraph, Graph &graph, int max_size, int eigen_vector_no);	//05-25-05
		void old_run(dict graph_dict);
		void run();	//05-25-05
		void walk_graph(std::ofstream &outf, Graph &subgraph, Graph &graph);	//05-25-05
		void output(std::string output_filename);	//05-25-05
		
		Graph g;
		//for input purpose
		std::map<int, vertexDescriptor> geneNoMap;
		//property_map of the vertex names
		std::map<std::string, vertexDescriptor> VertexNameMap;

		vertexNamePropertyMap vertex2name;
		int min_cluster_size;
		int eigen_vector_no;
		double connectivity_cutoff;
		std::string input_filename;
		std::string output_filename;
		int max_size;
		int cut_loop_num;
		int min_edge_weight;
		int matrix_format;
		//the weight map
		iterator_property_map<int*, EdgeIndexMap, int, int&> weight_pa;
		//store the final results
		boost::python::list good_clusters;
		//store the final results
		std::vector<Graph> good_clusters_vector;

};
