/*
*
*01-09-06
*      find edge items related to a fim dataset signature, a module of MpiFromDatasetSignatureToPattern.py
*/
#include <iostream>
#include <boost/dynamic_bitset.hpp>
#include <vector>
#include <boost/tokenizer.hpp>	//for tokenizer, parse input file
#include <fstream>
#include <boost/array.hpp>
#include <string>
#include <boost/python.hpp>


#include <boost/graph/adjacency_list.hpp>	//for graph definition
#include <boost/graph/connected_components.hpp>	//for connected_components
#include <boost/graph/subgraph.hpp>	//for boost::subgraph
#include <boost/graph/graph_traits.hpp>  // for boost::graph_traits
#include <boost/utility.hpp>             // for boost::tie

#include <getopt.h>	//to parse program options.

using namespace boost;
using namespace boost::python;
using namespace std;

typedef subgraph<adjacency_list<vecS, vecS, undirectedS,
	property<vertex_name_t, int>, property<edge_index_t, int> > > Graph;	//watch subgraph<>
typedef graph_traits < Graph >::vertex_descriptor vertexDescriptor;
typedef property_map<Graph, vertex_name_t>::type vertex_name_type;

//Mon Sep  5 11:38:25 2005
typedef property_map<Graph, boost::edge_index_t>::type EdgeIndexMap;
typedef graph_traits<Graph>::edge_descriptor edgeDescriptor;
typedef graph_traits<Graph>::vertex_iterator vertexIterator;

//Mon Jan  9 22:35:20 2006	for sorting in output_subgraph()
inline bool cmp_edge_array(boost::array<unsigned int, 2> e1, boost::array<unsigned int, 2> e2) { return e1[0] < e2[0]; }

class PostFim
{
	public:
		ofstream _out;		
		const int _no_cc;	// no connected components
		const int _no_of_datasets;
		const int _min_cluster_size;
		const std::string _node_outputfname;
		
		std::vector<boost::dynamic_bitset<> > pattern_bitset_vector;
		std::vector<int > freq_of_signature_vector;
		
		std::vector<unsigned int > edge_tuple_vector;
		std::vector<boost::dynamic_bitset<> > edge_bitset_vector;
		
		
		PostFim(int no_cc, int no_of_datasets, int min_cluster_size, std::string node_outputfname);
		~PostFim();
		void add_edge_sig_vector(boost::python::list edge_sig_list);
		void add_pattern_signature(boost::python::list pattern_sig_list);
		
		void patternFormation();
		void outputCcFromEdgeList(std::vector<int> &edge_id_vector, std::vector<unsigned int > &edge_tuple_vector,\
			int min_cluster_size, int no_cc);
		
		std::vector<Graph> cc2subgraph(Graph &graph);
		Graph init_graph_from_edge_tuple_vector(std::vector<int> &edge_id_vector, std::vector<unsigned int > &edge_tuple_vector);
		void output_subgraph(ofstream &out, Graph &subgraph, Graph &graph);
};
