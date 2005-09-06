/*
*04-29-05
*	a module to extract connected components(cc) from an edge list.
*
*/

#include <boost/graph/adjacency_list.hpp>	//for graph definition
#include <boost/graph/connected_components.hpp>	//for connected_components
#include <boost/graph/subgraph.hpp>	//for boost::subgraph
#include <boost/graph/graph_traits.hpp>  // for boost::graph_traits
#include <boost/utility.hpp>             // for boost::tie
#include <boost/python.hpp>	//for python module and dict, tuple, make_tuple

#include <iostream>			// for std::cout
#include <vector>			//for vector
#include <boost/property_map.hpp>	//09-05-05 for iterator_property_map
#include <boost/graph/betweenness_centrality.hpp>	//09-05-05 for betweenness_centrality
#include <boost/pending/indirect_cmp.hpp>	//Mon Sep  5 11:38:25 2005 for indirect_cmp
#include <boost/tuple/tuple.hpp>	//Mon Sep  5 19:03:02 2005 for boost::tuple

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
typedef	boost::iterator_property_map<vector<double>::iterator,\
	EdgeIndexMap, double, double&> EdgeCentralityMap;
typedef property_traits<EdgeCentralityMap>::value_type centrality_type;
//Mon Sep  5 18:45:37 2005
typedef graph_traits<Graph>::vertex_iterator vertexIterator;
typedef boost::tuple<boost::python::list, boost::python::list> twoListTuple;

class cc_from_edge_list
{
	public:
	void init_graph_from_edge_list(boost::python::list edge_list, Graph &graph);
	std::vector<Graph> cc2subgraph(Graph &graph);
	boost::python::list subgraph2list(Graph &subgraph, Graph &graph);
	std::vector<Graph> subgraph_components(Graph &subgraph, Graph &graph, std::vector<int> component, int no_of_components);
	std::vector<Graph> graph_components(Graph &graph, std::vector<int> component, int no_of_components);
	void run(boost::python::list edge_list);
	
	//for input purpose
	std::map<int, vertexDescriptor> geneNoMap;
	Graph g;
	std::vector<Graph> vector_subgraph;
	boost::python::list cc_list;
};

class ClusterByEBC:public cc_from_edge_list
{
	public:
		void cut_by_betweenness_centrality(Graph &graph, int size_cutoff, float conn_cutoff);
		void run(boost::python::list edge_list, int size_cutoff, float conn_cutoff);
		twoListTuple graph2list(Graph &graph);
		std::vector<Graph> good_subgraph_vector;
		boost::python::list cc_vertex_list;	//09-05-05	to return the vertex_list as well
};
