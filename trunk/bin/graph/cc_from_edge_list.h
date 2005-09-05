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

using namespace boost;
using namespace boost::python;
using namespace std;

typedef subgraph<adjacency_list<vecS, vecS, undirectedS,
	property<vertex_name_t, int>, property<edge_index_t, int> > > Graph;	//watch subgraph<>
typedef graph_traits < Graph >::vertex_descriptor vertexDescriptor;
typedef property_map<Graph, vertex_name_t>::type vertex_name_type;

class cc_from_edge_list
{
	public:
	void init_graph_from_edge_list(boost::python::list edge_list, Graph &graph);
	std::vector<Graph> cc2subgraph(Graph &graph);
	boost::python::list subgraph2list(Graph &subgraph, Graph &graph);
	void run(boost::python::list edge_list);
	
	//for input purpose
	std::map<int, vertexDescriptor> geneNoMap;
	Graph g;
	std::vector<Graph> vector_subgraph;
	boost::python::list cc_list;
};
