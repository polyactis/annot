/*
*10-09-05
*	a module to calculate the shortest path for a graph.
*
*/

#include <boost/graph/adjacency_list.hpp>	//for graph definition
#include <boost/graph/subgraph.hpp>	//for boost::subgraph
#include <boost/graph/graph_traits.hpp>  // for boost::graph_traits
#include <boost/graph/johnson_all_pairs_shortest.hpp>	//for boost::johnson_all_pairs_shortest_paths()

#include <boost/python.hpp>	//for python module and dict, tuple, make_tuple

#include <iostream>			// for std::cout
#include <vector>			//for vector

using namespace boost;
using namespace boost::python;
using namespace std;

typedef adjacency_list<vecS, vecS, undirectedS,
	property<vertex_name_t, int>, property< edge_weight_t, int, property< edge_weight2_t, int > > >  Graph;
typedef graph_traits < Graph >::vertex_descriptor vertexDescriptor;
typedef graph_traits<Graph>::edge_descriptor edgeDescriptor;

class johnson_sp
{
	public:
	void init_graph_from_vertex_edge_list(boost::python::list vertex_list, boost::python::list edge_list, Graph &graph);
	
	void calculate_sp(Graph &graph, boost::python::list &D_matrix);
	void run(boost::python::list vertex_list, boost::python::list edge_list);
	
	//for input purpose
	std::map<int, vertexDescriptor> geneNoMap;
	Graph g;
	boost::python::list D;
};
