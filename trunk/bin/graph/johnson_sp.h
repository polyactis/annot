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
//01-23-06
#include <ext/hash_map>	//for hash_map
#include <utility>	///for pair
#define hash_edge_name(o1, o2) ((o1<<30) + o2)
#include <boost/dynamic_bitset.hpp> //for dynamic_bitset

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
	johnson_sp();	//01-23-06
	johnson_sp(int no_of_datasets);	//01-23-06
	void add_edge_sig_vector(boost::python::list edge_sig_list);	//01-23-06
	
	std::pair<Graph, std::vector<boost::dynamic_bitset<> > > init_graph_from_vertex_edge_list(boost::python::list vertex_list, boost::python::list edge_list);
	boost::python::list calculate_sp(Graph &graph);
	boost::python::list py_shortest_distance(boost::python::list vertex_list, boost::python::list edge_list);
	boost::python::list py_recurrence_list();
	
	__gnu_cxx::hash_map<unsigned long, boost::dynamic_bitset<> > edge2bitset;
	const int _no_of_datasets;
	std::vector<boost::dynamic_bitset<> > _recurrence_bitset_vector;
};
