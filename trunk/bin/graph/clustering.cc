/*
*
*
*02-24-05
*	boost and boost::python conflict for the make_tuple function. use boost::python::make_tuple instead.*
*/
#include <boost/config.hpp>
#include <iostream>                      // for std::cout
#include <utility>                       // for std::pair
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
#include <boost/tokenizer.hpp>
#include <boost/tuple/tuple.hpp>

#include <gsl/gsl_math.h>		//for gsl_matrix, gsl_vector stuff
#include <gsl/gsl_eigen.h>		//for gsl_eigen stuff

#include <boost/python.hpp>		//for python module and dict, tuple, make_tuple

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
		clustering(double connectivity, int which_eigen_vector, int cluster_size);
		~clustering();
		void init_graph_from_dict(dict graph_dict, Graph &graph);
		dict graph2dict(Graph &subgraph, Graph &graph);
		gsl_matrix *graph2gsl_matrix(Graph &graph);
		gsl_vector *return_eigen_vector(gsl_matrix* graph_matrix, int which_eigen_vector);
		double connectivity_of_graph(Graph &graph);
		std::vector<Graph> subgraph_components(Graph &subgraph, Graph &graph, std::vector<int> component, int no_of_components);
		void cluster(Graph &graph);
		void run(dict graph_dict);
		
		Graph g;
		//for input purpose
		std::map<int, vertexDescriptor> geneNoMap;
		//property_map of the vertex names
		vertexNamePropertyMap vertex2name;
		int min_cluster_size;
		int eigen_vector_no;
		double connectivity_cutoff;
		//the weight map
		iterator_property_map<int*, EdgeIndexMap, int, int&> weight_pa;
		//store the final results
		list good_clusters;

};

clustering::clustering(double connectivity, int which_eigen_vector, int cluster_size)
{
	//parameter initialization
	eigen_vector_no = which_eigen_vector;		//second minimum
	connectivity_cutoff = connectivity;	
	min_cluster_size = cluster_size;
}

clustering::~clustering()
{
	//do nothing
	std::cout<<"End"<<std::endl;
}

void clustering::init_graph_from_dict(dict graph_dict, Graph &graph)
{
	vertex2name = get(vertex_name, graph);
	
	list edge_list = graph_dict.keys();
	int edge_length = extract<int>(edge_list.attr("__len__")());
	int weight_array[edge_length];
	
	for(int i=0; i<edge_length; i++)
	{
		int gene1 = extract<int>(edge_list[i][0]);
		int gene2 = extract<int>(edge_list[i][1]);
		//used in last to get the weight.
		boost::python::tuple tup = boost::python::make_tuple(gene1, gene2);
		std::map<int, vertexDescriptor>::iterator pos;
		bool inserted;
		vertexDescriptor u, v;
		tie(pos, inserted) = geneNoMap.insert(std::make_pair(gene1, vertexDescriptor()));
		if (inserted) {
			u = add_vertex(graph);
			vertex2name[u] = gene1;
			pos->second = u;
		} else
			u = pos->second;
		
		tie(pos, inserted) = geneNoMap.insert(std::make_pair(gene2, vertexDescriptor()));
		if (inserted) {
			v = add_vertex(graph);
			vertex2name[v] = gene2;
			pos->second = v;
		} else
			v = pos->second;
		
		graph_traits < Graph >::edge_descriptor e;
		tie(e, inserted) = add_edge(u, v, graph);
		if (inserted)
			weight_array[i] =extract<int>(graph_dict[tup]);
	}
	EdgeIndexMap edge_id = get(edge_index, graph);
	weight_pa = iterator_property_map<int*, EdgeIndexMap, int, int&>(weight_array, edge_id);
}

dict clustering::graph2dict(Graph &subgraph, Graph &graph)
{
	dict graph_dict;
	boost::property_map<Graph, vertex_index_t>::type
	vertex_id = get(vertex_index, graph);
	boost::graph_traits<Graph>::vertex_descriptor
	vertex_local, vertex_local1, vertex_global, vertex_global1;
	#if defined(DEBUG)
		std::cout << "edges(g) = ";
	#endif
	graph_traits<Graph>::edge_iterator ei, ei_end;
	for (tie(ei,ei_end) = edges(subgraph); ei != ei_end; ++ei)
	{
		vertex_local = source(*ei, subgraph);
		vertex_local1 = target(*ei, subgraph);
		vertex_global = subgraph.local_to_global(vertex_local);
		vertex_global1 = subgraph.local_to_global(vertex_local1);
		#if defined(DEBUG)
			std::cout << "(" << get(vertex_id, vertex_global)
			<< "," << get(vertex_id, vertex_global1) << ") ";
		#endif
		graph_dict[boost::python::make_tuple(get(vertex2name, vertex_global), get(vertex2name, vertex_global1))] = get(weight_pa, *ei);
		#if defined(DEBUG)
			std::cout << "[" << get(vertex2name, vertex_global)
			<< "," << get(vertex2name, vertex_global1) << "] " <<get(weight_pa, *ei)<<" weight";
		#endif
	}
	#if defined(DEBUG)
		std::cout << std::endl;
	#endif
	return graph_dict;
}

gsl_matrix* clustering::graph2gsl_matrix(Graph &graph)
{
	int dimension = num_vertices(graph);
	#if defined(DEBUG)
		std::cout<<"Dimension of the graph is "<<dimension<<std::endl;
	#endif
	gsl_matrix* m = gsl_matrix_calloc(dimension, dimension);	//calloc sets all elements to 0, different from alloc
	boost::property_map<Graph, vertex_index_t>::type
	vertex_id = get(vertex_index, graph);
	
	int index1,index2;
	graph_traits<Graph>::edge_iterator ei, ei_end;
	for (tie(ei, ei_end) = edges(graph); ei!=ei_end; ++ei)
	{
		index1 = get(vertex_id, source(*ei, graph));
		index2 = get(vertex_id, target(*ei, graph));
		gsl_matrix_set(m, index1, index2, 1.0);
		gsl_matrix_set(m, index2, index1, 1.0);	//undirected, symmetric
	}
	return m;
}


gsl_vector* clustering::return_eigen_vector(gsl_matrix* graph_matrix, int which_eigen_vector)
{
	gsl_vector *eval = gsl_vector_alloc (graph_matrix->size1);
	gsl_matrix *evec = gsl_matrix_alloc (graph_matrix->size1, graph_matrix->size2);

	gsl_eigen_symmv_workspace * w =
		gsl_eigen_symmv_alloc(graph_matrix->size1);
	gsl_eigen_symmv (graph_matrix, eval, evec, w);
	gsl_eigen_symmv_free (w);
	gsl_eigen_symmv_sort (eval, evec,
						  GSL_EIGEN_SORT_VAL_ASC);
	gsl_vector* evec_i = gsl_vector_alloc(graph_matrix->size1);
	if(gsl_matrix_get_col(evec_i, evec, which_eigen_vector))
	{
		std::cerr<<"error occured when get the eigenvector"<<std::endl;
		return NULL;
	};
	gsl_vector_free(eval);
	gsl_matrix_free(evec);
	return evec_i;
}

double clustering::connectivity_of_graph(Graph &graph)
{

	int no_of_vertices = num_vertices(graph);
	int no_of_edges = num_edges(graph);
	if(no_of_vertices<=1)
		return 0;
		//only one vertex, cause floating point exception
	double connectivity = 2*no_of_edges/(no_of_vertices*(no_of_vertices-1));
	return connectivity;
}

std::vector<Graph> clustering::subgraph_components(Graph &subgraph, Graph &graph, std::vector<int> component, int no_of_components)
{
	//initialize the vector_subgraph with the number of components
	std::vector<Graph> vector_subgraph(no_of_components);
	for(int i=0;i<no_of_components;i++)
		vector_subgraph[i] = graph.create_subgraph();
	//two vertex_descriptor
	boost::graph_traits<Graph>::vertex_descriptor
	vertex_local, vertex_global;
	//the vertex_index_global map is used to translate the global descriptor to the global index.
	boost::property_map<Graph, vertex_index_t>::type
	vertex_index_global;
	vertex_index_global = get(vertex_index, graph);
	for(int i=0; i<component.size(); i++)
	{
		int component_no = component[i];
		//i is the local index, get a descriptor from it
		vertex_local = vertex(i, subgraph);
		//find the global descriptor
		vertex_global = subgraph.local_to_global(vertex_local);
		//get the global index and add it to the subgraph
		add_vertex(get(vertex_index_global, vertex_global), vector_subgraph[component_no]);
	}
	
	return vector_subgraph;
}


void clustering::cluster(Graph &graph)
{
	gsl_matrix* m = graph2gsl_matrix(graph);
	#if defined(DEBUG)
		std::cout<<"The matrix is "<<std::endl;
		gsl_matrix_fprintf(stdout, m, "%g");		//check matrix
	#endif
	gsl_vector* evec_i = return_eigen_vector(m, eigen_vector_no);
	#if defined(DEBUG)
		gsl_vector_fprintf (stdout,evec_i, "%g");
	#endif
	int i;
	std::cout<<"Second minimum eigenvector: "<<std::endl;
	for(i=0;i<evec_i->size;++i)
		std::cout<<gsl_vector_get(evec_i, i)<<"\t";
	std::cout<<std::endl;
	//split the big graph based on eigenvector, >0 or <0
	std::vector<Graph> vector_subgraph(2);
	vector_subgraph[0] = graph.create_subgraph();
	vector_subgraph[1] = graph.create_subgraph();
	for(i=0;i<evec_i->size;++i)
	{
		if (gsl_vector_get(evec_i, i)<0)
			add_vertex(i, vector_subgraph[0]);
		else
			add_vertex(i, vector_subgraph[1]);
	}
	
	for(i=0; i<2; i++)
	{
		int num_vertices_of_subgraph = num_vertices(vector_subgraph[i]);
		if (num_vertices_of_subgraph < min_cluster_size)
			//stop here, too small, even it could be empty
			continue;
		//get all the components and check
		std::vector<int> component(num_vertices_of_subgraph);
		int no_of_components = connected_components(vector_subgraph[i], &component[0]);
		//the second parameter is g, not graph
		std::vector<Graph> vector_sub_subgraph = subgraph_components(vector_subgraph[i], g, component, no_of_components);
		
		std::vector<Graph>::iterator g_iterator;
		std::cout<<"No. of components in subgraph "<<i<<" is: "<<no_of_components<<std::endl;
		for(g_iterator=vector_sub_subgraph.begin();g_iterator!=vector_sub_subgraph.end();++g_iterator)
		{
			if(num_vertices(*g_iterator)>=min_cluster_size)
			{
				double connectivity = connectivity_of_graph(*g_iterator);
				if (connectivity>=connectivity_cutoff)
				{
					//the second parameter is g not graph
					good_clusters.append(graph2dict(*g_iterator, g));
				}
				else
					cluster(*g_iterator);
			}
		}
	}
}

void clustering::run(dict graph_dict)
{
	init_graph_from_dict(graph_dict, g);
	cluster(g);
}

BOOST_PYTHON_MODULE(clustering)
{
	class_<clustering>("clustering", init<double, int, int>())
		.def("init_graph_from_dict", &clustering::init_graph_from_dict)
		.def("run", &clustering::run)
		.def_readonly("good_clusters", &clustering::good_clusters)
	;
}
