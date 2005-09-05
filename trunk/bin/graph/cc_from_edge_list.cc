/*
*04-29-05
*	a module to extract connected components(cc) from an edge list.
*09-04-05
*	class definition and included header files go to cc_from_edge_list.h
*/

#include "cc_from_edge_list.h"

void cc_from_edge_list::init_graph_from_edge_list(boost::python::list edge_list, Graph &graph)
{
	vertex_name_type vertex2name = get(vertex_name, graph);
	int edge_length = boost::python::extract<int>(edge_list.attr("__len__")());
	for(int i=0; i<edge_length; i++)
	{
		int gene1 = boost::python::extract<int>(edge_list[i][0]);
		int gene2 = boost::python::extract<int>(edge_list[i][1]);
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
	}
}

std::vector<Graph> cc_from_edge_list::cc2subgraph(Graph &graph)
{
	std::vector<int> component(num_vertices(graph));
	int no_of_components = connected_components(graph, &component[0]);
	//initialize the vector_subgraph with the number of components
	std::vector<Graph> vector_subgraph(no_of_components);
	for(int i=0;i<no_of_components;i++)
		vector_subgraph[i] = graph.create_subgraph();
	//vertex_descriptor
	boost::graph_traits<Graph>::vertex_descriptor vertex_of_graph;

	for(int i=0; i<component.size(); i++)
	{
		int component_no = component[i];
		add_vertex(i, vector_subgraph[component_no]);
	}
	
	return vector_subgraph;
}

boost::python::list cc_from_edge_list::subgraph2list(Graph &subgraph, Graph &graph)
{
	boost::python::list subgraph_edge_list;	//store the egdes
	
	vertex_name_type vertex2name = get(vertex_name, graph);	//mapping
	
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
			std::cout << "(" << get(vertex2name, vertex_global)
			<< "," << get(vertex2name, vertex_global1) << ") ";
		#endif
		subgraph_edge_list.append(boost::python::make_tuple(get(vertex2name, vertex_global), get(vertex2name, vertex_global1)));

	}
	#if defined(DEBUG)
		std::cout << std::endl;
	#endif
	return subgraph_edge_list;
}

std::vector<Graph> cc_from_edge_list::subgraph_components(Graph &subgraph, Graph &graph, std::vector<int> component, int no_of_components)
/*
*09-05-05
*	copied from clustering.cc, dynamic .so can't find the definition of the function by just including the header file
*/
{
	int global_index;
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
		global_index = get(vertex_index_global, vertex_global);
		#if defined(DEBUG)
			std::cerr<<global_index<<" goes to component "<<component_no<<std::endl;
		#endif
		add_vertex(global_index, vector_subgraph[component_no]);
	}
	
	return vector_subgraph;
}

std::vector<Graph> cc_from_edge_list::graph_components(Graph &graph, std::vector<int> component, int no_of_components)
/*
*09-05-05
*	create component-graph-vector from a graph
*/
{
	vertex_name_type global_vertex2name = get(vertex_name, graph);
	std::map<int, vertexDescriptor> vertex_int2comp_graph_descriptor;	//to store the global vertex_index to local vertex_index
	std::vector<Graph> vector_graph(no_of_components);
	for(int i=0;i<no_of_components;i++)
	{
		Graph tmp_g;
		vector_graph[i] = tmp_g;
	}
	vertexDescriptor u, vertex1, vertex2;
	for(int i=0; i<component.size(); i++)
	{
		int component_no = component[i];
		u = add_vertex(vector_graph[component_no]);
		//preserve the name mapping
		vertex_name_type local_vertex2name = get(vertex_name, vector_graph[component_no]);
		local_vertex2name[u] = get(global_vertex2name, vertex(i,graph));
		//global to local vertex_index mapping
		vertex_int2comp_graph_descriptor[i] = u;
	}
	
	graph_traits<Graph>::edge_iterator ei, ei_end;
	for (tie(ei,ei_end) = edges(graph); ei != ei_end; ++ei)
	{
		vertex1 = source(*ei, graph);
		vertex2 = target(*ei, graph);
		int component_no = component[vertex1];
		#ifdef DEBUG
			int component_no2 = component[vertex2];
			if (component_no!=component_no2)
				std::cerr<<"Error: two vertices of one edge from two different component"<<std::endl;
		#endif
		add_edge(vertex_int2comp_graph_descriptor[vertex1], vertex_int2comp_graph_descriptor[vertex2], vector_graph[component_no]);
	}
	return vector_graph;
}


void cc_from_edge_list::run(boost::python::list edge_list)
{
	init_graph_from_edge_list(edge_list, g);
	vector_subgraph = cc2subgraph(g);
	std::vector<Graph>::iterator g_iterator;
	#ifdef DEBUG
		std::cout<<"No. of components in graph is: "<<vector_subgraph.size()<<std::endl;
	#endif
	for(g_iterator=vector_subgraph.begin();g_iterator!=vector_subgraph.end();++g_iterator)
	{
		boost::python::list subgraph_edge_list = subgraph2list(*g_iterator, g);
		cc_list.append(subgraph_edge_list);

	}
}

/*
*09-04-05
*	a clustering algorithm based  on edge betweenness_centrality
*
*/


void ClusterByEBC::cut_by_betweenness_centrality(Graph &graph, int size_cutoff, float conn_cutoff)
/*
*09-04-05
*	cut the graph with the edge of maximum betweenness_centrality
*/
{
	int no_of_vertices = num_vertices(graph);
	int no_of_edges = num_edges(graph);
	#ifdef DEBUG
		std::cerr<<"no_of_vertices: "<<no_of_vertices<<std::endl;
		std::cerr<<"no_of_edges: "<<no_of_edges<<std::endl;
	#endif
	if (no_of_vertices>=size_cutoff)
	{
		float connectivity = 2.0*no_of_edges/(no_of_vertices*(no_of_vertices-1));
		#ifdef DEBUG
			std::cerr<<"connectivity: "<<connectivity<<std::endl;
		#endif
		if (connectivity>=conn_cutoff)
		{
			#ifdef DEBUG
				std::cerr<<"good subgraph "<<std::endl;
			#endif
			good_subgraph_vector.push_back(graph);
		}
		else
		{
			vector<double> edge_centrality(no_of_edges);
			EdgeCentralityMap ec_map(edge_centrality.begin(), get(edge_index, graph));
				//"make_iterator_property_map(edge_centrality.begin(), get(edge_index, subgraph), double())" also works.
			indirect_cmp<EdgeCentralityMap, std::less<centrality_type> > cmp(ec_map);
			#ifdef DEBUG
				std::cerr<<"running brandes_betweenness_centrality... ";
			#endif
			brandes_betweenness_centrality(graph, edge_centrality_map(ec_map)); 
			#ifdef DEBUG
				std::cerr<<"done."<<std::endl;
			#endif
			edgeDescriptor e = *max_element(edges(graph).first, edges(graph).second, cmp);
			centrality_type max_centrality = get(ec_map, e);
			#ifdef DEBUG
				std::cerr<<"max_centrality is "<<max_centrality<<std::endl;
			#endif
			remove_edge(e, graph);
			#ifdef DEBUG
				std::cerr<<"after removal the subgraph has "<<num_edges(graph)<<" edges."<<std::endl;
			#endif
			std::vector<int> component(num_vertices(graph));
			int no_of_components = connected_components(graph, &component[0]);
			if (no_of_components==1)	//keep cutting
			{
				#ifdef DEBUG
					std::cerr<<"only one component, keep cutting "<<std::endl;
				#endif
				cut_by_betweenness_centrality(graph, size_cutoff, conn_cutoff);
			}
			else	//first get the connected_components into a vector_sub_subgraph and cut them separately
			{
				#ifdef DEBUG
					std::cerr<<"get the connected_components into a vector_graph and cut them separately "<<std::endl;
				#endif
				std::vector<Graph> vector_graph = graph_components(graph, component, no_of_components);
				std::vector<Graph>::iterator g_iterator;
				for(g_iterator=vector_graph.begin();g_iterator!=vector_graph.end();++g_iterator)
				{
					cut_by_betweenness_centrality(*g_iterator, size_cutoff, conn_cutoff);
				}
			}
		}
	}
	#ifdef DEBUG
	else
		std::cerr<<"Bad graph, too small,"<<no_of_vertices<<"vertices"<<std::endl;
	#endif
	
}

boost::python::list ClusterByEBC::graph2list(Graph &graph)
/*
*09-05-05
*	similar to subgraph2list()
*/
{
	boost::python::list graph_edge_list;	//store the egdes
	
	vertex_name_type vertex2name = get(vertex_name, graph);	//mapping
	
	boost::graph_traits<Graph>::vertex_descriptor
	vertex1, vertex2;
	#if defined(DEBUG)
		std::cout << "edges(g) = ";
	#endif
	graph_traits<Graph>::edge_iterator ei, ei_end;
	for (tie(ei,ei_end) = edges(graph); ei != ei_end; ++ei)
	{
		vertex1 = source(*ei, graph);
		vertex2 = target(*ei, graph);
		#if defined(DEBUG)
			std::cout << "(" << get(vertex2name, vertex1)
			<< "," << get(vertex2name, vertex2) << ") ";
		#endif
		graph_edge_list.append(boost::python::make_tuple(get(vertex2name, vertex1), get(vertex2name, vertex2)));

	}
	#if defined(DEBUG)
		std::cout << std::endl;
	#endif
	return graph_edge_list;
}

void ClusterByEBC::run(boost::python::list edge_list, int size_cutoff, float conn_cutoff)
/*
*09-04-05
*	override the run() of the parent class
*/
{
	#ifdef DEBUG
		std::cerr<<"starting ClusterByEBC...";
	#endif
	init_graph_from_edge_list(edge_list, g);
	//get the connected_components out of g
	std::vector<int> component(num_vertices(g));
	int no_of_components = connected_components(g, &component[0]);
	std::vector<Graph> vector_graph = graph_components(g, component, no_of_components);

	std::vector<Graph>::iterator g_iterator;
	#ifdef DEBUG
		std::cerr<<"No. of components in graph is: "<<vector_graph.size()<<std::endl;
	#endif
	for(g_iterator=vector_graph.begin();g_iterator!=vector_graph.end();++g_iterator)
	{
		cut_by_betweenness_centrality(*g_iterator, size_cutoff, conn_cutoff);
	}
	for(g_iterator=good_subgraph_vector.begin();g_iterator!=good_subgraph_vector.end();++g_iterator)
	{
		boost::python::list graph_edge_list = graph2list(*g_iterator);
		cc_list.append(graph_edge_list);
	}
}

BOOST_PYTHON_MODULE(cc_from_edge_list)
{
	class_<cc_from_edge_list>("cc_from_edge_list")
		.def("init_graph_from_edge_list", &cc_from_edge_list::init_graph_from_edge_list)
		.def("run", &cc_from_edge_list::run)
		.def_readonly("cc_list", &cc_from_edge_list::cc_list)
	;
	
	class_<ClusterByEBC, bases<cc_from_edge_list> >("ClusterByEBC")
		.def("run", &ClusterByEBC::run)
		.def_readonly("cc_list", &ClusterByEBC::cc_list)
	;
}
