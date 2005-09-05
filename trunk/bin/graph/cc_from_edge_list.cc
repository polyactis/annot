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

BOOST_PYTHON_MODULE(cc_from_edge_list)
{
	class_<cc_from_edge_list>("cc_from_edge_list")
		.def("init_graph_from_edge_list", &cc_from_edge_list::init_graph_from_edge_list)
		.def("run", &cc_from_edge_list::run)
		.def_readonly("cc_list", &cc_from_edge_list::cc_list)
	;
}
