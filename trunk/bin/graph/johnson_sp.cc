#include "johnson_sp.h"


void johnson_sp::init_graph_from_vertex_edge_list(boost::python::list vertex_list, boost::python::list edge_list, Graph &graph)
{
	#ifdef DEBUG
		std::cerr<<"Read in graph"<<std::endl;
	#endif
	int no_of_vertices = boost::python::extract<int>(vertex_list.attr("__len__")());
	vertexDescriptor u, v;
	for(int i=0; i<no_of_vertices; i++)
	{
		int gene_no = boost::python::extract<int>(vertex_list[i]);
		u = add_vertex(graph);
		geneNoMap[gene_no] = u;
	}
	int no_of_edges = boost::python::extract<int>(edge_list.attr("__len__")());
	property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, g);
	bool inserted;
	edgeDescriptor e;
	for(int i=0; i<no_of_edges; i++)
	{
		int gene1 = boost::python::extract<int>(edge_list[i][0]);
		int gene2 = boost::python::extract<int>(edge_list[i][1]);
		u = geneNoMap[gene1];
		v = geneNoMap[gene2];
		tie(e, inserted) = add_edge(u, v, graph);
		weightmap[e] = 1;	//weight is necessary, otherwise, johnson_all_pairs_shortest_paths regards it as infiniti
	}
	#ifdef DEBUG
		std::cerr<<"Graph readin done."<<std::endl;
	#endif
}


void johnson_sp::calculate_sp(Graph &graph, boost::python::list &D_matrix)
{
	#ifdef DEBUG
		std::cerr<<"Calculating SP..."<<std::endl;
	#endif
	int no_of_vertices = num_vertices(graph);
	std::vector<std::vector<int> > d_matrix;
	std::vector<int> one_row;
	for (int i=0; i<no_of_vertices; i++)
	{
		one_row.push_back(-1);
	}
	for (int i=0; i<no_of_vertices; i++)
	{
		d_matrix.push_back(one_row);
	}
	johnson_all_pairs_shortest_paths(graph, d_matrix);
	
	//put it into python list
	for(int i=0; i<d_matrix.size(); i++)
	{
		boost::python::list d_row;
		for(int j=0; j<d_matrix[i].size(); j++)
		{
			#ifdef DEBUG
				std::cerr<<d_matrix[i][j];
			#endif
			d_row.append(d_matrix[i][j]);
		}
		#ifdef DEBUG
			std::cerr<<std::endl;
		#endif
		D_matrix.append(d_row);
	}
	#ifdef DEBUG
		std::cerr<<"SP done."<<std::endl;
	#endif
}


void johnson_sp::run(boost::python::list vertex_list, boost::python::list edge_list)
{
	init_graph_from_vertex_edge_list(vertex_list, edge_list, g);
	calculate_sp(g,D);
}

BOOST_PYTHON_MODULE(johnson_sp)
{
	class_<johnson_sp>("johnson_sp")
		.def("run", &johnson_sp::run)
		.def_readonly("D", &johnson_sp::D)
	;

}
