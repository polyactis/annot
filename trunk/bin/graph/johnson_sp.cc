#include "johnson_sp.h"

johnson_sp::johnson_sp()
	:_no_of_datasets(0)
{
}

johnson_sp::johnson_sp(int no_of_datasets)
	:_no_of_datasets(no_of_datasets)
{
}

void johnson_sp::add_edge_sig_vector(boost::python::list edge_sig_list)
/*
01-23-06
*/
{
	boost::dynamic_bitset<> sig_bitset(_no_of_datasets);
	//int len = extract<int>(edge_sig_list.attr("__len__")());
	unsigned long edge_tuple[2] = {0,0};
	for(int i=0; i<_no_of_datasets+2;i++)
	{
		if ((i==0)||(i==1))	//first 2 digits are for edge tuple
			edge_tuple[i] = extract<unsigned long>(edge_sig_list[i]);
		else
			sig_bitset[i-2] = extract<int>(edge_sig_list[i]);
	}
	//in ascending order, edge names in *.sig_vector are actually in reverse order, due to Haiyan's remnancy
	if (edge_tuple[0]<edge_tuple[1])
		edge2bitset[hash_edge_name(edge_tuple[0], edge_tuple[1])] = sig_bitset;
	else
		edge2bitset[hash_edge_name(edge_tuple[1], edge_tuple[0])] = sig_bitset;
}

std::pair<Graph, std::vector<boost::dynamic_bitset<> > > johnson_sp::init_graph_from_vertex_edge_list(boost::python::list vertex_list, boost::python::list edge_list)
/*
01-23-06
	move the graph from input argument to return argument
	all gene_no-related int changed to unsigned long
	add codes to fill recurrence_bitset_vector
*/
{
	#ifdef DEBUG
		std::cerr<<"Read in graph"<<std::endl;
	#endif
	Graph graph;
	std::vector<boost::dynamic_bitset<> > recurrence_bitset_vector;	//01-23-06
	std::map<unsigned long, vertexDescriptor> geneNoMap;
	int no_of_vertices = boost::python::extract<int>(vertex_list.attr("__len__")());
	vertexDescriptor u, v;
	for(int i=0; i<no_of_vertices; i++)
	{
		unsigned long gene_no = boost::python::extract<unsigned long>(vertex_list[i]);
		u = add_vertex(graph);
		geneNoMap[gene_no] = u;
	}
	int no_of_edges = boost::python::extract<int>(edge_list.attr("__len__")());
	property_map<Graph, edge_weight_t>::type weightmap = get(edge_weight, graph);
	bool inserted;
	edgeDescriptor e;
	for(int i=0; i<no_of_edges; i++)
	{
		unsigned long gene1 = boost::python::extract<unsigned long>(edge_list[i][0]);
		unsigned long gene2 = boost::python::extract<unsigned long>(edge_list[i][1]);
		u = geneNoMap[gene1];
		v = geneNoMap[gene2];
		tie(e, inserted) = add_edge(u, v, graph);
		weightmap[e] = 1;	//weight is necessary, otherwise, johnson_all_pairs_shortest_paths regards it as infiniti
		if (_no_of_datasets!=0)
			recurrence_bitset_vector.push_back(edge2bitset[hash_edge_name(gene1, gene2)]);	//gene1 and gene2 are already in ascending order
	}
	#ifdef DEBUG
		std::cerr<<"Graph readin done."<<std::endl;
	#endif
	return std::make_pair(graph, recurrence_bitset_vector);
}


boost::python::list johnson_sp::calculate_sp(Graph &graph)
/*
01-23-06
	move the D from input argument to return argument
*/
{
	#ifdef DEBUG
		std::cerr<<"Calculating SP..."<<std::endl;
	#endif
	boost::python::list D_matrix;
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
	return D_matrix;
}


boost::python::list johnson_sp::py_shortest_distance(boost::python::list vertex_list, boost::python::list edge_list)
/*
01-23-06
	named from run()
	lots of changes following others
*/
{
	std::pair<Graph, std::vector<boost::dynamic_bitset<> > > graph_recurrence_bitset_vector;
	graph_recurrence_bitset_vector = init_graph_from_vertex_edge_list(vertex_list, edge_list);
	_recurrence_bitset_vector = graph_recurrence_bitset_vector.second;
	return calculate_sp(graph_recurrence_bitset_vector.first);
}

boost::python::list johnson_sp::py_recurrence_list()
/*
01-23-06
	calculate the recurrence_list from _recurrence_bitset_vector
	
	must be run after init_graph_from_vertex_edge_list()
*/
{
	#ifdef DEBUG
		std::cerr<<"calculating py_recurrence_list..."<<std::endl;
	#endif
	boost::python::list recurrence_list;
	for (int i=0; i<_no_of_datasets; i++)
	{
		float recurrence = 0.0;
		for (int j=0; j< _recurrence_bitset_vector.size(); j++)
			recurrence += _recurrence_bitset_vector[j][i];	//watch the order of i and j
		recurrence /= _recurrence_bitset_vector.size();	//divided by no_of_edges
		recurrence_list.append(recurrence);
	}
	return recurrence_list;
	#ifdef DEBUG
		std::cerr<<"done."<<std::endl;
	#endif
}

BOOST_PYTHON_MODULE(johnson_sp)
{
	class_<johnson_sp>("johnson_sp")
		.def(init<int>())
		.def("add_edge_sig_vector", &johnson_sp::add_edge_sig_vector)
		.def("py_shortest_distance", &johnson_sp::py_shortest_distance)
		.def("py_recurrence_list", &johnson_sp::py_recurrence_list)
	;

}
