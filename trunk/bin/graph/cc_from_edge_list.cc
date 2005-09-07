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

ClusterByEBC::ClusterByEBC()
	:_input_filename(""),_min_edge_weight(0),_size_cutoff(5),_conn_cutoff(0.2),_output_filename(""),_format_type(1),_offset(0)
/*
*09-06-05
*	do nothing
*/
{
}

ClusterByEBC::ClusterByEBC(boost::python::list edge_list, int size_cutoff, float conn_cutoff)\
	:_input_filename(""),_min_edge_weight(0), _size_cutoff(size_cutoff),_conn_cutoff(conn_cutoff),_output_filename(""),_format_type(3),_offset(0)
/*
*09-06-05
	Constructor for python to call
	Watch: _format_type(3)
*/
{
	_edge_list = edge_list;
}

ClusterByEBC::ClusterByEBC(std::string input_filename, int min_edge_weight, int size_cutoff, \
	float conn_cutoff, int format_type, int offset, std::string output_filename):\
	_input_filename(input_filename),_min_edge_weight(min_edge_weight),_size_cutoff(size_cutoff),\
	_conn_cutoff(conn_cutoff), _format_type(format_type), _offset(offset),  _output_filename(output_filename)
/*
09-06-05
	Normal constructor
*/
{
}


ClusterByEBC::~ClusterByEBC()
/*
*09-06-06
*	clear the good_subgraph_vector
*/
{
	good_subgraph_vector.clear();
}

void ClusterByEBC::init_graph_from_file(const std::string &input_filename, Graph &graph, const int &min_edge_weight)
/*09-06-05
*	declare weight_array and edge_weight to float
*/
{
	#ifdef DEBUG
		std::cerr<<"Read in graph from "<<input_filename<<" ...";
	#endif
	std::ifstream datafile(input_filename.c_str());
	std::vector<float> weight_array;	//a local weight_array
	std::map<std::string, vertexDescriptor> VertexNameMap;	//to check whether one vertex is already inserted
	vertex_name_type vertex2name = get(vertex_name, graph);
	char_separator<char> sep(" \t");		//blank or '\t' is the separator
	for (std::string line; std::getline(datafile, line);) {
		char_tokenizer line_toks(line, sep);
		char_tokenizer::iterator tokenizer_iter = line_toks.begin();
		if(*tokenizer_iter=="t")
			continue;	//skip the whole line
		if(*tokenizer_iter=="v")
			continue;	//skip the whole line
		if(*tokenizer_iter=="e")
			*tokenizer_iter++;	//skip 'e'
		
		std::string gene1 = *tokenizer_iter++;
		std::string gene2 = *tokenizer_iter++;
		std::string edge_weight_string = *tokenizer_iter;
		float edge_weight = atof(edge_weight_string.c_str());
		if (edge_weight>=min_edge_weight)
		{
			std::map<std::string, vertexDescriptor>::iterator pos;
			bool inserted;
			vertexDescriptor u, v;
			tie(pos, inserted) = VertexNameMap.insert(std::make_pair(gene1, vertexDescriptor()));
			if (inserted) {
				u = add_vertex(graph);
				vertex2name[u] = atoi(gene1.c_str());
				pos->second = u;
			} else
				u = pos->second;
			tie(pos, inserted) = VertexNameMap.insert(std::make_pair(gene2, vertexDescriptor()));
			if (inserted) {
				v = add_vertex(graph);
				vertex2name[v] = atoi(gene2.c_str());
				pos->second = v;
			} else
				v = pos->second;
	
			graph_traits < Graph >::edge_descriptor e;
			tie(e, inserted) = add_edge(u, v, graph);
			if (inserted)
				weight_array.push_back(edge_weight);
		}
	}
	datafile.close();
	#ifdef DEBUG
		std::cerr<<"Done."<<std::endl;
	#endif

}

void ClusterByEBC::init_graph_from_file(const std::string &input_filename, Graph &graph, const int &size_cutoff, const int &offset)
/*
*09-06-05
*	construct the graph from MpiFromDatasetSignatureToPattern.py output file, given an offset
*		size_cutoff is not used.
*/
{
	#ifdef DEBUG
		std::cerr<<"Read in graph from "<<input_filename<<" ...";
	#endif
	std::ifstream datafile(input_filename.c_str());
	std::vector<float> weight_array;	//a local weight_array
	std::map<std::string, vertexDescriptor> VertexNameMap;	//to check whether one vertex is already inserted
	vertex_name_type vertex2name = get(vertex_name, graph);
	char_separator<char> line_sep("\t");		//'\t' is the separator
	char_separator<char> gene_no_list_sep("[], ", "", boost::drop_empty_tokens);	//[ ] , or blank is separator and drop the empty token
	std::string line;
	//skip some lines based on offset
	for (int no_of_lines_skipped=0;no_of_lines_skipped<=offset;std::getline(datafile, line),no_of_lines_skipped++)
	{
		#ifdef DEBUG
			std::cerr<<"lines skipped is "<<no_of_lines_skipped<<std::endl;
		#endif
	}
	//line is the last fetched
	char_tokenizer line_toks(line, line_sep);
	char_tokenizer::iterator tok_iter = line_toks.begin();
	std::string vertex_list_string = *tok_iter++;
	std::string edge_list_string = *tok_iter++;

	//parse the vertex_list_string
	char_tokenizer vertex_list_string_tok(vertex_list_string, gene_no_list_sep);
	vertexDescriptor u, v;
	for(char_tokenizer::iterator tok_iter = vertex_list_string_tok.begin(); tok_iter != vertex_list_string_tok.end(); ++tok_iter)
	{
		std::string gene_no = *tok_iter;
		#ifdef DEBUG
			std::cerr<<"Gene_no is "<<gene_no<<std::endl;
		#endif
		u = add_vertex(graph);
		VertexNameMap.insert(std::make_pair(gene_no, u));
		vertex2name[u] = atoi(gene_no.c_str());
	}
	//parse the edge_list_string
	char_tokenizer edge_list_string_tok(edge_list_string, gene_no_list_sep);
	for(char_tokenizer::iterator tok_iter = edge_list_string_tok.begin(); tok_iter != edge_list_string_tok.end();)	//Watch: no tok_iter++
	{
		std::string gene1 = *tok_iter++;
		std::string gene2 = *tok_iter++;
		#ifdef DEBUG
			std::cerr<<"gene1 and gene2 is "<<gene1<<" and "<<gene2<<std::endl;
		#endif
		u = VertexNameMap[gene1];
		v = VertexNameMap[gene2];
		add_edge(u, v, graph);
	}
	datafile.close();
	#ifdef DEBUG
		std::cerr<<"Done."<<std::endl;
	#endif
}

void ClusterByEBC::reindex_edge(Graph &graph)
/*
09-06-05
	reindex the edges to be sure they are in the range[0, num_edges(graph))
*/
{
	EdgeIndexMap edge2index = get(edge_index, graph);
	graph_traits<Graph>::edge_iterator ei, ei_end;
	int no_of_edges = 0;
	for (tie(ei,ei_end) = edges(graph); ei != ei_end; ++ei)
	{
		edge2index[*ei] = no_of_edges;
		no_of_edges++;
	}
}

void ClusterByEBC::cut_by_betweenness_centrality(Graph &graph, const int &size_cutoff, const float &conn_cutoff)
/*
*09-04-05
*	cut the graph with the edge of maximum betweenness_centrality
09-07-05
	a bug arises in calling brandes_betweenness_centrality(), Doug's email solved it.
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
			reindex_edge(graph);	//09-07-05	fix an important bug here.
			#ifdef DEBUG
				std::cerr<<"after removal the subgraph has "<<num_edges(graph)<<" edges."<<std::endl;
			#endif
			std::vector<int> component(num_vertices(graph));
			int no_of_components = connected_components(graph, &component[0]);
			#ifdef DEBUG
				std::cerr<<"no_of_components: "<<no_of_components<<std::endl;
			#endif
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

twoListTuple ClusterByEBC::graph2list(Graph &graph)
/*
*09-05-05
*	similar to subgraph2list()
*09-05-05
*	transform the vertex set as well.
*/
{
	boost::python::list graph_edge_list;	//store the egdes
	boost::python::list graph_vertex_list;	//store the vertices
	
	vertex_name_type vertex2name = get(vertex_name, graph);	//mapping
	
	std::pair<vertexIterator, vertexIterator> vp;
	#ifdef DEBUG
		std::cout<<"vertices(g) = ";
	#endif
	for (vp = vertices(graph); vp.first != vp.second; ++vp.first)
	{
		graph_vertex_list.append(get(vertex2name, *vp.first));
		#ifdef DEBUG
			std::cout << get(vertex2name, *vp.first) <<  " ";
		#endif
	}
	#if defined(DEBUG)
		std::cout << std::endl;
	#endif
	
	vertexDescriptor vertex1, vertex2;
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
	return boost::make_tuple(graph_vertex_list,graph_edge_list);
}

void ClusterByEBC::output_graph(std::ofstream &outf, Graph &graph)
/*
09-06-05
	output the graph to output_filename
*/
{
	vertex_name_type vertex2name = get(vertex_name, graph);	//mapping
	
	std::pair<vertexIterator, vertexIterator> vp;
	outf<<"[";
	int no_of_vertices=0;
	for (vp = vertices(graph); vp.first != vp.second; ++vp.first)
	{
		no_of_vertices++;
		if (no_of_vertices==num_vertices(graph))	//the difference from the middle to the end
			outf<<get(vertex2name, *vp.first);
		else
			outf << get(vertex2name, *vp.first) <<  ", ";
	}
	outf<<"]\t";
	int no_of_edges =0;
	vertexDescriptor vertex1, vertex2;
	outf<<"[";
	graph_traits<Graph>::edge_iterator ei, ei_end;
	for (tie(ei,ei_end) = edges(graph); ei != ei_end; ++ei)
	{
		no_of_edges++;
		vertex1 = source(*ei, graph);
		vertex2 = target(*ei, graph);
		outf << "[" << get(vertex2name, vertex1)
			<< ", " << get(vertex2name, vertex2) << "]";
		if (no_of_edges!=num_edges(graph))
			outf<<", ";
	}
	outf<<"]"<<std::endl;
}

void ClusterByEBC::run()
/*
*09-04-05
*	override the run() of the parent class
*09-05-05
*	append the cc_vertex_list
	
	--init_graph_from_file()
	or
	--init_graph_from_file()(overloaded)
	or
	--init_graph_from_edge_list()

	--connected_components()
	--graph_components()
	(loop)
		--cut_by_betweenness_centrality()
	(loop)
		--graph2list()
		or
		--output_graph()
*/
{
	if (_format_type==1 && _input_filename!="")
		init_graph_from_file(_input_filename, g, _min_edge_weight);
	else
	{
		if (_format_type==2 &&_input_filename!="")
			init_graph_from_file(_input_filename, g, _size_cutoff, _offset);
		else
		{
			if (_format_type==3 && _input_filename=="")
				init_graph_from_edge_list(_edge_list, g);
			else
			{
				std::cerr<<"Exit: format_type "<<_format_type<<" and input_filename "\
					<<_input_filename<<" combination error."<<std::endl;
				exit(3);
			}
		}
	}
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
		cut_by_betweenness_centrality(*g_iterator, _size_cutoff, _conn_cutoff);
	}
	
	if (_output_filename=="")
	{
		for(g_iterator=good_subgraph_vector.begin();g_iterator!=good_subgraph_vector.end();++g_iterator)
		{
			twoListTuple vertex_edge_tuple = graph2list(*g_iterator);
			cc_vertex_list.append(vertex_edge_tuple.get<0>());
			cc_list.append(vertex_edge_tuple.get<1>());
		}
	}
	else
	{		
		std::ofstream outf(_output_filename.c_str());
		for(g_iterator=good_subgraph_vector.begin();g_iterator!=good_subgraph_vector.end();++g_iterator)
			output_graph(outf, *g_iterator);
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
		.def(init<std::string, int, int, float, boost::python::optional<int, int, std::string> >())
		.def(init<boost::python::list, int, float>())
		.def("run", &ClusterByEBC::run)
		.def_readonly("cc_list", &ClusterByEBC::cc_list)
		.def_readonly("cc_vertex_list", &ClusterByEBC::cc_vertex_list)
	;
}

void print_usage(FILE* stream, char* program_name)
{
	assert(stream !=NULL);
        fprintf(stream,"Usage: %s options -i INPUTFILE\n",program_name);
	fprintf(stream,"\t-h  --help	Display the usage infomation.\n"\
		"\t-i ..., --input=...	INPUTFILE\n"\
		"\t-o ..., --output=...	if not given, no output\n"\
		"\t-s ..., --min_size=...	Min graph size, 5(default)\n"\
		"\t-d ..., --density_cutoff=...	density cutoff, 0.2(default)\n"\
		"\t-e ..., --min_edge_weight=...	minimum edge weight, 0(default).\n"\
		"\t-m ..., --format_type=...	the format type, 1(gspan format, default), 2, 3\n"\
		"\t-f ..., --offset=...	the offset into the inputfile, 0(default)\n"\
		"\tFor long option, = or ' '(blank) is same.\n"\
		"\tFormat 2 is MpiFromDatasetSignatureToPattern.py output format.\n"\
		"\tFormat 3 is to get edge_list from python.\n");
	exit(3);
}


int main(int argc, char* argv[])
{
	int next_option;
	const char* const short_options="hi:o:s:d:e:m:f:";
	const struct option long_options[]={
	  {"help",0,NULL,'h'},
	  {"input", 1, NULL, 'i'},
	  {"output",1,NULL,'o'},
	  {"min_size",1,NULL,'s'},
	  {"density_cutoff", 1, NULL, 'd'},
	  {"min_edge_weight", 1, NULL, 'e'},
	  {"format_type",1,NULL,'m'},
	  {"offset",1,NULL,'f'},
	  {NULL,0,NULL,0}
	};
	
	char* program_name=argv[0];	
	std::string input_filename = "";
	std::string output_filename = "";
	int min_size = 5;
	double density_cutoff = 0.2;
	int min_edge_weight = 0;
	int format_type = 1;
	int offset = 0;

	do
	{
		next_option=getopt_long(argc,argv,short_options,long_options,NULL);
		switch(next_option)
		{
		case 'h':
			print_usage(stdout,0);
	  		exit(1);
		case 'i':
			input_filename = optarg;
			break;
		case 'o':
			output_filename = optarg;
			break;
		case 's':
			min_size = atoi(optarg);
			break;
		case 'd':
			density_cutoff = atof(optarg);
			break;
		case 'e':
			min_edge_weight = atoi(optarg);
			break;
		case 'm':
			format_type = atoi(optarg);
			break;
		case 'f':
			offset = atoi(optarg);
			break;
		case '?':
			print_usage(stderr, program_name);
		case -1:
			break;
		default:
			abort();
		}
	}while(next_option!=-1);
	
	//ifstream inf(argv[1]);
	//ofstream outf(argv[2], ios::app | ios::out);

	if (input_filename!="")
	{
		ClusterByEBC instance(input_filename, min_edge_weight, min_size, density_cutoff, format_type, offset, output_filename);
		instance.run();
	}
	else
		print_usage(stderr, program_name);
}
