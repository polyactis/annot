/*
*
*01-09-06
*      find edge items related to a fim dataset signature, a module of MpiFromDatasetSignatureToPattern.py
*/


#include "PostFim.h"

PostFim::PostFim(int no_cc, int no_of_datasets, int min_cluster_size, std::string node_outputfname)\
	:_no_cc(no_cc), _no_of_datasets(no_of_datasets), _min_cluster_size(min_cluster_size),_node_outputfname(node_outputfname)
{
	#if defined(DEBUG)
		std::cerr<<"PostFim Begins"<<std::endl;
	#endif
	_out.open(node_outputfname.c_str());
}

PostFim::~PostFim()
{
	_out.close();
	#if defined(DEBUG)
		std::cerr<<"PostFim Exists"<<std::endl;
	#endif
}

void PostFim::add_edge_sig_vector(boost::python::list edge_sig_list)
{
	boost::dynamic_bitset<> sig_bitset(_no_of_datasets);
	//int len = extract<int>(edge_sig_list.attr("__len__")());
	for(int i=0; i<_no_of_datasets+2;i++)
	{
		if ((i==0)||(i==1))	//first 2 digits are for edge tuple
			edge_tuple_vector.push_back(extract<unsigned int>(edge_sig_list[i]));
		else
			sig_bitset[i-2] = extract<int>(edge_sig_list[i]);
	}
	edge_bitset_vector.push_back(sig_bitset);
}

void PostFim::add_pattern_signature(boost::python::list pattern_sig_list)
{
	boost::dynamic_bitset<> sig_bitset(_no_of_datasets);
	int len = extract<int>(pattern_sig_list.attr("__len__")());
	for(int i=0; i<len-1;i++)
		sig_bitset[extract<int>(pattern_sig_list[i])-1] = 1;	//in pattern_sig_list, it starts from 1
	pattern_bitset_vector.push_back(sig_bitset);
	freq_of_signature_vector.push_back(extract<int>(pattern_sig_list[len-1]));
}


void PostFim::patternFormation()
{
	std::vector<int> edge_id_vector;
	for (int i = 0; i<pattern_bitset_vector.size(); i++)
	{
		for (int j = 0; j<edge_bitset_vector.size(); j++)
		{
			if ((pattern_bitset_vector[i] & edge_bitset_vector[j])==pattern_bitset_vector[i])
				edge_id_vector.push_back(j);
		}
		if (freq_of_signature_vector[i]==edge_id_vector.size())
		{
			outputCcFromEdgeList(edge_id_vector, edge_tuple_vector, _min_cluster_size, _no_cc);
		}
		else
		{
			std::cerr<<"Warning: signature "<<pattern_bitset_vector[i]<<" frequency "<<freq_of_signature_vector[i]\
				<<" doesn't  match no of edges "<<edge_id_vector.size()<<std::endl;
		}
		edge_id_vector.clear();	//clear it for the next pattern
	}
	
	//clear all patterns, pattern_bitset_vector and freq_of_signature_vector
	pattern_bitset_vector.clear();
	freq_of_signature_vector.clear();
	
}

void PostFim::outputCcFromEdgeList(std::vector<int> &edge_id_vector, \
	std::vector<unsigned int > &edge_tuple_vector, int min_cluster_size, int no_cc)
/*
*01-09-06
*	most functions used here copied from  cc_from_edge_list
*	no_cc is not handled right now
*	--init_graph_from_edge_tuple_vector()
*	--cc2subgraph()
*	(loop)
*		--output_subgraph()
*/
{
	#ifdef DEBUG
		std::cerr<<"Get cc from edge list and output..."<<std::endl;
	#endif
	Graph g = init_graph_from_edge_tuple_vector(edge_id_vector, edge_tuple_vector);
	std::vector<Graph> vector_subgraph = cc2subgraph(g);
	std::vector<Graph>::iterator g_iterator;
	for(g_iterator=vector_subgraph.begin();g_iterator!=vector_subgraph.end();++g_iterator)
		if (num_vertices(*g_iterator)>=min_cluster_size)
			output_subgraph(_out, *g_iterator, g);
	#if defined(DEBUG)
		std::cerr << "done" <<std::endl;
	#endif
}

std::vector<Graph> PostFim::cc2subgraph(Graph &graph)
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

Graph PostFim::init_graph_from_edge_tuple_vector(std::vector<int> &edge_id_vector, std::vector<unsigned int > &edge_tuple_vector)
/*
*01-09-06
*	for PostFim.cc
*
*/
{
	Graph graph;
	std::map<int, vertexDescriptor> gene_no2vertexDescriptor;
	std::map<int, vertexDescriptor>::iterator pos;
	bool inserted;
	vertexDescriptor u, v;
	unsigned int gene1, gene2;
	vertex_name_type vertex2name = get(vertex_name, graph);
	for (int i=0; i<edge_id_vector.size(); i++)
	{
		tie(pos, inserted) = gene_no2vertexDescriptor.insert(std::make_pair(gene1, vertexDescriptor()));
		if (inserted) {
			u = add_vertex(graph);
			vertex2name[u] = gene1;
			pos->second = u;
		} else
			u = pos->second;
		
		tie(pos, inserted) = gene_no2vertexDescriptor.insert(std::make_pair(gene2, vertexDescriptor()));
		if (inserted) {
			v = add_vertex(graph);
			vertex2name[v] = gene2;
			pos->second = v;
		} else
			v = pos->second;
		add_edge(u, v, graph);
	}
	return graph;
}

void PostFim::output_subgraph(ofstream &outf, Graph &subgraph, Graph &graph)
/*
*01-09-06
*	for PostFim.cc
*	both vertex list and edge list are sorted
*/
{
	#if defined(DEBUG)
		std::cerr << "outputting...";
	#endif
	vertex_name_type vertex2name = get(vertex_name, graph);	//mapping
	
	boost::graph_traits<Graph>::vertex_descriptor
		vertex_local, vertex_local1, vertex_global, vertex_global1;
	int vg_name, vg_name1;
	std::pair<vertexIterator, vertexIterator> vp;
	int no_of_vertices = num_vertices(subgraph);
	std::vector<unsigned int> vertex_vector;
	//put vertices into vertex_array
	for (vp = vertices(subgraph); vp.first != vp.second; ++vp.first)
	{
		vertex_global = subgraph.local_to_global(*vp.first);
		vertex_vector.push_back(get(vertex2name, vertex_global));
	}
	sort(vertex_vector.begin(), vertex_vector.end());	//sort it in ascending order
	
	outf<<"[";
	for (int i=0; i<no_of_vertices; i++)
	{
		outf<<vertex_vector[i];
		if (i!= no_of_vertices-1)
			outf << ", ";	//the difference from the middle to the end
	}		
	outf<<"]\t";
	
	const int no_of_edges = num_edges(subgraph);
	std::vector<boost::array<unsigned int, 2> > edge_array_vector;
	boost::array<unsigned int, 2> edge_array;
	graph_traits<Graph>::edge_iterator ei, ei_end;
	//put edges into edge_array_vector
	for (tie(ei,ei_end) = edges(subgraph); ei != ei_end; ++ei)
	{
		vertex_local = source(*ei, subgraph);
		vertex_local1 = target(*ei, subgraph);
		vertex_global = subgraph.local_to_global(vertex_local);
		vertex_global1 = subgraph.local_to_global(vertex_local1);
		vg_name = get(vertex2name, vertex_global);
		vg_name1 = get(vertex2name, vertex_global1);
		if (vg_name<vg_name1)	//in ascending order
		{
			edge_array[0] = vg_name;
			edge_array[1] = vg_name1;
		}
		else
		{
			edge_array[0] = vg_name1;
			edge_array[1] = vg_name;
		}
		edge_array_vector.push_back(edge_array);
	}
	
	sort(edge_array_vector.begin(), edge_array_vector.end(), cmp_edge_array);	//sort the edges
	outf<<"[";
	for (int i=0; i<no_of_edges; i++)
	{
		outf << "(" << edge_array_vector[i][0] << ", " << edge_array_vector[i][1] << ")";
		if (i != no_of_edges-1)	//the difference from the middle to the end
			outf<<", ";
	}
	outf<<"]"<<std::endl;
	#if defined(DEBUG)
		std::cerr << "done" <<std::endl;
	#endif
}


BOOST_PYTHON_MODULE(PostFim)
{
	class_<PostFim, boost::noncopyable>("PostFim", init<int, int, int, std::string>())
		.def("add_edge_sig_vector", &PostFim::add_edge_sig_vector)
		.def("add_pattern_signature", &PostFim::add_pattern_signature)
		.def("patternFormation", &PostFim::patternFormation)
	;
}

/*

void print_usage(FILE* stream,int exit_code)
{
	assert(stream !=NULL);
        fprintf(stream,"Usage: PostFim options inputfile\n");
	fprintf(stream,"\t-h  --help	Display the usage infomation.\n"\
		"\t-o ..., --output=...	Write output to file, gph_output(default)\n"\
		"\t-n ..., --name=...	Same as output_filename(default)\n"\
		"\t-p .., --p_value_cut_off=...	p_value significance cutoff,0.01(default)\n"\
		"\t-c ..., --cor_cut_off=...	correlation cutoff, 0.6(default)\n"\
		"\t\tif p_value_cut_off=0, this cut_off is used instead.\n"\
		"\t-t ..., --top_percentage=...	0.01(default).\n"\
		"\t\t if p_value_cut_off=0 and cor_cut_off=0, top_percentage is used to select edges.\n"\
		"\t-d ..., --max_degree=...	maximum degree of freedom(#columns-2), 10000,(default).\n"\
		"\t-l, --leave_one_out	leave_one_out.\n"\
		"\tFor long option, = or ' '(blank) is same.\n"\
		"\tLine tokenizer is one space, tab, or \\r\n");
	exit(3);
}

int main(int argc, char* argv[])
{
	int next_option;
	const char* const short_options="ho:i:s:";
	const struct option long_options[]={
	  {"help",0,NULL,'h'},
	  {"output",1,NULL,'o'},
	  {"input",1,NULL,'i'},
	  {"sig_vector", 1, NULL, 's'},
	  {NULL,0,NULL,0}
	};
	
	char* output_filename = "PostFim.output";
	char* input_filename = NULL;
	char* sig_vector_fname = NULL;

	do
	{
		next_option=getopt_long(argc,argv,short_options,long_options,NULL);
		switch(next_option)
		{
		case 'h':
			print_usage(stdout,0);
	  		exit(1);
		case 'o':
			output_filename=optarg;
			break;
		case 'i':
			input_filename = optarg;
			break;
		case 's':
			sig_vector_fname = optarg;
			break;
		case '?':
			print_usage(stderr,-1);
		case -1:
			break;
		default:
			abort();
		}
	}while(next_option!=-1);

	if (input_filename!=NULL && sig_vector_fname!=NULL)
	{
		
		graph_construct instance(input_filename, sig_vector_fname, output_filename);
		instance.run();
		
	}
	else
		print_usage(stderr, 1);
}
*/
