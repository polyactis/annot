/*
*
*
*02-24-05
*	boost and boost::python conflict for the make_tuple function. use boost::python::make_tuple instead.*
*05-25-05
*	modify it to be a standalone program and do normalized cut in a way similar to netmine(copath).
*	cut the graph until its size is less than a threshold.
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
#include <boost/tokenizer.hpp>	//for tokenizer, parse input file
#include <boost/tuple/tuple.hpp>

#include <gsl/gsl_math.h>		//for gsl_matrix, gsl_vector stuff
#include <gsl/gsl_eigen.h>		//for gsl_eigen stuff

#include <boost/python.hpp>		//for python module and dict, tuple, make_tuple

#include <getopt.h>	//05-25-05	to parse program options.

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
		clustering(std::string input_filename_option, std::string output_filename_option, int max_size_option,
			int cut_loop_num_option, double density_cutoff_option, int min_edge_weight_option);	//05-25-05
		~clustering();
		void init_graph_from_dict(dict graph_dict, Graph &graph);
		void init_graph_from_file(std::string input_filename, Graph &graph, int min_edge_weight);	//05-25-05
		dict graph2dict(Graph &subgraph, Graph &graph);
		gsl_matrix *graph2gsl_matrix(Graph &graph);
		gsl_vector *return_eigen_vector(gsl_matrix* graph_matrix, int which_eigen_vector);
		double connectivity_of_graph(Graph &graph);
		std::vector<Graph> subgraph_components(Graph &subgraph, Graph &graph, std::vector<int> component, int no_of_components);
		void cluster(Graph &graph);
		void normalized_cut(Graph &graph, int max_size, int eigen_vector_no);	//05-25-05
		void old_run(dict graph_dict);
		void run();	//05-25-05
		void walk_graph(Graph &subgraph, Graph &graph, vertexNamePropertyMap vertex2name);	//05-25-05
		void output();	//05-25-05
		
		Graph g;
		//for input purpose
		std::map<int, vertexDescriptor> geneNoMap;
		//property_map of the vertex names
		std::map<std::string, vertexDescriptor> VertexNameMap;

		vertexNamePropertyMap vertex2name;
		int min_cluster_size;
		int eigen_vector_no;
		double connectivity_cutoff;
		std::string input_filename;
		std::string output_filename;
		int max_size;
		int cut_loop_num;
		int min_edge_weight;
		//the weight map
		iterator_property_map<int*, EdgeIndexMap, int, int&> weight_pa;
		//store the final results
		list good_clusters;
		//store the final results
		std::vector<Graph> good_clusters_vector;

};

clustering::clustering(double connectivity, int which_eigen_vector, int cluster_size)
{
	//parameter initialization
	eigen_vector_no = which_eigen_vector;		//second minimum
	connectivity_cutoff = connectivity;	
	min_cluster_size = cluster_size;
}

clustering::clustering(std::string input_filename_option, std::string output_filename_option, int max_size_option,
			int cut_loop_num_option, double density_cutoff_option, int min_edge_weight_option)
{
	input_filename = input_filename_option;
	output_filename = output_filename_option;
	max_size = max_size_option;
	cut_loop_num = cut_loop_num_option;
	connectivity_cutoff = density_cutoff_option;
	min_edge_weight = min_edge_weight_option;
	eigen_vector_no = 1;
}

clustering::~clustering()
{
	//do nothing
	#if defined(DEBUG)
		std::cout<<"clustering exits"<<std::endl;
	#endif
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

void clustering::init_graph_from_file(std::string input_filename, Graph &graph, int min_edge_weight)
{
	std::cerr<<"Read in graph from "<<input_filename<<std::endl;
	std::ifstream datafile(input_filename.c_str());
	std::vector<int> weight_array;	//05-25-05	a local weight_array
	vertex2name = get(vertex_name, graph);
	for (std::string line; std::getline(datafile, line);) {
		char_separator<char> sep(" \t");		//05-25-05	blank or '\t' is the separator
		typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
		tokenizer line_toks(line, sep);
		tokenizer::iterator i = line_toks.begin();
		if(*i=="t")
			continue;	//skip the whole line
		if(*i=="v")
			continue;	//skip the whole line
		if(*i=="e")
			*i++;	//skip 'e'
		
		std::string gene1 = *i++;
		std::string gene2 = *i++;
		std::string edge_weight_string = *i;
		int edge_weight = atoi(edge_weight_string.c_str());
		if (edge_weight>min_edge_weight)
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
	std::cerr<<"Done."<<std::endl;
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
/*
*05-26-05
*	transformed into a laplacian matrix
*/
{
	int dimension = num_vertices(graph);
	#if defined(DEBUG)
		std::cout<<"Dimension of the graph is "<<dimension<<std::endl;
	#endif
	gsl_matrix* m = gsl_matrix_calloc(dimension, dimension);	//calloc sets all elements to 0, different from alloc
	boost::property_map<Graph, vertex_index_t>::type
	vertex_id = get(vertex_index, graph);
	vertexDescriptor vertex1, vertex2;
	int index1,index2;
	int degree1, degree2;
	graph_traits<Graph>::edge_iterator ei, ei_end;
	for (tie(ei, ei_end) = edges(graph); ei!=ei_end; ++ei)
	{
		vertex1 = source(*ei, graph);
		vertex2 = target(*ei, graph);
		degree1 = degree(vertex1, graph);
		degree2 = degree(vertex2, graph);
		index1 = get(vertex_id, vertex1);
		index2 = get(vertex_id, vertex2);
		gsl_matrix_set(m, index1, index1, degree1);
		gsl_matrix_set(m, index2, index2, degree2);
		gsl_matrix_set(m, index1, index2, -1.0);
		gsl_matrix_set(m, index2, index1, -1.0);	//undirected, symmetric
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
	#if defined(DEBUG)
		int i;
		std::cout<<"The "<<which_eigen_vector<<"th eigenvalue: "<<gsl_vector_get(eval, which_eigen_vector)<<std::endl;
		std::cout<<"The "<<which_eigen_vector<<"th eigenvector: "<<std::endl;
		for(i=0;i<evec_i->size;++i)
			std::cout<<gsl_vector_get(evec_i, i)<<"\t";
		std::cout<<std::endl;
	#endif
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

void clustering::normalized_cut(Graph &graph, int max_size, int eigen_vector_no)
/*
*05-25-05	normalized_cut, similar to cluster()
*/
{
	int i;
	gsl_matrix* m = graph2gsl_matrix(graph);
	#if defined(DEBUG)
		std::cout<<"The matrix is "<<std::endl;
		gsl_matrix_fprintf(stdout, m, "%g");		//check matrix
	#endif
	gsl_vector* evec_i = return_eigen_vector(m, eigen_vector_no);

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
		//get all the components and check
		std::vector<int> component(num_vertices_of_subgraph);
		int no_of_components = connected_components(vector_subgraph[i], &component[0]);
		//the second parameter is g, not graph
		std::vector<Graph> vector_sub_subgraph = subgraph_components(vector_subgraph[i], g, component, no_of_components);
		
		std::vector<Graph>::iterator g_iterator;
		#if defined(DEBUG)
			std::cout<<"No. of components in subgraph "<<i<<" is: "<<no_of_components<<std::endl;
		#endif
		for(g_iterator=vector_sub_subgraph.begin();g_iterator!=vector_sub_subgraph.end();++g_iterator)
		{
			if(num_vertices(*g_iterator)>max_size)
			{
				normalized_cut(*g_iterator, max_size, eigen_vector_no);
			}
			else
			{
				good_clusters_vector.push_back(*g_iterator);
			}
		}
	}
}


void clustering::walk_graph(Graph &subgraph, Graph &graph, vertexNamePropertyMap vertex2name)
/*
*05-25-05
*	copied from bgl_test.cc, modified, combine walk_edges and walk_vertices
*/
{
	
	boost::property_map<Graph, vertex_index_t>::type
	vertex_id = get(vertex_index, graph);
	boost::graph_traits<Graph>::vertex_descriptor
	vertex_local, vertex_local1, vertex_global, vertex_global1;
	
	std::cout<<"no. of vertices: "<<num_vertices(subgraph)<<std::endl;
	
	typedef graph_traits<Graph>::vertex_iterator vertex_iter;
	std::pair<vertex_iter, vertex_iter> vp;
	std::cout << "vertices(g) = ";
	for (vp = vertices(subgraph); vp.first != vp.second; ++vp.first)
	{
		vertex_global = subgraph.local_to_global(*vp.first);
		std::cout << get(vertex2name, vertex_global) <<  " ";
	}
	std::cout << std::endl;
	
	std::cout << "no of edges: "<<num_edges(subgraph)<<std::endl;
	graph_traits<Graph>::edge_iterator ei, ei_end;
	std::cout << "edges(g)= ";
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
		std::cout << "[" << get(vertex2name, vertex_global)
		<< "," << get(vertex2name, vertex_global1) << "] ";
	}
	std::cout << std::endl<<std::endl;
	
}

void clustering::output()
/*
*05-25-05
*	copied from bgl_test.cc
*/
{
	std::cerr<<"Outputting subgraphs..."<<std::endl;
	std::vector<Graph>::iterator g_iterator;
	for(g_iterator=good_clusters_vector.begin();g_iterator!=good_clusters_vector.end();++g_iterator)
		walk_graph(*g_iterator, g, vertex2name);
	std::cerr<<"Done"<<std::endl;
}

void clustering::old_run(dict graph_dict)
{
	init_graph_from_dict(graph_dict, g);
	cluster(g);
}

void clustering::run()
{
	init_graph_from_file(input_filename, g, min_edge_weight);
	
	//get all the components and check
	std::vector<int> component(num_vertices(g));
	int no_of_components = connected_components(g, &component[0]);
	//the second parameter is g, not graph
	std::vector<Graph> vector_sub_subgraph = subgraph_components(g, g, component, no_of_components);
	
	std::vector<Graph>::iterator g_iterator;
	#if defined(DEBUG)
		std::cout<<"No. of components in the g is: "<<no_of_components<<std::endl;
	#endif
	for(g_iterator=vector_sub_subgraph.begin();g_iterator!=vector_sub_subgraph.end();++g_iterator)
	{
		if(num_vertices(*g_iterator)>max_size)
		{
			normalized_cut(*g_iterator, max_size, eigen_vector_no);
		}
		else
		{
			good_clusters_vector.push_back(*g_iterator);
		}
	}
	output();
}

BOOST_PYTHON_MODULE(clustering)
{
	class_<clustering>("clustering", init<double, int, int>())
		.def("init_graph_from_dict", &clustering::init_graph_from_dict)
		.def("run", &clustering::run)
		.def_readonly("good_clusters", &clustering::good_clusters)
	;
}

void print_usage(FILE* stream, char* program_name)
{
	assert(stream !=NULL);
        fprintf(stream,"Usage: %s options -i INPUTFILE\n",program_name);
	fprintf(stream,"\t-h  --help	Display the usage infomation.\n"\
		"\t-i ..., --input=...	INPUTFILE\n"\
		"\t-o ..., --output=...	Write output to file, INPUTFILE.1st(default)\n"\
		"\t-s ..., --max_size=...	Min graph size, 200(default)\n"\
		"\t-r .., --cut_loop_num=...	cut_loop_num, 2(default)\n"\
		"\t-d ..., --density_cutoff=...	density cutoff, 0.2(default)\n"\
		"\t\tif max_size=0, this cut_off is used instead.\n"\
		"\t-e ..., --min_edge_weight=...	minimum edge weight, 5(default).\n"\
		"\tFor long option, = or ' '(blank) is same.\n");
	exit(3);
}


int main(int argc, char* argv[])
{
	int next_option;
	const char* const short_options="hi:o:s:r:d:e:";
	const struct option long_options[]={
	  {"help",0,NULL,'h'},
	  {"input", 1, NULL, 'i'},
	  {"output",1,NULL,'o'},
	  {"max_size",1,NULL,'s'},
	  {"cut_loop_num", 1, NULL, 'r'},
	  {"density_cutoff", 1, NULL, 'd'},
	  {"min_edge_weight", 1, NULL, 'e'},
	  {NULL,0,NULL,0}
	};
	
	char* program_name=argv[0];	
	std::string input_filename = "";
	std::string output_filename = "";
	int max_size = 200;
	int cut_loop_num = 2;
	double density_cutoff = 0.2;
	int min_edge_weight = 5;

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
			max_size = atoi(optarg);
			break;
		case 'r':
			cut_loop_num = atoi(optarg);
			break;
		case 'd':
			density_cutoff = atof(optarg);
			break;
		case 'e':
			min_edge_weight = atoi(optarg);
			break;
		case '?':
			print_usage(stderr, program_name);
		case -1:
			break;
		default:
			abort();
		}
	}while(next_option!=-1);
	
	if (output_filename == "")
	{
		output_filename = input_filename+".1st";
	}
	//ifstream inf(argv[1]);
	//ofstream outf(argv[2], ios::app | ios::out);

	if (input_filename!="")
	{
		clustering instance(input_filename, output_filename, max_size, cut_loop_num, density_cutoff, min_edge_weight);
		instance.run();
	}
	else
		print_usage(stderr, program_name);
}
