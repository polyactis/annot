/*
*09-04-05
*	class definition and included header files go to clustering.h
*/
#include "clustering.h"

clustering::clustering()
{
}

clustering::clustering(double connectivity, int which_eigen_vector, int cluster_size)
{
	//parameter initialization
	eigen_vector_no = which_eigen_vector;		//second minimum
	connectivity_cutoff = connectivity;	
	min_cluster_size = cluster_size;
}

clustering::clustering(std::string input_filename_option, std::string output_filename_option, int max_size_option,
			int cut_loop_num_option, double density_cutoff_option, int min_edge_weight_option, int matrix_format_option)
{
	input_filename = input_filename_option;
	output_filename = output_filename_option;
	max_size = max_size_option;
	cut_loop_num = cut_loop_num_option;
	connectivity_cutoff = density_cutoff_option;
	min_edge_weight = min_edge_weight_option;
	matrix_format = matrix_format_option;
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
	std::cerr<<"Read in graph from "<<input_filename<<" ...";
	std::ifstream datafile(input_filename.c_str());
	std::vector<int> weight_array;	//05-25-05	a local weight_array
	vertex2name = get(vertex_name, graph);
	for (std::string line; std::getline(datafile, line);) {
		char_separator<char> sep(" \t");		//05-25-05	blank or '\t' is the separator
		typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
		tokenizer line_toks(line, sep);
		tokenizer::iterator tokenizer_iter = line_toks.begin();
		if(*tokenizer_iter=="t")
			continue;	//skip the whole line
		if(*tokenizer_iter=="v")
			continue;	//skip the whole line
		if(*tokenizer_iter=="e")
			*tokenizer_iter++;	//skip 'e'
		
		std::string gene1 = *tokenizer_iter++;
		std::string gene2 = *tokenizer_iter++;
		std::string edge_weight_string = *tokenizer_iter;
		int edge_weight = atoi(edge_weight_string.c_str());
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
	std::cerr<<"Done."<<std::endl;
}

void clustering::init_graph_from_file_matrix(std::string input_filename, Graph &graph, int min_edge_weight)
/*
*05-26-05
*	add this function to read in matrix format input file.
*/
{
	std::cerr<<"Read in graph from matrix_file "<<input_filename<<"...";
	std::ifstream datafile(input_filename.c_str());
	std::vector<int> weight_array;	//05-25-05	a local weight_array
	vertex2name = get(vertex_name, graph);
	int i=0;
	int j=0;
	for (std::string line; std::getline(datafile, line);) {
		char_separator<char> sep(" \t");		//05-25-05	blank or '\t' is the separator
		typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
		tokenizer line_toks(line, sep);
		j=0;
		for (tokenizer::iterator tokenizer_iter = line_toks.begin(); tokenizer_iter!=line_toks.end();++tokenizer_iter)
		{
			if(i<j)	//05-26-05 only triangle
			{
				int edge_weight = atoi((*tokenizer_iter).c_str());
				if(edge_weight>=min_edge_weight)
				{
					std::map<int, vertexDescriptor>::iterator pos;
					bool inserted;
					vertexDescriptor u, v;
					tie(pos, inserted) = geneNoMap.insert(std::make_pair(i, vertexDescriptor()));
					if (inserted) {
						u = add_vertex(graph);
						vertex2name[u] = i;
						pos->second = u;
					} else
						u = pos->second;
					
					tie(pos, inserted) = geneNoMap.insert(std::make_pair(j, vertexDescriptor()));
					if (inserted) {
						v = add_vertex(graph);
						vertex2name[v] = j;
						pos->second = v;
					} else
						v = pos->second;
					
					graph_traits < Graph >::edge_descriptor e;
					tie(e, inserted) = add_edge(u, v, graph);
					if (inserted)
						weight_array.push_back(edge_weight);
				}
			}
			
			j++;
		}
		
		i++;
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
*05-26-05
*	transformed into a normalized laplacian matrix
*/
{
	int dimension = num_vertices(graph);
	#if defined(DEBUG)
		std::cerr<<"Dimension of the graph is "<<dimension<<std::endl;
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
		gsl_matrix_set(m, index1, index1, 1.0);
		gsl_matrix_set(m, index2, index2, 1.0);
		gsl_matrix_set(m, index1, index2, -1.0*pow(degree1,-0.5)*pow(degree2,-0.5));
		gsl_matrix_set(m, index2, index1, -1.0*pow(degree1,-0.5)*pow(degree2,-0.5));	//undirected, symmetric
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

gsl_vector* clustering::return_eigen_vector_by_arpack(Graph &graph, int which_eigen_vector)
/*
*05-30-05
*	use AREig() to get eigenvector for large matrix fastly.
*	see areig.h for more information
*	it's a laplacian matrix, different from bgl_test.cc
*/
{
	int n;
	int nnz;
	int* irow;
	int* pcol;
	double* A;
	
	n = num_vertices(graph);
	nnz = num_edges(graph)+n;
	A = new double[nnz];
	irow = new int[nnz];
	pcol = new int[n+1];
	
	std::vector<vertexDescriptor> adjacency_vector;
	boost::graph_traits<Graph>::adjacency_iterator a_it, a_it_end;
	boost::graph_traits<Graph>::vertex_descriptor s, t;
	std::vector<vertexDescriptor>::iterator av_it;
	int j=0;
	int i=0;
	pcol[0] = 0;	//first one is always the first element of the first column.
	int degree1, degree2;
	for(i=0; i<n; i++)
	{
		s = vertex(i, graph);
		adjacency_vector.clear();
		for(tie(a_it, a_it_end)=adjacent_vertices(s, graph); a_it!=a_it_end; ++a_it)
		{
			adjacency_vector.push_back(*a_it);
		}
		std::sort(adjacency_vector.begin(), adjacency_vector.end());
		#if defined(DEBUG)
			std::cout<<"The sorted adjacency list of vertex "<<s<<" is: ";
			std::copy(adjacency_vector.begin(), adjacency_vector.end(), std::ostream_iterator<vertexDescriptor>(std::cout, " "));
			std::cout<<std::endl;
		#endif
		
		//the diagonal element
		irow[j] = i;
		A[j++] = 1;
		
		av_it = adjacency_vector.begin();
		for(av_it = adjacency_vector.begin(); av_it!=adjacency_vector.end(); ++av_it)
		{
			if(*av_it>i)
			{
				//if(*av_it<pcol[i])
				//	pcol[i] = int(*av_it);
				t = vertex(*av_it, graph);
				degree1 = degree(s, graph);
				degree2 = degree(t, graph);
				irow[j] = int(*av_it);
				A[j++] = -1.0*pow(degree1,-0.5)*pow(degree2,-0.5);	//normalized laplacian
			}
		}
		pcol[i+1] = j;
		
	}
	#if defined(DEBUG)
		std::cout<<"The number of non-zero elements in A is "<<j<<"."<<std::endl;
		std::copy(A, A+nnz, std::ostream_iterator<double>(std::cout, " "));
		std::cout<<std::endl;
		std::cout<<"Elements of irow is ";
		std::copy(irow, irow+nnz, std::ostream_iterator<int>(std::cout, " "));
		std::cout<<std::endl;
		std::cout<<"Elements of pcol is ";
		std::copy(pcol, pcol+n+1, std::ostream_iterator<int>(std::cout, " "));
		std::cout<<std::endl;
	#endif
	
	//prepare to do AREig()
	double* EigVal;
	double* EigVec;
	EigVal = new double[n+1];
	EigVec = new double[n*(which_eigen_vector+1)*2+1];	//including the number of eigenvectors before which_eigen_vector, 
		//It's given in complex eigenvectors with two consecutive vectors.
	char uplo='L';
	int nev = which_eigen_vector+1;	//which_eigen_vector starts from 0.
	char* which="SM";
	int ncv = nev*2+1;
	double tol = 0.0001;	//Stopping criterion (relative accuracy of Ritz values), 0.01 or 0.00001 doesn't matter.
	int maxit = 100000;
	
	int nconv = AREig(EigVal, EigVec, n, nnz, A, irow, pcol, uplo, nev, which, ncv, tol, maxit);
	gsl_vector* evec_i = gsl_vector_alloc(n);
	for(i=0; i<n; i++)
		gsl_vector_set(evec_i, i, EigVec[n*which_eigen_vector+i]);	//it's in value(not absolute) ascending order).
	
	#if defined(DEBUG)
		std::cout<<"No of converged eigen values is "<<nconv<<std::endl;
		std::copy(EigVal, EigVal+nconv, std::ostream_iterator<double>(std::cout, " "));
		std::cout<<std::endl;
		//output the eigen vector
		std::cout<<"The "<<which_eigen_vector<<"th eigenvector: "<<std::endl;
		for(i=0;i<evec_i->size;++i)
			std::cout<<gsl_vector_get(evec_i, i)<<"\t";
		std::cout<<std::endl;
	#endif

	delete[] EigVal;
	delete[] EigVec;
	delete[] A;
	delete[] irow;
	delete[] pcol;
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


void clustering::normalized_cut(std::ofstream &outf, Graph &subgraph, Graph &graph, int max_size, int eigen_vector_no)
/*
*05-25-05	normalized_cut, similar to cluster()
*05-26-05	fix an important bug
*05-30-05	replace gsl functions with return_eigen_vector_by_arpack()
*/
{
	int i, global_index;
	/*
	gsl_matrix* m = graph2gsl_matrix(subgraph);	//05-26-05	transform subgraph to matrix
	#if defined(DEBUG)
		std::cerr<<"The matrix is "<<std::endl;
		gsl_matrix_fprintf(stdout, m, "%g");		//check matrix
	#endif
	gsl_vector* evec_i = return_eigen_vector(m, eigen_vector_no);
	*/
	gsl_vector* evec_i = return_eigen_vector_by_arpack(subgraph, eigen_vector_no);	
	//two vertex_descriptor
	boost::graph_traits<Graph>::vertex_descriptor
	vertex_local, vertex_global;
	//the vertex_index_global map is used to translate the global descriptor to the global index.
	boost::property_map<Graph, vertex_index_t>::type
	vertex_index_global;
	vertex_index_global = get(vertex_index, graph);
	
	//split the big graph based on eigenvector, >0 or <0
	std::vector<Graph> vector_subgraph(2);
	vector_subgraph[0] = graph.create_subgraph();	//05-26-05 always use the ancestor graph to create subgraph
	vector_subgraph[1] = graph.create_subgraph();
	
	for(i=0;i<evec_i->size;++i)
	{
		if (gsl_vector_get(evec_i, i)<0)
		{
			//i is the local index, get a descriptor from it
			vertex_local = vertex(i, subgraph);
			//find the global descriptor
			vertex_global = subgraph.local_to_global(vertex_local);
			//get the global index and add it to the subgraph
			global_index = get(vertex_index_global, vertex_global);
			#if defined(DEBUG)
				std::cerr<<global_index<<" goes to the first subgraph"<<std::endl;
			#endif
			add_vertex(global_index, vector_subgraph[0]);
		}
		else
		{
			//i is the local index, get a descriptor from it
			vertex_local = vertex(i, subgraph);
			//find the global descriptor
			vertex_global = subgraph.local_to_global(vertex_local);
			//get the global index and add it to the subgraph
			global_index = get(vertex_index_global, vertex_global);
			#if defined(DEBUG)
				std::cerr<<global_index<<" goes to the second subgraph"<<std::endl;
			#endif
			add_vertex(global_index, vector_subgraph[1]);
		}
	}
	
	for(i=0; i<2; i++)
	{
		int num_vertices_of_subgraph = num_vertices(vector_subgraph[i]);
		//get all the components and check
		std::vector<int> component(num_vertices_of_subgraph);
		int no_of_components = connected_components(vector_subgraph[i], &component[0]);
		//the second parameter is g, not graph
		std::vector<Graph> vector_sub_subgraph = subgraph_components(vector_subgraph[i], graph, component, no_of_components);
		
		std::vector<Graph>::iterator g_iterator;
		#if defined(DEBUG)
			std::cerr<<"No. of components in subgraph "<<i<<" is: "<<no_of_components<<std::endl;
		#endif
		for(g_iterator=vector_sub_subgraph.begin();g_iterator!=vector_sub_subgraph.end();++g_iterator)
		{
			if(num_vertices(*g_iterator)>max_size)
			{
				normalized_cut(outf, *g_iterator, graph, max_size, eigen_vector_no);
			}
			else
			{
				walk_graph(outf, *g_iterator, graph);
				//good_clusters_vector.push_back(*g_iterator);
			}
		}
	}
}


void clustering::walk_graph(std::ofstream &outf, Graph &subgraph, Graph &graph)
/*
*05-25-05
*	copied from bgl_test.cc, modified, combine walk_edges and walk_vertices
*/
{
	#if defined(DEBUG)
		std::cerr<<"Outputting subgraph...";
	#endif
	
	boost::property_map<Graph, vertex_index_t>::type
	vertex_id = get(vertex_index, graph);
	vertexNamePropertyMap vertex2name = get(vertex_name, graph);
	boost::graph_traits<Graph>::vertex_descriptor
	vertex_local, vertex_local1, vertex_global, vertex_global1;
	
	outf<<"no. of vertices: "<<num_vertices(subgraph)<<std::endl;
	
	typedef graph_traits<Graph>::vertex_iterator vertex_iter;
	std::pair<vertex_iter, vertex_iter> vp;
	outf << "vertices(g) = ";
	for (vp = vertices(subgraph); vp.first != vp.second; ++vp.first)
	{
		vertex_global = subgraph.local_to_global(*vp.first);
		outf << get(vertex2name, vertex_global) <<  " ";
	}
	outf << std::endl;
	
	outf << "no of edges: "<<num_edges(subgraph)<<std::endl;
	graph_traits<Graph>::edge_iterator ei, ei_end;
	outf << "edges(g)= ";
	for (tie(ei,ei_end) = edges(subgraph); ei != ei_end; ++ei)
	{
		vertex_local = source(*ei, subgraph);
		vertex_local1 = target(*ei, subgraph);
		vertex_global = subgraph.local_to_global(vertex_local);
		vertex_global1 = subgraph.local_to_global(vertex_local1);
		#if defined(DEBUG)
			outf << "(" << get(vertex_id, vertex_global)
			<< "," << get(vertex_id, vertex_global1) << ") ";
		#endif
		outf << "[" << get(vertex2name, vertex_global)
		<< "," << get(vertex2name, vertex_global1) << "] ";
	}
	outf << std::endl<<std::endl;
	
	#if defined(DEBUG)
		std::cerr<<"Done."<<std::endl;
	#endif
}

void clustering::output(std::string output_filename)
/*
*05-25-05
*	copied from bgl_test.cc
*/
{
	std::cerr<<"Outputting subgraphs..."<<std::endl;
	std::ofstream outf(output_filename.c_str());
	std::vector<Graph>::iterator g_iterator;
	for(g_iterator=good_clusters_vector.begin();g_iterator!=good_clusters_vector.end();++g_iterator)
		walk_graph(outf, *g_iterator, g);
	outf.close();
	std::cerr<<"Done"<<std::endl;
}

void clustering::old_run(dict graph_dict)
{
	init_graph_from_dict(graph_dict, g);
	//cluster(g);
}

void clustering::run()
{
	if(matrix_format)
		init_graph_from_file_matrix(input_filename, g, min_edge_weight);
	else
		init_graph_from_file(input_filename, g, min_edge_weight);
	std::ofstream outf(output_filename.c_str());
	eigen_vector_no = 1;
	//normalized_cut(outf, g, g, max_size, eigen_vector_no);
	
	//05-26-05 g itself might have several connected components
	//get all the components and check
	std::vector<int> component(num_vertices(g));
	int no_of_components = connected_components(g, &component[0]);
	//the second parameter is g, not graph
	std::vector<Graph> vector_sub_subgraph = subgraph_components(g, g, component, no_of_components);
	
	std::vector<Graph>::iterator g_iterator;
	#if defined(DEBUG)
		std::cerr<<"No. of components in the g is: "<<no_of_components<<std::endl;
	#endif
	for(g_iterator=vector_sub_subgraph.begin();g_iterator!=vector_sub_subgraph.end();++g_iterator)
	{
		if(num_vertices(*g_iterator)>max_size)
		{
			normalized_cut(outf, *g_iterator, g, max_size, eigen_vector_no);
		}
		else
		{
			walk_graph(outf, *g_iterator, g);
			//good_clusters_vector.push_back(*g_iterator);
		}
	}
	
	outf.close();
	//output(output_filename);
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
		"\t-i ..., --input=...	INPUTFILE(gspan or haiyan's edge format)\n"\
		"\t-o ..., --output=...	Write output to file, INPUTFILE.1st(default)\n"\
		"\t-s ..., --max_size=...	Min graph size, 200(default)\n"\
		"\t-r .., --cut_loop_num=...	cut_loop_num, 2(default)\n"\
		"\t-d ..., --density_cutoff=...	density cutoff, 0.2(default)\n"\
		"\t\tif max_size=0, this cut_off is used instead.\n"\
		"\t-e ..., --min_edge_weight=...	minimum edge weight, 5(default).\n"\
		"\t-m, --matrix_format	the inputfile is in matrix format.\n"\
		"\tFor long option, = or ' '(blank) is same.\n");
	exit(3);
}


int main(int argc, char* argv[])
{
	int next_option;
	const char* const short_options="hi:o:s:r:d:e:m";
	const struct option long_options[]={
	  {"help",0,NULL,'h'},
	  {"input", 1, NULL, 'i'},
	  {"output",1,NULL,'o'},
	  {"max_size",1,NULL,'s'},
	  {"cut_loop_num", 1, NULL, 'r'},
	  {"density_cutoff", 1, NULL, 'd'},
	  {"min_edge_weight", 1, NULL, 'e'},
	  {"matrix_format",0,NULL,'m'},
	  {NULL,0,NULL,0}
	};
	
	char* program_name=argv[0];	
	std::string input_filename = "";
	std::string output_filename = "";
	int max_size = 200;
	int cut_loop_num = 2;
	double density_cutoff = 0.2;
	int min_edge_weight = 5;
	int matrix_format = 0;

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
		case 'm':
			matrix_format = 1;
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
		clustering instance(input_filename, output_filename, max_size, cut_loop_num, density_cutoff, min_edge_weight, matrix_format);
		instance.run();
	}
	else
		print_usage(stderr, program_name);
}
