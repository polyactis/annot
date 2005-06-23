/*
*
*03-02-05
*       add a feature to get a significance flag for each edge
*/


#include "graph_modeling.h"

const int GENE_CUT_OFF = 8;
const int JK_CUT_OFF = 7;
const float COR_CUT_OFF = 0.6;

vector<float> cor_cut_off_vector;

void python_call(char* outf_name, vector<int> edge_vector, vector<string> expr_array, int no_genes)
{
	graph_construct instance(outf_name, edge_vector);
	instance.gene_array_fill(expr_array, no_genes);
	instance.edge_construct_no_cut_off();
}

edge ind_min_cor(vf v1, vf v2)
{
	edge edge_data_to_return;
	edge edge_data_tmp;
	edge_data_to_return.value = 1.1;
	edge_data_to_return.degree = 0;
	//03-02-05	significance flag default to 0
	edge_data_to_return.significance = 0;
	for(int i=0; i<v1.size(); i++)
	{
		edge_data_tmp = ind_cor(v1, v2, i);
		if(abs(edge_data_tmp.value)<abs(edge_data_to_return.value) && (edge_data_tmp.degree+2)>=JK_CUT_OFF)
			//04-02-05 if no_of_valids < JK_CUT_OFF, should not modify the value to return
			edge_data_to_return = edge_data_tmp;
	}
	if((edge_data_to_return.degree+2)>=JK_CUT_OFF && \
		edge_data_to_return.value >= cor_cut_off_vector[edge_data_to_return.degree-1] && edge_data_to_return.value<=1.0)
	{
		//03-02-05	valid data points = degree + 2
		edge_data_to_return.significance = 1;
	}
	return edge_data_to_return;
}

edge ind_cor(vf v1, vf v2, int position)
{
	/*
	*05-12-05
	*	add gsl_isnan() to check whether it's NAN, for graph_construct class usage
	*/
	float xx = 0.0;
	float yy = 0.0;
	float xy = 0.0;
	float mean_x = 0.0;
	float mean_y = 0.0;
	int no_of_valids=0;
	edge edge_data;

	for (int i=0; i<v1.size(); i++)
	{
		//100000000 is regarded as NAN
		if ( (i == position) || ( v1[i]==100000000 )|| ( v2[i]==100000000 )|| ( gsl_isnan(v1[i]) )|| ( gsl_isnan(v2[i]) ) )
			continue;
		else
		{
			mean_x += v1[i];
			mean_y += v2[i];
			no_of_valids++;
		}
	}
	
	mean_x /= no_of_valids;
	mean_y /= no_of_valids;
	
	for (int i=0; i<v1.size(); i++)
	{
		if ( (i == position) || ( v1[i]==100000000 )|| ( v2[i]==100000000) || ( gsl_isnan(v1[i]) )|| ( gsl_isnan(v2[i]) ) )
			continue;
		else
		{
			xy += (v1[i]-mean_x) * (v2[i]-mean_y);
			xx += (v1[i]-mean_x) * (v1[i]-mean_x) ;
			yy += (v2[i]-mean_y)  * (v2[i]-mean_y) ;
		}
	}
	
	if(xx==0.0 || yy==0.0)
		//all NAN
		edge_data.value = 1.1;
	else
		edge_data.value = xy/(sqrt(xx*yy));
	//correlation is modeled as a t distribution of n-2 degree.
	edge_data.degree = no_of_valids-2;
	edge_data.significance = 0;
	return edge_data;
}

//04-25-05	calculate the euclidean distance
edge euc_dist(vf v1, vf v2)
{
	float xx = 0.0;
	int no_of_valids=0;
	edge edge_data;

	for (int i=0; i<v1.size(); i++)
	{
		//100000000 is regarded as NAN
		if ( ( v1[i]==100000000 )|| ( v2[i]==100000000 ) )
			continue;
		else
		{
			xx += (v1[i] - v2[i])*(v1[i]-v2[i]);
			no_of_valids++;
		}
	}
	
	if(no_of_valids==0)
		//all NAN
		edge_data.value = -1;
	else
		edge_data.value = sqrt(xx);
	//correlation is modeled as a t distribution of n-2 degree.
	edge_data.degree = no_of_valids-2;
	edge_data.significance = 0;
	return edge_data;
}

//03-02-05	fill the global cor_cut_off_vector with cutoff values, if p_value_cut_off is ==0, use cor_cut_off_given instead.
void cor_cut_off_vector_construct(double p_value_cut_off, double cor_cut_off_given)
{
	double cor_cut_off, t;
	for(int i=1; i<400; i++)
	{
		if(p_value_cut_off != 0)
		{
			t = gsl_cdf_tdist_Qinv(p_value_cut_off, i);
			//convert the t to the correlation, see log file.
			cor_cut_off = sqrt(gsl_pow_2(t)/(gsl_pow_2(t)+i));
			cor_cut_off_vector.push_back((float)cor_cut_off);
		}
		else
			cor_cut_off_vector.push_back((float)cor_cut_off_given);

	}
}

//05-27-05
vector<float> cor_cut_off_vector_return(double p_value_cut_off, double cor_cut_off_given)
{
	cor_cut_off_vector_construct(p_value_cut_off, cor_cut_off_given);
	return cor_cut_off_vector;
}
	

graph_construct::graph_construct(char* outf_name, vector<int> edge_vector)
{
	out.open(outf_name);
	edge_tuple_vector = edge_vector;
	no_of_01 = 0;
}

graph_construct::graph_construct(char* inf_name, char* outf_name, char* g_name, double p_value_cut_off_given, \
	double cor_cut_off_given, float top_percentage_given, int max_degree_given, bool leave_one_out_given)
{
	in.open(inf_name);
	out.open(outf_name);
	graph_name = g_name;
	p_value_cut_off = p_value_cut_off_given;
	cor_cut_off = cor_cut_off_given;
	top_percentage = top_percentage_given;
	max_degree = max_degree_given;
	leave_one_out = leave_one_out_given;
	//ios::app | ios::out);
	no_of_01 = 0;
	//histogram = gsl_histogram_alloc (50);
	//gsl_histogram_set_ranges_uniform (histogram, -5.0, 5.0);
}

graph_construct::~graph_construct()
{
	in.close();
	out.close();
	//gsl_histogram_free (histogram);
}

int graph_construct::input(float top_percentage)
{	
	#if defined(DEBUG)
		std::cerr<<"Read in the data...";
	#endif
	string line;
	while(getline(in, line))
	{
		split(line);
	}
	no_of_genes = gene_array.size();
	no_of_cols = gene_array[0].size();
	#if defined(DEBUG)
		std::cerr<<no_of_genes<<" genes."<<endl;
	#endif
	return int(no_of_genes*no_of_genes*top_percentage);
}

vector<float> graph_construct::cor_cut_off_array_construct(double p_value_cut_off, double cor_cut_off_given, int max_degree)
{
	/*
	*04-30-05
	*	max_degree to control the length of the cor_cut_off_array
	*05-12-05
	*	return cor_cut_off_array for global cor_cut_off_vector
	*/
	#if defined(DEBUG)
		std::cerr<<"constructing cor_cut_off_array...";
	#endif
	double cor_cut_off, t;
	for(int i=1; i<max_degree; i++)
	{
		if(p_value_cut_off != 0)
		{
			//cor_vector.clear();
			t = gsl_cdf_tdist_Qinv(p_value_cut_off, i);
			//convert the t to the correlation, see log file.
			cor_cut_off = sqrt(gsl_pow_2(t)/(gsl_pow_2(t)+i));
			cor_cut_off_array.push_back((float)cor_cut_off);
		}
		else	
			cor_cut_off_array.push_back((float)cor_cut_off_given);
		
		/*old method
		cor_list = general_split(line, '\t');
		for(i=0; i<cor_list.size(); i++)
			cor_vector.push_back(atof(cor_list[i].c_str()));
		cor_cut_off_array.push_back(cor_vector);
		*/

	}
	#if defined(DEBUG)
		std::cerr<<"Done."<<std::endl;
	#endif
	return cor_cut_off_array;
}

void graph_construct::gene_label2index_setup(vector<string> label_vector)
{
	//typedef pair <std::string, int> string_int_Pair;
	no_of_genes = label_vector.size();
	for(int i=0; i<label_vector.size();i++)
	{
		//gene_label2index.insert(string_int_Pair(label_vector[i], i));
		gene_label2index[label_vector[i]] = i;
		cout<<label_vector[i]<<'\t'<<gene_label2index.size()<<endl;
	}
}

void graph_construct::gene_array_fill(vector<string> expr_array, int no_genes)
{
	vf gene_vector;
	no_of_cols = expr_array.size()/no_genes;
	no_of_genes = no_genes;
	//initialize the gene_vector
	for(int i=0; i<no_of_cols; i++)
		gene_vector.push_back(GSL_NAN);
	//initialize the gene_array
	for(int i=0; i<no_of_genes; i++)
	{
		gene_array.push_back(gene_vector);
	}
	int x_dim = 0;
	int y_dim = 0;
	for(int i=0; i<expr_array.size();i++)
	{
		x_dim = i/no_of_cols;
		y_dim = (int)fmod((float)i, (float)no_of_cols);
		if(expr_array[i] != "NA")
			gene_array[x_dim][y_dim] = atof(expr_array[i].c_str());
	}
}

void graph_construct::split(string line)
{
	vf gene_vector;
	string gene_label;
	bit_vector bv;

	int no_of_tabs = 0;
	int no_of_nas = 0;
	string tmp = "";

	for(int i=0; i< line.size(); i++)
	{
		if (line[i]=='\t')
		{
			no_of_tabs++;
			if (no_of_tabs==1)
				gene_label = tmp;
			else
			{
				if (tmp=="NA")
				{
					gene_vector.push_back(GSL_NAN);
					bv.push_back(true);
					no_of_nas++;
				}
				else
				{
					gene_vector.push_back(atof(tmp.c_str()) );
					bv.push_back(false);
				}
			}
			tmp="";
		}
		else
			tmp+=line[i];
	}
	//the last field
	if (tmp=="NA")
	{
		gene_vector.push_back(GSL_NAN);
		bv.push_back(true);
		no_of_nas++;
	}
	else
		if (tmp!="")  //a row is ended with 0.432\t\n.
		{
			gene_vector.push_back(atof(tmp.c_str())  );
			bv.push_back(false);
		}
	
	//discard some genes with too many missing values.

	if ( gene_vector.size() -  no_of_nas  >=GENE_CUT_OFF)
	{
		gene_labels_vector.push_back(gene_label);
		mask_vector.push_back(bv);
		gene_array.push_back(gene_vector);
	}
}

void graph_construct::edge_construct(bool leave_one_out, int top_number)
{
	/*04-30-05
	*	flag leave_one_out to control leave_one_out or not
	*05-12-05
	*	use ind_min_cor(), as min_cor() is outdated. And ind_min_cor() is modified to be closest to graph.cc.
	*06-22-05
	*	top_number is used to filter edges,
	*/
	#if defined(DEBUG)
		std::cerr<<"Constructing edges...";
	#endif
	edge edge_data;
	edge_string_cor top_element;
	out<<"t\t#\t"<<graph_name<<endl;
	for (int i=0; i<no_of_genes; i++)
	{
		for (int j=i+1; j<no_of_genes; j++)
		{
			if (leave_one_out)
				edge_data = ind_min_cor(gene_array[i], gene_array[j]);
			else
				edge_data = cor(gene_array[i], gene_array[j], -1);	//leave_one_out position=-1 means no leave_one_out
			if ((edge_data.degree+2)>=JK_CUT_OFF && edge_data.value<=1.0)	//enough valid pairs and value is under 1.0
			{
				if (top_number<=0)	//06-22-05	top_number<=0 means no top selection
				{
					if(edge_data.value >= cor_cut_off_array[edge_data.degree-1])
					{
						no_of_01++;
						out<<"e\t"<<gene_labels_vector[i]<<'\t'<<gene_labels_vector[j]<<'\t'<<edge_data.value<<endl;
					}
				}
				else	//06-22-05	if top_number not zero, we do top selection
				{
					if (edge_pq.size()<top_number)
					{
						no_of_01++;	//this counter is nothing
						edge_pq.push(boost::make_tuple(gene_labels_vector[i], gene_labels_vector[j], edge_data.value));
					}
					else
					{
						top_element = edge_pq.top();
						if (edge_data.value>top_element.get<2>())
						{
							edge_pq.pop();	//throw away the smallest
							edge_pq.push(boost::make_tuple(gene_labels_vector[i], gene_labels_vector[j], edge_data.value));
						}
					}
				}
		}
		}
	}
	#if defined(DEBUG)
		std::cerr<<"Done."<<endl;
	#endif
}

vf graph_construct::edge_construct_no_cut_off()
{
	vf data_vector;
	edge edge_data;
	int gene_one, gene_two;
	int variant_value;
	for(int i=0; i<edge_tuple_vector.size(); i=i+2)
	{
		gene_one = edge_tuple_vector[i];
		gene_two  = edge_tuple_vector[i+1];
		edge_data = ind_min_cor(gene_array[gene_one], gene_array[gene_two]);
		//haiyan's program needs the value*1000 and only the integer part
		variant_value = (short)(edge_data.value*1000);
		data_vector.push_back(variant_value);
		out<<variant_value<<endl;
	}
	return data_vector;
}

edge graph_construct::min_cor(vf v1, vf v2)
{
	edge edge_data_to_return;
	edge edge_data_tmp;
	edge_data_to_return.value = 10.0;
	edge_data_to_return.degree = 0;
	for(int i=0; i<v1.size(); i++)
	{
		edge_data_tmp = cor(v1, v2, i);
		if(abs(edge_data_tmp.value)<abs(edge_data_to_return.value))
			edge_data_to_return = edge_data_tmp;
	}
	return edge_data_to_return;
}

edge graph_construct::d_distance(vf v1, vf v2)
{
	int no_of_valids=0;
	edge edge_data;
	double* d_distance_vector = (double *)calloc(v1.size(), sizeof(double));
	for (int i=0; i<v1.size(); i++)
	{
		if ( ( gsl_isnan(v1[i]) )|| ( gsl_isnan(v2[i]) ) )
			continue;
		else
		{
			d_distance_vector[i] = v1[i] - v2[i];
			no_of_valids++;
		}
	}
	double var = gsl_stats_variance(d_distance_vector, 1, no_of_valids)/gsl_pow_2(gsl_stats_mean(d_distance_vector, 1, no_of_valids))	;
	edge_data.value = (float)var;
	//s^2 is of n-1 degree 
	edge_data.degree = no_of_valids-1;
	//gsl_histogram_increment(histogram, var);
	free(d_distance_vector);
	return edge_data;
}

edge graph_construct::cor(vf v1, vf v2, int position)
{
	float xx = 0.0;
	float yy = 0.0;
	float xy = 0.0;
	float mean_x = 0.0;
	float mean_y = 0.0;
	int no_of_valids=0;
	edge edge_data;

	for (int i=0; i<v1.size(); i++)
	{
		if ( (i == position) || ( gsl_isnan(v1[i]) )|| ( gsl_isnan(v2[i]) ) )
			continue;
		else
		{
			mean_x += v1[i];
			mean_y += v2[i];
			no_of_valids++;
		}
	}
	
	mean_x /= no_of_valids;
	mean_y /= no_of_valids;
	
	for (int i=0; i<v1.size(); i++)
	{
		if ( (i == position) || ( gsl_isnan(v1[i]) )|| ( gsl_isnan(v2[i]) ) )
			continue;
		else
		{
			xy += (v1[i]-mean_x) * (v2[i]-mean_y);
			xx += (v1[i]-mean_x) * (v1[i]-mean_x) ;
			yy += (v2[i]-mean_y)  * (v2[i]-mean_y) ;
		}
	}
	
	if(xx==0.0 || yy==0.0)
		//all NAN
		edge_data.value = 1.1;
	else
		edge_data.value = xy/(sqrt(xx*yy));
	//correlation is modeled as a t distribution of n-2 degree.
	edge_data.degree = no_of_valids-2;
	return edge_data;
}

vector<string> graph_construct::general_split(string line, char ch)
{
	vector<string> string_list;
	int delimiter_pos = line.find_first_of(ch);
	while(delimiter_pos!=-1)
	{
		string_list.push_back(line.substr(0, delimiter_pos));
		line.erase(0, (delimiter_pos+1));
		delimiter_pos = line.find_first_of(ch);
	}
	string_list.push_back(line);
	return string_list;
}


void graph_construct::output()
{
	#if defined(DEBUG)
		std::cerr<<no_of_01<<std::endl;
	#endif
	//gsl_histogram_fprintf (stdout, histogram, "%g", "%g");
	
	//06-22-05	output the edge_pq
	#if defined(DEBUG)
		std::cerr<<"Outputing edge_pq with "<<edge_pq.size()<<" edges ...";
	#endif
	edge_string_cor top_element;
	while(!edge_pq.empty())	//don't use edge_pq.size() as it decreases everytime.
	{
		top_element = edge_pq.top();
		out<<"e\t"<<top_element.get<0>()<<'\t'<<top_element.get<1>()<<'\t'<<top_element.get<2>()<<endl;
		edge_pq.pop();	//throw away the first element
	}
	#if defined(DEBUG)
		std::cerr<<"Done."<<std::endl;
	#endif

}

void graph_construct::run()
{
	cor_cut_off_vector = cor_cut_off_array_construct(p_value_cut_off, cor_cut_off, max_degree);
	int top_number = input(top_percentage);
	int top_number_to_be_passed = 0;
	if (p_value_cut_off==0 && cor_cut_off==0)
	{
		#if defined(DEBUG)
			std::cerr<<"Top number is "<<top_number<<std::endl;
		#endif
		top_number_to_be_passed = top_number;
	}
	edge_construct(leave_one_out, top_number_to_be_passed);	//0 means no top_number selection, non 0 means top_number selection
	output();
}

const char* program_name;

void print_usage(FILE* stream,int exit_code)
{
	assert(stream !=NULL);
        fprintf(stream,"Usage: %s options inputfile\n",program_name);
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
		"\tFor long option, = or ' '(blank) is same.\n");
	exit(3);
}


int main(int argc, char* argv[])
{
	int next_option;
	const char* const short_options="ho:n:p:c:t:d:l";
	const struct option long_options[]={
	  {"help",0,NULL,'h'},
	  {"output",1,NULL,'o'},
	  {"name",1,NULL,'n'},
	  {"p_value_cut_off", 1, NULL, 'p'},
	  {"cor_cut_off", 1, NULL, 'c'},
	  {"top_percentage", 1, NULL, 't'},
	  {"max_degree", 1, NULL, 'd'},
	  {"leave_one_out", 0, NULL, 'l'},
	  {NULL,0,NULL,0}
	};
	
	char* output_filename = "gph_output";
	char* name = NULL;
	program_name=argv[0];
	double p_value_cut_off = 0.01;
	double cor_cut_off = 0.6;
	float top_percentage = 0.01;
	int max_degree = 10000;
	bool leave_one_out=false;

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
		case 'n':
			name = optarg;
			break;
		case 'p':
			p_value_cut_off = atof(optarg);
			break;
		case 'c':
			cor_cut_off = atof(optarg);
			break;
		case 't':
			top_percentage = atof(optarg);
			break;
		case 'd':
			max_degree = atoi(optarg);
			break;
		case 'l':
			leave_one_out = true;
			break;
		case '?':
			print_usage(stderr,-1);
		case -1:
			break;
		default:
			abort();
		}
	}while(next_option!=-1);
	if (name==NULL)
	{
		name = output_filename;
	}
	//ifstream inf(argv[1]);
	//ofstream outf(argv[2], ios::app | ios::out);

	if (optind < argc)
	{
		
		graph_construct instance(argv[optind], output_filename, name, p_value_cut_off, cor_cut_off, top_percentage, \
			max_degree, leave_one_out);
		instance.run();
		/*
		//testing
		vector<int> i_vector;
		i_vector.push_back(0);
		i_vector.push_back(1);
		graph_construct instance(output_filename, i_vector);
		vector<string> s_vector;
		vector<float> fv;
		s_vector.push_back("1.4320");
		s_vector.push_back("2.4310");
		s_vector.push_back("3.0");
		s_vector.push_back("4.0");
		s_vector.push_back("2.10");
		s_vector.push_back("5.3120");
		s_vector.push_back("1.23");
		s_vector.push_back("8.32");
		//instance.gene_label2index_setup(s_vector);
		instance.gene_array_fill(s_vector, 2);
		fv=instance.edge_construct_no_cut_off();
		vector<float>::iterator fv_iter;
		for(fv_iter=fv.begin(); fv_iter!=fv.end(); fv_iter++)
			cout<<*fv_iter<<endl;
		*/
		
	}
	else
		print_usage(stderr, 1);
}
