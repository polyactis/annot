#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>
#include <string>
#include <map>		//for std::map
#include <queue>	//for priority_queue
#include <getopt.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics_double.h>
#include <algorithm>	//08-24-05 for std::copy
#include <iterator>	//08-24-05 for std::ostream_iterator
#include <boost/tokenizer.hpp>	//08-24-05 for tokenizer, parse input file
#include "boost/tuple/tuple.hpp"	//for boost::tuple	06-22-05

using namespace std;

typedef vector<float> vf;
typedef boost::tuple<std::string, std::string, float> edge_string_cor;

class cmp_edge	//06-22-05
{
public:
    bool operator() (const edge_string_cor & s1, const edge_string_cor & s2) const
    {
        return s1.get<2>() > s2.get<2>();	//compare the 2th element float value correlation
    }
};

struct edge
{
	//could be correlatoin, variance, distance
	float value;
	//degree of freedom, discard those invalid columns
	int degree;
	//03-02-05 significane flag
	int significance;
};

//once I thought wrapping a C++ class with a c function would solve the problem linking with python. But it's not.
void python_call(char* outf_name, vector<int> edge_vector, vector<string> expr_array, int no_genes);

/*
 * python interface
 * 	ind_min_cor
 * 		--ind_cor
 * */
edge ind_min_cor(vf v1, vf v2);
edge ind_cor(vf v1, vf v2, int position);
void set_jk_cut_off(int jk_cut_off);	//07-03-05
void cor_cut_off_vector_construct(double p_value_cut_off, double cor_cut_off_given);

edge euc_dist(vf v1, vf v2);

//05-27-05 declare to let python have the cor_cut_off_vector.
vector<float> cor_cut_off_vector_return(double p_value_cut_off, double cor_cut_off_given);

class graph_construct
{
	/*one sequence to call
	 *	graph_construct
	 *	cor_cut_off_array_construct
	 *	input
	 *		--split
	 *	edge_construct
	 *		--min_cor
	 *			--cor
	 *	output
	 *
	 * another sequence to call
	 * 	gene_array_fill
	 * 		--general_split
	 * 	edge_construct_no_cut_off
	 * 		--min_cor
	 * 			--cor
	 * */
	ifstream in;
	ofstream out;
	vector<vf> gene_array;
	vector<float> cor_cut_off_array;
	
	vector<bit_vector> mask_vector;
	vector<string> gene_labels_vector;
	int no_of_genes;
	int no_of_cols;
	int no_of_01;
	//gsl_histogram * histogram;
	public:
		graph_construct(char* inf_name, char* outf_name, char* g_name, double p_value_cut_off_given, \
			double cor_cut_off_given, float top_percentage_given, int max_degree_given, bool leave_one_out_given);
		graph_construct(char* outf_name, vector<int> edge_vector);
		~graph_construct();
		vector<float> cor_cut_off_array_construct(double p_value_cut_off, double cor_cut_off_given, int max_degree);
		int input(float top_percentage);
		void edge_construct(bool leave_one_out, int top_number);	//04-30-05	flag leave_one_out to control whether leave_one_out or not.
			//06-22-05	top_number used to get the top edges
		void output();
		edge min_cor(vf v1, vf v2);
		edge cor(vf v1, vf v2, int position);
		edge d_distance(vf v1, vf v2);
		vector<string> general_split(string line, char ch);
		void split(string line);
		
		void gene_label2index_setup(vector<string> label_vector);
		void gene_array_fill(vector<string> expr_array, int no_of_genes);
		vf edge_construct_no_cut_off();
	
		void run();	//06-22-05
	
		string graph_name;
		double p_value_cut_off;	//06-22-05	parameters saved in class variablee
		double cor_cut_off;
		int max_degree;
		bool leave_one_out;
		float top_percentage;
		
		std::map<std::string, int > gene_label2index;	//06-22-05	start to use map<>, not hash_map<>
		std::priority_queue<edge_string_cor, vector<edge_string_cor>, cmp_edge > edge_pq;
		vector<int> edge_tuple_vector;

};
