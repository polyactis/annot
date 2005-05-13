#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>
#include <string>
#include <hash_map.h>
#include <getopt.h>
#include <assert.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_statistics_double.h>
using namespace std;

typedef vector<float> vf;

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
void cor_cut_off_vector_construct(double p_value_cut_off, double cor_cut_off_given);

edge euc_dist(vf v1, vf v2);

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
		graph_construct(char* inf_name, char* outf_name, char* g_name);
		graph_construct(char* outf_name, vector<int> edge_vector);
		~graph_construct();
		vector<float> cor_cut_off_array_construct(double p_value_cut_off, double cor_cut_off_given, int max_degree);
		void input();
		void edge_construct(bool leave_one_out);	//04-30-05	flag leave_one_out to control whether leave_one_out or not.
		void output();
		edge min_cor(vf v1, vf v2);
		edge cor(vf v1, vf v2, int position);
		edge d_distance(vf v1, vf v2);
		vector<string> general_split(string line, char ch);
		void split(string line);
		
		//
		void gene_label2index_setup(vector<string> label_vector);
		void gene_array_fill(vector<string> expr_array, int no_of_genes);
		vf edge_construct_no_cut_off();
	
		string graph_name;
		hash_map<const char*, int > gene_label2index;
		vector<int> edge_tuple_vector;


};
