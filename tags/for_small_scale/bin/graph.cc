#include <iostream>
#include <fstream>
#include <vector>
#include <bitset>
#include <cmath>
#define GENE_CUT_OFF 8
#define JK_CUT_OFF 7
#define COR_CUT_OFF 0.6
using namespace std;

typedef vector<float> vf;

//Usage:	'program_name' arg1 arg2 arg3
//			arg1 is the dataset source file
//			arg2 is the output file to hold the result
//			arg3 is the name of the dataset(whatever you like).
//			example: ./program_name dataset1 dataset1.gph dataset1

class graph_construct
{
	vector<vf> gene_array;
	vector<bit_vector> mask_vector;
	vector<string> gene_labels_vector;
	char* graph_name;
	int no_of_genes;
	int no_of_cols;
	ifstream in;
	ofstream out;
	public:
		graph_construct(char* inf_name, char* outf_name, char* g_name);
		~graph_construct();
		void input();
		void edge_construct();
		void output();
		float cor(vf v1, vf v2, int position);
		bit_vector bv_or(bit_vector bv1, bit_vector bv2);
		int bv_count(bit_vector bv);
		void split(string line);
};

graph_construct::graph_construct(char* inf_name, char* outf_name, char* g_name)
{
	in.open(inf_name);
	out.open(outf_name);
	graph_name = g_name;
	//ios::app | ios::out);

}

graph_construct::~graph_construct()
{
	in.close();
	out.close();
}

void graph_construct::input()
{
	string line;
	while(getline(in, line))
	{
		split(line);
	}
	no_of_genes = gene_array.size();
	no_of_cols = gene_array[0].size();
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
					gene_vector.push_back(NAN);
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
		gene_vector.push_back(NAN);
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

void graph_construct::edge_construct()
{
	int i,j,k;
	float tmp_cor, min_cor = 1.1;
	
	out<<"t\t#\t"<<graph_name<<endl;
	for ( i=0; i<no_of_genes; i++)
		for ( j=i+1; j<no_of_genes; j++)
		{	
			min_cor = 1.1;
			for (k=0; k<no_of_cols; k++)
			{
				tmp_cor = cor(gene_array[i], gene_array[j], k);
				//cerr<<tmp_cor<<'\t';
				//if #valid quantities shared by [i] and [j] is <7, it will return 2.0. 
				if ( abs(tmp_cor) < abs(min_cor) )
					min_cor = tmp_cor;
			}
			//cerr<<endl;
			if ( (abs(min_cor) >= COR_CUT_OFF) && (abs(min_cor)<=1.0) )
				out<<"e\t"<<gene_labels_vector[i]<<'\t'<<gene_labels_vector[j]<<'\t'<<min_cor<<endl;
		}			
}

void graph_construct::output()
{
	int i;
	out<<no_of_cols<<'\t'<<no_of_genes<<endl;
		
	for ( i=0; i<no_of_genes; i++)
	{
		out<<gene_labels_vector[i]<<'\t';
		for (int j=0; j<no_of_cols; j++)
			out<<mask_vector[i][j];
		out<<endl;
		for (int j=0; j<no_of_cols; j++)
			out<<gene_array[i][j]<<'\t';
	
	out<<endl;
	}

}


float graph_construct::cor(vf v1, vf v2, int position)
{
	float xx = 0.0;
	float yy = 0.0;
	float xy = 0.0;
	float mean_x = 0.0;
	float mean_y = 0.0;
	int no_of_valids=0;
	
	for (int i=0; i<no_of_cols; i++)
		if ( (i == position) || ( isnan(v1[i]) )|| ( isnan(v2[i]) ) )
			continue;
		else
		{
			mean_x += v1[i];
			mean_y += v2[i];
			no_of_valids++;
		}
	
	if (no_of_valids<JK_CUT_OFF)
	{	
		return 2.0;
	}
	
	mean_x /= no_of_valids;
	mean_y /= no_of_valids;
	
	for (int i=0; i<no_of_cols; i++)
		if ( (i == position) || ( isnan(v1[i]) )|| ( isnan(v2[i]) ) )
			continue;
		else
		{
			xy += (v1[i]-mean_x) * (v2[i]-mean_y);
			xx += (v1[i]-mean_x) * (v1[i]-mean_x) ;
			yy += (v2[i]-mean_y)  * (v2[i]-mean_y) ;
		}
	return xy/(sqrt(xx)*sqrt(yy));
}

bit_vector graph_construct::bv_or(bit_vector bv1, bit_vector bv2)
{
}

int graph_construct::bv_count(bit_vector bv)
{
	int count=0;
	for (int i=0; i<bv.size(); i++)
		bv[i]?count++:count;
	return count;
}

int main(int argc, char* argv[])
{
	//ifstream inf(argv[1]);
	//ofstream outf(argv[2], ios::app | ios::out);
	graph_construct instance(argv[1], argv[2], argv[3] );
	instance.input();
	instance.edge_construct();
	//instance.output();
}
