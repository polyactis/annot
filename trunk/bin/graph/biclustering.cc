/*
*04-12-05
*	biclustering algorithm modelled after Biclustering.java by Cheng & Church 2000
*
*/

#include <iostream>			// for std::cout
#include <vector>			//for vector
#include <string>			//for string
#include <gsl/gsl_rng.h>		//for random number generator
#include <boost/python.hpp>	//for python module and dict, tuple, make_tuple

using namespace std;
using namespace boost::python;

typedef vector<double> vd;
typedef vector<bool> vb;
typedef vector<float> vf;

//data structure to hold the information
struct bicluster
{
	boost::python::list row_index_list;
	boost::python::list column_index_list;
	boost::python::list consensus_list;
	double score;
};

class biclustering
{
	public:
		biclustering(double mxScore, int mnHeight, int mnWidth, int bThreshold);
		~biclustering();
		
		void data_read_in(boost::python::numeric::array& list_2d);
		void scoring();
		bicluster getbicluster();
		boost::python::list return_matrix_data();
	
		int numberOfColumns;
		int numberOfRows;
		vector<vf> matrix;
		double maxScore;
		int minHeight;
		int minWidth;
		int batchThreshold;
		vb remainingR;
		vb remainingC;
		vd rowMean;
		vd columnMean;
		vd rowScore;
		vd columnScore;
		double mean;
		int smWidth;
		int smHeight;
		double HScore;
	
		gsl_rng * r;	//random number generator
};

biclustering::biclustering(double mxScore, int mnHeight, int mnWidth, int bThreshold)
{
	//parameter initialize
	maxScore = mxScore;
	minHeight = mnHeight;
	minWidth = mnWidth;
	batchThreshold = bThreshold;
	
	//initialize the random number generator.
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
	r = gsl_rng_alloc (T);
}

biclustering::~biclustering()
{
	//do nothing
	std::cout<<"Biclustering c++ exits"<<endl;
}

void biclustering::data_read_in(boost::python::numeric::array& list_2d)
{
	//convert the 2 dimensional list of python into matrix
	boost::python::tuple shape = boost::python::extract<tuple>(list_2d.attr("shape"));
	numberOfRows = boost::python::extract<int>(shape[0]);
	if (numberOfRows == 0)
	{
		cout<<"No data"<<endl;
		numberOfColumns = 0;
	}
	else
		numberOfColumns = boost::python::extract<int>(shape[1]);
	cout<<"Number of rows: "<<numberOfRows<<endl;
	cout<<"Number of columns: "<<numberOfColumns<<endl;
	
	//fill matrix
	vf row_vector;	//temporarily hold the data before pushed into matrix
	float value;
	for(int i=0; i<numberOfRows; i++)
	{
		for(int j=0; j<numberOfColumns; j++)
		{
			value = boost::python::extract<float>(list_2d[i][j]);
			/*
			if ((value-1.1)<0.0001)
				value = float(gsl_rng_uniform_int(r, 1600)-800);
			else
			{
				value = atof(value_string.c_str());
				//see log on 04-14-05 in section of c++/c. tricks in the data type coersion.
				#if defined(DEBUG)
					std::cout << "Original: "<<value_string<<"Middle: "<<atof(value_string.c_str())<<" New: "<<value<<endl;
				#endif
			}
			*/
			row_vector.push_back(value);
		}
		matrix.push_back(row_vector);
		row_vector.clear();	//clear it for the next time
	}
	
	//initialize other sequence structures
	for(int i=0; i<numberOfRows; i++)
	{
		remainingR.push_back(true);
		rowMean.push_back(0);
		rowScore.push_back(0);
	}
	for(int j=0; j<numberOfColumns; j++)
	{
		remainingC.push_back(true);
		columnMean.push_back(0);
		columnScore.push_back(0);
	}
}

void biclustering::scoring()
{
	/*
	*copied from Biclustering.java, no change
	*04-26-05 add some DEBUG macros to control compilation
	*/
	mean = 0;
	for (int j = 0; j < numberOfColumns; j++)
		if (remainingC[j])
			columnMean[j] = 0;	//initialization
	for (int i = 0; i < numberOfRows; i++)
		if (remainingR[i])
		{
			rowMean[i] = 0;	//initialization
			for (int j = 0; j < numberOfColumns; j++)
				if (remainingC[j])
				{
					rowMean[i] += matrix[i][j];
					columnMean[j] += matrix[i][j];
				}
			mean += rowMean[i];
			rowMean[i] /= smWidth;	//smWidth determined in getBicluster()
			#if defined(DEBUG)
				cout<<"rowMean at "<<i<<" is "<<rowMean[i]<<"; incremental mean is "<<mean<<endl;
			#endif
		}
	for (int j = 0; j < numberOfColumns; j++)
		if (remainingC[j])
		{
			columnMean[j] /= smHeight;	//smHeight determined in getBicluster()
			#if defined(DEBUG)
				cout<<"colMean at "<<j<<" is "<<columnMean[j]<<endl;
			#endif
		}
	#if defined(DEBUG)
		cout<<"mean before dividing is "<<mean<<endl;
		cout<<"smWidth is "<<smWidth<<endl;
		cout<<"smHeight is "<<smHeight<<endl;
	#endif
	mean /= smWidth * smHeight;
	#if defined(DEBUG)
		cout<<"mean is "<<mean<<endl;
	#endif
	HScore = 0;
	for (int j = 0; j < numberOfColumns; j++)
		if (remainingC[j])
			columnScore[j] = 0;	//initialization
	for (int i = 0; i < numberOfRows; i++)
		if (remainingR[i])
		{
			rowScore[i] = 0;	//initialization
			for (int j = 0; j < numberOfColumns; j++)
				if (remainingC[j])
				{
					double r = matrix[i][j] - rowMean[i] - columnMean[j] + mean;
					#if defined(DEBUG)
						cout<<"r score at "<<i<<" and "<<j<<" is "<<r<<endl;
					#endif
					r = r * r;
					rowScore[i] += r;
					columnScore[j] += r;
				}
			HScore += rowScore[i];
			rowScore[i] /= smWidth;
			#if defined(DEBUG)
				cout<<"rowScore at "<<i<<" is "<<rowScore[i]<<endl;
			#endif
		}
	#if defined(DEBUG)
		cout<<"HScore before dividing is "<<HScore<<endl;
	#endif
	HScore /= smWidth * smHeight;
	for (int j = 0; j < numberOfColumns; j++)
		if (remainingC[j])
			columnScore[j] /= smHeight;
}

bicluster biclustering::getbicluster()
/*
*04-13-05
*	return structure bicluster, not tuple.
*/
{
	bicluster bicluster_to_return;
	for (int i = 0; i < numberOfRows; i++)
		remainingR[i] = true;
	for (int j = 0; j < numberOfColumns; j++)
		remainingC[j] = true;
	smWidth = numberOfColumns;
	smHeight = numberOfRows;
	scoring();
	int index = 0;
	while ((HScore > maxScore) && (index > -1))
	{	//if no more columns or rows can be removed,
		if (smHeight > batchThreshold)
		{			//This algorithm 2 of the paper. Multiple node deletion, but kind of different,
					//no batchThreshold in the paper's algorithm
			for (int i = 0; i < numberOfRows; i++)
				if (remainingR[i] && (rowScore[i] > HScore))
				{
					remainingR[i] = false;
					smHeight--;
				}
		}else
		{	//find the maximum from rowScore and columnScore
			double ms = 0;
			index = -1;
			bool row = true;
			if (smHeight > minHeight)
			{
				for (int i = 0; i < numberOfRows; i++)
					if (remainingR[i] && (rowScore[i] > ms))
					{
						ms = rowScore[i];
						index = i;
					}
			}
			if (smWidth > minWidth)
			{
				for (int i = 0; i < numberOfColumns; i++)
					if (remainingC[i] && (columnScore[i] > ms))
					{
						ms = columnScore[i];
						index = i;
						row = false;
					}
			}
			if (index > -1)
				if (row)
				{
					remainingR[index] = false;
					smHeight--;
				}else
				{
					remainingC[index] = false;
					smWidth--;
				}
		}
		scoring();
	}
	
	list row_index_list;
	list column_index_list;
	//put random numbers in the block(cluster) already discovered
	for (int i = 0; i < numberOfRows; i++)
		if (remainingR[i])
		{
			row_index_list.append(i);	//fill the row index of the cluster
			for (int j = 0; j < numberOfColumns; j++)
				if (remainingC[j])
					matrix[i][j] = int(gsl_rng_uniform_int(r, 1600)-800);
			}
	//fill the column index of the cluster
	for (int j=0; j<numberOfColumns; j++)
		if(remainingC[j])
		{
			column_index_list.append(j);
			bicluster_to_return.consensus_list.append(columnMean[j]);
		}

	bicluster_to_return.row_index_list = row_index_list;
	bicluster_to_return.column_index_list = column_index_list;
	bicluster_to_return.score = HScore;
	return bicluster_to_return;
}

boost::python::list biclustering::return_matrix_data()
{
	list list_2d;
	for(int i=0; i<numberOfRows; i++)
	{
		list list_row;
		for(int j=0; j<numberOfColumns; j++)
		{
			list_row.append(matrix[i][j]);
		}
		list_2d.append(list_row);
	}
	return list_2d;
}

BOOST_PYTHON_MODULE(biclustering)
{
	boost::python::class_<biclustering>("biclustering", init<double, int, int, int>())
		.def("data_read_in", &biclustering::data_read_in)
		.def("scoring", &biclustering::scoring)
		.def("getbicluster", &biclustering::getbicluster)
		.def("return_matrix_data", &biclustering::return_matrix_data)
		.def_readonly("maxScore", &biclustering::maxScore)
		.def_readonly("minHeight", &biclustering::minHeight)
		.def_readonly("minWidth", &biclustering::minWidth)
		.def_readonly("batchThreshold", &biclustering::batchThreshold)
	;
	boost::python::class_<bicluster>("bicluster")
		.def_readonly("row_index_list", &bicluster::row_index_list)
		.def_readonly("column_index_list", &bicluster::column_index_list)
		.def_readonly("consensus_list", &bicluster::consensus_list)
		.def_readonly("score", &bicluster::score)
	;
	def("set_module_and_type", &numeric::array::set_module_and_type);
}
