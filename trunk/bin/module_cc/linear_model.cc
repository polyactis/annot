/*
* 02-28-05 linear model fitting module for python, boost.python + gsl
*
*/
#include <boost/python.hpp>
#include <gsl/gsl_multifit.h>		//for gsl linear model stuff
using namespace boost::python;

class linear_model
{
	public:
	linear_model();
	~linear_model();
	void prepare_data(list y_list, list x_2d_list);
	void run();
	list coefficients();
	tuple chisq_return();
	
	double chisq;
	gsl_matrix *X, *cov;
	gsl_vector *y, *c;
	gsl_multifit_linear_workspace * work;
};

linear_model::linear_model()
{
}

linear_model::~linear_model()
{
	gsl_matrix_free(X);
	gsl_matrix_free(cov);
	gsl_vector_free(y);
	gsl_vector_free(c);
	gsl_multifit_linear_free (work);
}

void linear_model::prepare_data(list y_list, list x_2d_list)
{
	/*
	*02-28-05
	*	initialize all the matrices and vectors
	*	transfer the values from python list to matrices
	*/
	int y_dimension = extract<int>(y_list.attr("__len__")());
	int x_dimension = extract<int>(x_2d_list[0].attr("__len__")());
	X = gsl_matrix_alloc(y_dimension, x_dimension);
	cov = gsl_matrix_alloc(x_dimension, x_dimension);
	y = gsl_vector_alloc(y_dimension);
	c = gsl_vector_alloc(x_dimension);
	work = gsl_multifit_linear_alloc (y_dimension, x_dimension );
	for (int i=0; i<y_dimension; i++)
	{
		gsl_vector_set(y, i, extract<double>(y_list[i]));
		for(int j=0; j<x_dimension; j++)
		{
			gsl_matrix_set(X, i, j, extract<double>(x_2d_list[i][j]));
		}
	}
}

void linear_model::run()
{
	gsl_multifit_linear (X, y, c, cov, &chisq, work);
}

list linear_model::coefficients()
{
	/*
	*02-28-05
	*	return a python list of coefficients
	*
	*/
	list coeff_list;
	int size = c->size;
	for (int i=0; i<size; i++)
		coeff_list.append(gsl_vector_get(c,i));
	return coeff_list;
}

tuple linear_model::chisq_return()
{
	/*
	*02-28-05
	*	return the chisq in a tuple
	*/
	return boost::python::make_tuple(chisq);
}

BOOST_PYTHON_MODULE(linear_model)
{
	class_<linear_model>("linear_model")
		.def("prepare_data", &linear_model::prepare_data)
		.def("run", &linear_model::run)
		.def("coefficients", &linear_model::coefficients)
		.def("chisq_return", &linear_model::chisq_return)
	;
}
