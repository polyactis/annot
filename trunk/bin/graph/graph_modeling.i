%module graph_modeling
%{
#include "graph_modeling.h"
%}

%include "std_vector.i"
%include "std_string.i"
namespace std {
	%template(IntVector) vector<int>;
	%template(FloatVector) vector<float>;
	%template(DoubleVector) vector<double>;
	%template(StringVector) vector<string>;
}
%include "graph_modeling.h"
