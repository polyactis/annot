BoostPythonTarget	= clustering.so biclustering.so cc_from_edge_list.so
BoostPythonBinary	= clustering
BoostPythonLib  = -lboost_python

swig_include	= -I/usr/include/python2.3
cxx_swig_targets	= _graph_modeling.so
swig_py		= graph_modeling.py
cxx_wrap	= graph_modeling_wrap.cc


c_flags		= -O2
cxx_flags	= -O2 -fPIC
shared_lib_flags	= -shared -fPIC

.SUFFIXES:	.i .cc .c .so .o _wrap.o .py
.PHONY:	swig

gsl_lib		= -lgsl -lgslcblas
gsl_c_targets	= 
gsl_cxx_targets	= graph_modeling

c_compiler	= gcc
cxx_compiler	= g++


all:	$(gsl_cxx_targets) $(gsl_c_targets) $(BoostPythonTarget) $(cxx_swig_targets) $(BoostPythonBinary)

swig:	$(cxx_swig_targets)

$(gsl_cxx_targets):	%:	%.o
	$(cxx_compiler) $(cxx_flags) $(gsl_lib) $< -o $@

.cc.o:
	$(cxx_compiler) $(cxx_flags) $(swig_include) -c $< -o $@

$(gsl_c_targets):	%:	%.o
	$(c_compiler) $(c_flags) $(gsl_lib) $< -o $@

.c.o:
	$(c_compiler) $(c_flags) -c $< -o $@


#for swig
$(cxx_wrap):	%_wrap.cc:	%.i %.cc
	swig -c++ -python -o $@ $<

%_wrap.o:	%_wrap.cc
	$(cxx_compiler) $(swig_include) $(cxx_flags) -c $< -o $@

$(cxx_swig_targets):	_%.so:	%.o %_wrap.o
	$(cxx_compiler)  $(gsl_lib) $(shared_lib_flags) $^ -o $@


#for boost.python target(shared libraries)
$(BoostPythonTarget):	%.so:	%.o
	$(cxx_compiler) $(gsl_lib) $(BoostPythonLib) $(shared_lib_flags) $< -o $@

#for boost python binary 
$(BoostPythonBinary):	%:	%.o
	$(cxx_compiler) $(cxx_flags) $(gsl_lib) $(BoostPythonLib) $< -o $@



clean:
	-rm -f *.class *.o $(gsl_cxx_targets) $(gsl_c_targets) $(cxx_swig_targets) $(cxx_wrap) *.pyc $(swig_py) $(BoostPythonTarget) $(BoostPythonBinary)
