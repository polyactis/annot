BoostPythonTarget	= clustering.so biclustering.so
BoostPythonLib  = -lboost_python

swig_py		= graph_modeling.py
cxx_wrap	= graph_modeling_wrap.cc
cxx_swig_targets	= _graph_modeling.so

swig_include	= -I/usr/include/python2.3
swig_flags	= -shared -fPIC

.SUFFIXES:	.i .cc .c .so .o _wrap.o .py
.PHONY:	swig

gsl_lib		= -lgsl -lgslcblas
gsl_cxx_targets	= graph_modeling
cxx_flags	= -O2 -fPIC
cxx_compiler	= g++
gsl_c_targets	= 
c_flags		= -O2
c_compiler	= gcc


all:	$(gsl_cxx_targets) $(gsl_c_targets) $(BoostPythonTarget) $(cxx_swig_targets)

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
	$(cxx_compiler)  $(gsl_lib) $(swig_flags) $^ -o $@


#for boost.python
$(BoostPythonTarget):	%.so:	%.o
	$(cxx_compiler) $(gsl_lib) $(BoostPythonLib) $(swig_flags) $< -o $@
clean:
	-rm -f *.class *.o $(gsl_cxx_targets) $(gsl_c_targets) $(cxx_swig_targets) $(cxx_wrap) *.pyc $(swig_py) $(BoostPythonTarget)
