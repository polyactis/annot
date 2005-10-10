ARPACKPP_INC	= ~/script/hhu_clustering/include/arpack++/include/

ARPACK_LIB	= ~/script/hhu_clustering/lib/libarpack_i686.a
SUPERLU_LIB	= ~/script/hhu_clustering/lib/superlu_linux.a
LAPACK_LIB	= ~/script/hhu_clustering/lib/liblapack.a
BLAS_LIB	= ~/script/hhu_clustering/lib/libblas.a
FORTRAN_LIB	= ~/script/hhu_clustering/lib/libg2c.a




BoostPythonTarget	= clustering.so biclustering.so cc_from_edge_list.so johnson_sp.so
BoostPythonBinary	= clustering cc_from_edge_list
BoostPythonLib  = -lboost_python

SwigInclude	= /usr/include/python2.3
CxxSwigTargets	= _graph_modeling.so
SwigPy		= graph_modeling.py
CxxWrap	= graph_modeling_wrap.cc


CFlags		= -O3
CxxFlags	= -O3 -fPIC
SharedLibFlags	= -shared -fPIC

.SUFFIXES:	.i .cc .c .so .o _wrap.o .py
.PHONY:	swig

GslLib		= -lgsl -lgslcblas
GslCTargets	= 
GslCxxTargets	= graph_modeling

CCompiler	= gcc
CxxCompiler	= g++

IncludeDirs	= -I $(ARPACKPP_INC) -I $(SwigInclude)

Libs	= $(GslLib) $(LEDA_LIB) $(ARPACK_LIB) $(SUPERLU_LIB) $(LAPACK_LIB) $(BLAS_LIB) $(FORTRAN_LIB)

all:	$(GslCxxTargets) $(GslCTargets) $(BoostPythonTarget) $(CxxSwigTargets) $(BoostPythonBinary)

swig:	$(CxxSwigTargets)

$(GslCxxTargets):	%:	%.o
	-$(CxxCompiler) $(CxxFlags) $(GslLib) $< -o $@

.cc.o:
	-$(CxxCompiler) $(CxxFlags) $(IncludeDirs) -c $< -o $@

$(GslCTargets):	%:	%.o
	-$(CCompiler) $(CFlags) $(GslLib) $< -o $@

.c.o:
	-$(CCompiler) $(CFlags) $(IncludeDirs) -c $< -o $@


#for swig
$(CxxWrap):	%_wrap.cc:	%.i %.cc
	swig -c++ -python -o $@ $<

%_wrap.o:	%_wrap.cc
	-$(CxxCompiler) $(IncludeDirs) $(CxxFlags) -c $< -o $@

$(CxxSwigTargets):	_%.so:	%.o %_wrap.o
	-$(CxxCompiler)  $(GslLib) $(SharedLibFlags) $^ -o $@


#for boost.python target(shared libraries)
$(BoostPythonTarget):	%.so:	%.o
	-$(CxxCompiler) $(GslLib) $(BoostPythonLib) $(SharedLibFlags) $< -o $@

#for boost python binary 
$(BoostPythonBinary):	%:	%.o
	-$(CxxCompiler) $(CxxFlags) $(BoostPythonLib) $< $(Libs) -o $@



clean:
	-rm -f *.class *.o $(GslCxxTargets) $(GslCTargets) $(CxxSwigTargets) $(CxxWrap) *.pyc $(SwigPy) $(BoostPythonTarget) $(BoostPythonBinary)
