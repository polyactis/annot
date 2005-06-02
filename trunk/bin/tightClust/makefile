CFlags	=

IncludeDirs	=

Libs	= -lm

CCompiler	= gcc

.c.o:
	$(CCompiler) -c $(IncludeDirs) $(CFlags) $<

TARGET	= tightClust

all:	$(TARGET)

OBJS	= cluster.o com.o linpack.o ranlib.o tightClust.o

$(TARGET): $(OBJS)
	$(CCompiler) -o $@ $(OBJS) $(Libs)

clean:
	-rm -f $(TARGET)
	-rm -f *.o *.bak
	-rm -f core
