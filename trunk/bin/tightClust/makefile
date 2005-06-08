CFLAGS	=

LDFLAGS	= -lm

CC	= gcc

CXX	= g++

.c.o:
	$(CC) -c $(CFLAGS) $<

TARGET	= tightClust

all:	$(TARGET)

OBJS	= cluster.o com.o linpack.o ranlib.o tightClust.o

$(TARGET): $(OBJS)
	$(CC) -o $@ $(OBJS) $(LDFLAGS)

clean:
	-rm -f $(TARGET)
	-rm -f *.o *.bak
	-rm -f core
