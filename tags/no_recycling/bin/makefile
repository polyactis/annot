all:	graph.bin

CFLAGS = -O2
.SUFFIXES:	.bin	.cc	.o

.cc.o:
	-g++ $(CFLAGS) -c  $< -o $@

.o.bin:
	-g++ $(CFLAGS) $< -o $@

clean:
	-rm *.bin *.o
