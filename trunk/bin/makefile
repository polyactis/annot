SUBDIRS	= module_cc graph
.PHONY:	clean_subdirs subdirs $(SUBDIRS)

all:	graph.bin subdirs

CFLAGS = -O2
.SUFFIXES:	.bin	.cc	.o

.cc.o:
	-g++ $(CFLAGS) -c  $< -o $@

.o.bin:
	-g++ $(CFLAGS) $< -o $@

subdirs:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir; \
	done
	
clean_subdirs:
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean; \
	done

clean:	clean_subdirs
	-rm *.bin *.o
