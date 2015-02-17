COMPILER = icc

ifeq ($(COMPILER),gcc)
CC=g++-4.9
CCFLAGS=-O3 -std=c++11 -fopenmp -I/usr/local/include/eigen3
endif

ifeq ($(COMPILER),icc)
CC=icc
CCFLAGS=-O3 -std=c++11 -openmp -I/usr/local/include/eigen3
endif


%.o: %.cc variates.h graph.h
	$(CC) $(CCFLAGS) -c $<

kuramoto: graph.o kuramoto.o
	$(CC) $(CCFLAGS) $^ -o $@

clean:
	-rm -f kuramoto *.o
