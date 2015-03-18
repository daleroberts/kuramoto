COMPILER=gcc

ifeq ($(COMPILER),gcc)
CC=g++-4.9
CCFLAGS=-O3 -std=c++11 -fopenmp -I/usr/local/include/eigen3 -I.
endif

ifeq ($(COMPILER),icc)
CC=icc
CCFLAGS=-O3 -std=c++11 -openmp -I/usr/local/include/eigen3 -I.
endif

all: kuramoto kuramoto_onepath

%.o: %.cc variates.h graph.h statistics.h
	$(CC) $(CCFLAGS) -c $<

kuramoto: graph.o statistics.o kuramoto.o
	$(CC) $(CCFLAGS) $^ -o $@

kuramoto_onepath: graph.o statistics.o kuramoto_onepath.o
	$(CC) $(CCFLAGS) $^ -o $@

dist: dist.o
	$(CC) $(CCFLAGS) $^ -o $@

clean:
	-rm -f kuramoto dist kuramoto_onepath *.o
