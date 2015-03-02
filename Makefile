COMPILER=gcc

ifeq ($(COMPILER),gcc)
CC=g++-4.9
CCFLAGS=-O3 -std=c++11 -fopenmp -I/usr/local/include/eigen3 -I.
endif

ifeq ($(COMPILER),icc)
CC=icc
CCFLAGS=-O3 -std=c++11 -openmp -I/usr/local/include/eigen3 -I.
endif

all: kuramoto

%.o: %.cc variates.h graph.h statistics.h
	$(CC) $(CCFLAGS) -c $<

kuramoto: graph.o statistics.o kuramoto.o
	$(CC) $(CCFLAGS) $^ -o $@

test_stat: test_stat.o statistics.o
	$(CC) $(CCFLAGS) $^ -o $@

clean:
	-rm -f kuramoto test_stat test *.o
