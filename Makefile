COMPILER = gcc

ifeq ($(COMPILER),gcc)
CC=g++
CCFLAGS=-O3 -std=c++11 -I/usr/local/include/eigen3 -I.
endif

ifeq ($(COMPILER),icc)
CC=icc
CCFLAGS=-O3 -std=c++11 -openmp -I/usr/local/include/eigen3 -I.
endif

all: kuramoto

%.o: %.cc variates.h graph.h
	$(CC) $(CCFLAGS) -c $<

kuramoto: graph.o kuramoto.o
	$(CC) $(CCFLAGS) $^ -o $@

clean:
	-rm -f kuramoto *.o
