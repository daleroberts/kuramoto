BC?=release
CC=g++-4.9
SHELL=/bin/bash
CFLAGS=-Wall -fopenmp -lm -std=c++11 -I. -I/usr/local/include/eigen3
 
ifeq ($(BC),debug)
	CFLAGS += -g3
else
	CFLAGS += -O2
endif

OBJS=$(patsubst %.cc,%.o,$(wildcard *.cc))
DEPS=$(OBJS:.o=.d)

EXEC=kuramoto_mpi

all: kuramoto_mpi kuramoto

kuramoto: kuramoto.o graph.o statistics.o
	$(CC) $(CFLAGS) $^ -o $@

kuramoto_mpi: kuramoto_mpi.o graph.o statistics.o
	$(CC) -std=c++11 `mpic++ -showme:link` -lboost_serialization-mt -lboost_mpi-mt $^ -o $@

-include $(DEPS)

%.o: %.cc
	$(CC) $(CFLAGS) -MMD -c $< -o $@

clean:
	rm -fr $(EXEC) $(OBJS) $(DEPS)
