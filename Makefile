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

EXEC=kuramoto

all: $(EXEC)

kuramoto: kuramoto.o graph.o statistics.o

test_reduce: test_reduce.cc
	$(CC) -std=c++11 `mpic++ -showme:compile` -lmpi_cxx -lmpi -lboost_serialization-mt -lboost_mpi-mt $^ -o $@

test_mpi: test_mpi.cc statistics.o
	$(CC) -std=c++11 `mpic++ -showme:compile` -lmpi_cxx -lmpi -lboost_serialization-mt -lboost_mpi-mt $^ -o $@

$(EXEC):
	$(CC) $(CFLAGS) $^ -o $@

-include $(DEPS)

%.o: %.cc
	$(CC) $(CFLAGS) -MMD -c $< -o $@

clean:
	rm -fr $(EXEC) $(OBJS) $(DEPS)
