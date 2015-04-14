BC?=release
CC=g++
SHELL=/bin/bash
CFLAGS=-Wall -lm -std=c++11 -I. -I${HOME}/.local/include
 
ifeq ($(BC),debug)
	CFLAGS += -g3
else
	CFLAGS += -O2
endif

OBJS=$(patsubst %.cc,%.o,$(wildcard *.cc))
DEPS=$(OBJS:.o=.d)
EXEC=kuramoto_mpi kuramoto

all: $(EXEC)

kuramoto: kuramoto.o graph.o statistics.o
	$(CC) $(CFLAGS) -fopenmp $^ -o $@

kuramoto_mpi: kuramoto_mpi.o graph.o statistics.o
	$(CC) -L${HOME}/.local/lib `mpic++ -showme:link` -lboost_serialization -lboost_mpi $^ -o $@

-include $(DEPS)

%.o: %.cc
	$(CC) $(CFLAGS) -MMD -c $< -o $@

clean:
	rm -fr $(EXEC) $(OBJS) $(DEPS)
