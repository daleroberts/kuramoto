BC?=release
CC=g++
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

$(EXEC):
	$(CC) $(CFLAGS) $^ -o $@

-include $(DEPS)

%.o: %.cc
	$(CC) $(CFLAGS) -MMD -c $< -o $@

clean:
	rm -fr $(EXEC) $(OBJS) $(DEPS)
