CXX?=g++
CXXFLAGS=-Wall -std=c++11 -Wno-deprecated-declarations
LIBS=-lmpi -lboost_serialization -lboost_mpi

OBJS=$(patsubst %.cc,%.o,$(wildcard *.cc))
DEPS=$(OBJS:.o=.d)
EXEC=kuramoto kuramoto_onepath

all: $(EXEC)

kuramoto_onepath: kuramoto_onepath.o graph.o statistics.o variates.o
	$(CXX) $(CXXFLAGS) $(LIBS) $^ -o $@

kuramoto: kuramoto.o graph.o statistics.o variates.o
	$(CXX) $(CXXFLAGS) $(LIBS) `mpic++ -showme:link` $^ -o $@

-include $(DEPS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -MMD -c $< -o $@

clean:
	@rm -fr $(OBJS) $(DEPS) $(EXEC)
