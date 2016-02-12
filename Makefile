# detect compiler
ifeq ($(shell which icpc &>/dev/null; echo $$?),0)
CXX?=icpc
CXXFLAGS=-Wall -std=c++11 -lboost_serialization-mt -lboost_mpi-mt -L/apps/boost/1.59.0/lib
else
CXX?=g++
CXXFLAGS=-Wall -std=c++11 -I/apps/eigen/3.2.1/include/eigen3
LIBS=-L/apps/openmpi/1.10.0/lib -L/apps/boost/1.59.0/lib -lm -lmpi -lboost_serialization-mt -lboost_mpi-mt
endif

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
	@rm *.o
	@rm -fr $(OBJS) $(DEPS)
	@rm $(EXEC)
