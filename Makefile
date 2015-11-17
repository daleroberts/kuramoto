# detect compiler
ifeq ($(shell which icpc &>/dev/null; echo $$?),0)
CXX?=icpc
CXXFLAGS=-Wall -std=c++11 -lboost_serialization-mt -lboost_mpi-mt -L/apps/boost/1.57.0/lib
else
CXX?=g++-5
CXXFLAGS=-Wall -std=c++11 -I/usr/local/include/eigen3
LDFLAGS=-lm -lboost_serialization-mt -lboost_mpi-mt -L/usr/local/lib
endif

OBJS=$(patsubst %.cc,%.o,$(wildcard *.cc))
DEPS=$(OBJS:.o=.d)
EXEC=kuramoto_mpi

all: kuramoto_mpi kuramoto_onepath

kuramoto_omp: kuramoto_omp.o graph.o statistics.o variates.o
	$(CXX) $(CXXFLAGS) $^ -o $@

kuramoto_onepath: kuramoto_onepath.o graph.o statistics.o variates.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@

kuramoto_mpi: kuramoto_mpi.o graph.o statistics.o variates.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) `mpic++ -showme:link` $^ -o $@

-include $(DEPS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -MMD -c $< -o $@

clean:
	@rm kuramoto_mpi kuramoto_onepath
	@rm -fr $(OBJS) $(DEPS) *.o*
