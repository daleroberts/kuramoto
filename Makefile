# detect compiler
ifeq ($(shell which icpc &>/dev/null; echo $$?),0)
CXX=icpc
CXXFLAGS=-Wall -std=c++11 -O2 -lboost_serialization-mt -lboost_mpi-mt -L/apps/boost/1.57.0/lib
else
CXX=g++
CXXFLAGS=-Wall -std=c++11 -O2 -lm -lboost_serialization-mt -lboost_mpi-mt -L/apps/boost/1.57.0/lib
endif

OBJS=$(patsubst %.cc,%.o,$(wildcard *.cc))
DEPS=$(OBJS:.o=.d)
EXEC=kuramoto_mpi

all: kuramoto_mpi

kuramoto_omp: kuramoto_omp.o graph.o statistics.o
	$(CXX) $(CXXFLAGS) $^ -o $@

kuramoto_onepath: kuramoto_onepath.o graph.o statistics.o
	$(CXX) $(CXXFLAGS) $^ -o $@

kuramoto_mpi: kuramoto_mpi.o graph.o statistics.o
	$(CXX) $(CXXFLAGS) `mpic++ -showme:link` $^ -o $@

-include $(DEPS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -MMD -c $< -o $@

clean:
	@find . -executable -type f -delete
	@rm -fr $(OBJS) $(DEPS) *.o*
