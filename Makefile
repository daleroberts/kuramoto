# detect compiler
ifeq ($(shell which icpc &>/dev/null; echo $$?),0)
CXX?=icpc
CXXFLAGS=-Wall -std=c++11 -lboost_serialization-mt -lboost_mpi-mt -L/apps/boost/1.57.0/lib
else
CXX?=g++
CXXFLAGS=-Wall -std=c++11 -lm -lboost_serialization-mt -lboost_mpi-mt -L/usr/local/lib -I/usr/local/include/eigen3 -L/apps/boost/1.57.0/lib
endif

OBJS=$(patsubst %.cc,%.o,$(wildcard *.cc))
DEPS=$(OBJS:.o=.d)
EXEC=kuramoto_mpi

all: kuramoto_mpi

kuramoto_omp: kuramoto_omp.o graph.o statistics.o variates.o
	$(CXX) $(CXXFLAGS) $^ -o $@

kuramoto_onepath: kuramoto_onepath.o graph.o statistics.o variates.o
	$(CXX) $(CXXFLAGS) $^ -o $@

kuramoto_x: kuramoto_x.o graph.o statistics.o variates.o
	$(CXX) $(CXXFLAGS) $^ -o $@

kuramoto_drift: kuramoto_drift.o graph.o statistics.o variates.o
	$(CXX) $(CXXFLAGS) $^ -o $@

kuramoto_mpi: kuramoto_mpi.o graph.o statistics.o variates.o
	$(CXX) $(CXXFLAGS) `mpic++ -showme:link` $^ -o $@

dist: dist.o variates.o
	$(CXX) $(CXXFLAGS) $^ -o $@

test_dist: test_dist.o variates.o
	$(CXX) $(CXXFLAGS) $^ -o $@

test_sum: test_sum.o graph.o
	$(CXX) $(CXXFLAGS) $^ -o $@

-include $(DEPS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) -MMD -c $< -o $@

clean:
	@find . \( -executable -type f ! -iname *.sh \) -delete
	@rm -fr $(OBJS) $(DEPS) *.o*
