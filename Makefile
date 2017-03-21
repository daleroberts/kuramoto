CXX=g++ -std=c++11

CXXFLAGS=-Wall -Wno-deprecated-declarations
MPIFLAGS=$(shell mpic++ -showme:compile)
EIGENFLAGS=$(shell pkg-config eigen3 --cflags)

MPILIBS=$(shell mpic++ -showme:link)
LIBS=-lboost_serialization -lboost_mpi

OBJS=$(patsubst %.cc,%.o,$(wildcard *.cc))
DEPS=$(patsubst %.cc,%.d,$(wildcard *.cc))

kuramoto: kuramoto.o graph.o statistics.o variates.o
	$(CXX) $(MPILIBS) $(LIBS) $^ -o $@

clean:
	@rm -fr *.d $(OBJS) $(EXEC)
	
-include $(DEPS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(MPIFLAGS) $(EIGENFLAGS) -MMD -c $< -o $@
