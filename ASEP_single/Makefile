CXX = g++
CXXFLAGS = -lm -lgsl -lgslcblas

all: ASEP_single

clean:
	rm -f *.o ASEP_single

ASEP_single_core.o: ASEP_single_core.cpp ASEP_single_core.h
	$(CXX) $(CXXFLAGS) -c ASEP_single_core.cpp

ASEP_single_main.o: ASEP_single_main.cpp ASEP_single_core.h timer.h
	$(CXX) $(CXXFLAGS) -c ASEP_single_main.cpp

ASEP_single: ASEP_single_main.o ASEP_single_core.o
	$(CXX) -o ASEP_single ASEP_single_main.o ASEP_single_core.o $(CXXFLAGS) 