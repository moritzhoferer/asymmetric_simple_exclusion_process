CXX = g++
CXXFLAGS = -lm -lgsl -lgslcblas

all: constructor

clean:
	rm -f *.o Parameters* option.txt constructor

constructor: constructor.o
	$(CXX) -o constructor constructor.o $(CXXFLAGS)

constructor.o: constructor.cpp
	$(CXX) $(CXXFLAGS) -c constructor.cpp
