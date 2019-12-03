CXX= g++
CXXFLAGS =  -lm -lgsl -lgslcblas

all: constructor_Ensemble

clean:
	rm -f *.o Parameters* option.txt constructor_Ensemble

constructor_Ensemble: constructor_Ensemble.o
	$(CXX) -o constructor_Ensemble constructor_Ensemble.o $(CXXFLAGS)

constructor_Ensemble.o: constructor_Ensemble.cpp
	$(CXX) $(CXXFLAGS) -c constructor_Ensemble.cpp
