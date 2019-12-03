CXX = g++
CXXFLAGS = -lm -lgsl -lgslcblas

all: ASEP_ensemble

clean:
	rm -f *.o ASEP_ensemble

ASEP_ensemble: ASEP_Ensemble_main.o ASEP_Ensemble_core.o
	$(CXX) -o ASEP_ensemble ASEP_Ensemble_main.o ASEP_Ensemble_core.o $(CXXFLAGS)

ASEP_Ensemble_main.o: ASEP_Ensemble_main.cpp ASEP_Ensemble_core.h timer.h
	$(CXX) $(CXXFLAGS) -c ASEP_Ensemble_main.cpp

ASEP_Ensemble_core.o: ASEP_Ensemble_core.cpp ASEP_Ensemble_core.h timer.h
	$(CXX) $(CXXFLAGS) -c ASEP_Ensemble_core.cpp
