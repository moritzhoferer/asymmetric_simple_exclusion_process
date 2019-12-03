//============================================================================
// Name        : ASEP_Moritz.cpp
// Version     : 1.0
// Author      : Moritz Hoferer
// Date		   : 2015/08/07
// First run   : 2015/08/23
// Checked	   : In progress
// Copyright   : no
// Description : Simulation of ASEP with different possibilities of actions, system
//				 configurations, measurements:
//					- hopping right and left
//					- being pushed and repelled
//					- attachment and detachment of particles
//					- open and periodic boundary conditions
//============================================================================


#include <iostream>
#include <string>
#include <iomanip>

#include "ASEP_single_core.h"
#include "timer.h"

using namespace std;

int main(int argc, char *argv[]) {
	double Parameters[21];

	// initialization of the timer
	timer control;
	// start the timer
	control.start();

	// Initialization of the output files
	ofstream file0, file1,file2;
	string str, identifier, out0, out1, out2;

	// define name of input file
	str = "Parameters";
	str.append(argv[1]);
	str.append(".txt");

	// read in the control parameter for lattice configuration and particle properties
	read_in_parameters(str, Parameters);
	for(unsigned int i = 0; i < 21; i++){cout << Parameters[i] << endl;	}

	// initialization of the test system
	Lattice sys(Parameters);

	// construct an identifier
	// has to be adjusted for different measurement
	if(bool(Parameters[7])){
		identifier = "rho_0_";
		identifier.append(do_to_str(Parameters[19]));
	}
	else{
		identifier = "in_left_";
		identifier.append(do_to_str(Parameters[15]));
		identifier.append("_out_right_");
		identifier.append(do_to_str(Parameters[18]));
	}
	// construction of the names of the output files (use initial density as identifier)
	out0 = "parameter_"+ identifier+".txt";
	out1 = "density_"+ identifier+".txt";
	out2 = "current_"+ identifier+".txt";

	// prints the initial state in the console
	sys.print_state();

	// evolution of the simulation step by step. Check every time, if t_max > t_now, otherwise stop simulation.
	while (sys.not_finished()){
		sys.simulate();
		// Remove this options for better performance
		//sys.print_state();
		//sys.print_density();
		//sys.print_current();
	}

	// determines the final spatial density rho
	sys.get_final_rho();

	// open the output file
	file0.open(out0.c_str());
	if(sys.recording_density()){file1.open(out1.c_str());}
	if(sys.recording_current()) {file2.open(out2.c_str());}

	// write the lattice parameters and particle properties into the output files
	sys.write_parameter(file0);
	if(sys.recording_density()) {
		sys.write_parameter(file1);
		sys.write_density(file1);}
	if(sys.recording_current()) {
		sys.write_parameter(file2);
		sys.write_current(file2);}

	// close the output files
	file0.close();
	if(sys.recording_density()) {file1.close();}
	if(sys.recording_current()) {file2.close();}

	// output in the console after finishing the simulation
	cout << "Number of iteration steps:" << sys.get_count() << endl;
	sys.print_state();

	// stop the time measurement
	control.stop();
	control.check();

	return 0;
}
