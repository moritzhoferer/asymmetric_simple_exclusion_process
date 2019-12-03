//============================================================================
// Name        : ASEP_Moritz.cpp
// Version     : 2.0
// Author      : Moritz Hoferer
// Date		   : 2015/09/08
// First run   : 2015/09/--
// Checked	   : In progress
// Copyright   : no
// Description : Simulation of ASEP with different possibilities of actions, system
//				 configurations, measurements:
//					- hopping right and left
//					- being pushed and repelled
//					- attachment and detachment of particles
//					- open and periodic boundary conditions
//				 Extension:
//					- Simulation of several systems simultaneously
//					- Measurement of a density profile over the whole ensemble at one/several time point(s)
//					- Initial density step with rho_0 (left) and rho_1 (right)
//============================================================================

#include <iostream>
#include <string>
#include <iomanip>

#include "ASEP_Ensemble_core.h"
#include "timer.h"

using namespace std;

int main(int argc, char *argv[]) {
	double Parameters[22];

	// initialization of the timer
	timer control;
	// start the timer
	control.start();

	// Initialization of the output files
	ofstream file1,file2;
	string str, identifier, out1, out2;

	// define name of input file
	str = "Parameters";
	str.append(argv[1]);
	str.append(".txt");

	// read in the control parameter for lattice configuration and particle properties
	read_in_parameters(str, Parameters);
	for(unsigned int i = 0; i < 22; i++){cout << Parameters[i] << endl;	}

	// initialization of the test system
	Lattice sys(Parameters);

	// initialization of the measurement
	Measurement erecord(Parameters);

	// Introduction of the ensemble
	// In this for-loop the system is simulated and the data is recorded
	for(unsigned int j = 0; j < erecord.get_ensemble(); j++){

		// reset the system and time t_now
		sys.reset();

		// prints the initial state in the console
		//sys.print_state();

		// set the count of the time points to zero
		int i = 0;

		// evolution of the simulation step by step. Check every time, if t_max > t_now, otherwise stop simulation.
		while (erecord.not_finished(sys)){

			// simulated the next time step depending on all possible actions and the rates for the actions
			sys.simulate_time();

			// Remove this options for better performance
			//sys.print_state();

			if(sys.get_time() > erecord.get_time_point(i) ){

				// insert here function to extract the density profile form the sys object into the erecord object
				erecord.record(sys,i);

				// raise the time point counter by one to record at the next time point the density profile
				i++;
			}

			// select a random action and execute it
			sys.simulate_movement();

		}
		// determines the final spatial density rho
		sys.get_final_rho();
	}

	// in this for-loop the data is written in separated files
	for(int i = 0; i <= erecord.get_numb_time_points(); i++){

		// construct an identifier
		// has to be adjusted for different measurement
		identifier = "rho_left_";
		identifier.append(do_to_str(Parameters[5]));
		identifier.append("_rho_right_");
		identifier.append(do_to_str(Parameters[6]));
		// extension of the identifier to get separate output files for the different time points
		identifier.append("_time_");
		identifier.append(do_to_str(erecord.get_time_point(i)));

		// construction of the names of the output files (use initial density as identifier)
		out1 = "density_"+ identifier+".txt";
		out2 = "current_"+ identifier+".txt";


		// open the output file
		if(erecord.recording_density()){file1.open(out1.c_str());}
		if(erecord.recording_current()) {file2.open(out2.c_str());}

		// write the lattice parameters and particle properties into the output files
		if(erecord.recording_density()) {
			sys.write_parameter(file1);
			erecord.write_parameter(file1);
			erecord.write_density(file1,i);
		}

		if(erecord.recording_current()) {
			sys.write_parameter(file2);
			erecord.write_parameter(file2);
			erecord.write_current(file2, i);
		}

		// close the output files
		if(erecord.recording_density()) {file1.close();}
		if(erecord.recording_current()) {file2.close();}

	}

	// output in the console after finishing the simulation
	//cout << "Number of iteration steps:" << sys.get_count() << endl;
	//sys.print_state();

	// stop the time measurement
	control.stop();
	control.check();

	return 0;
}
