//============================================================================
// Name        : constructor_Ensemble.cpp
// Author      : Moritz Hoferer
// Date		   : 2015/09/12
// Version     : 1.0
// Copyright   : n/a
// Description : Constructs a set of parameter files for the ASEP_Ensemble
//============================================================================

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

string do_to_str(const double &d){
	stringstream ss;
	ss << d;
	return(ss.str());
}

int main(){

	//gsl random number generator
	gsl_rng * r;
	r = gsl_rng_alloc(gsl_rng_mt19937);
	gsl_rng_set(r, 8);

	int i, ensemble, sites, rec_current, rec_density, boundary_cond, obc_option, empty_sys, systems, steps1, steps2;
	long int seed;
	double t_max, delta_t, r_right, r_left, r_push, r_repel, r_on, r_off, r_in_left, r_out_left, r_in_right, r_out_right, rho_0, rho_1, min1, min2, max1, max2, delta1, delta2;

	ofstream file,file2;
	string out,out2;

	// system configuration
	ensemble		= 1000;			// number of systems in the ensemble
	sites 			= 3000;			// number of sites in the system
	boundary_cond	= 1; 			// pbc (1) OR opc (0)
	obc_option 		= 0;			// obc_option:	(1) normal obc w/o kicking out at the edge
									//				(2) normal obc w/  kicking out at the edge
									//				(3) pseudo pbc
	empty_sys		= 0;			// empty system (1) OR system with density rho_0 and rho_1 (0)
	rho_0			= .125;			// density in the left half of the lattice
	rho_1			= .375;			// density in the right half the lattice
	systems			= 1;			// number of system, that are simulated parallel (NOT ENSEMBLE)
									// Creates several times the same file with different seed


	// particle properties
	r_right			= 1.0;
	r_left			= 0.0;
	r_push			= 1.0;
	r_repel			= 0.0;
	r_on			= 0.0;
	r_off			= 0.0;
	r_in_left		= 0.0;
	r_out_left		= 0.0;
	r_in_right		= 0.0;
	r_out_right		= 0.0;

	// simulation adjustment
	// yes (1) OR no (0)
	rec_current		= 0;		// record on (1) OR record off (0)
	rec_density		= 1;		// record on (1) OR record off (0)
	t_max			= 300;
	delta_t			= 300;

	i = 0;

	// first variable parameter
	// rho_0
	min1 	= 0.25;
	max1 	= 0.75;
	steps1 	= 2;
	delta1 	= (max1 - min1) / steps1;

	// second variable parameter
	// rho_1
	min2 	= 0.25;
	max2	= 0.75;
	steps2	= 2;
	delta2 = (max2 -min2) / steps2;

	out2 = "option.txt";
	file2.open(out2.c_str());

	// for loop second variable parameter
	for (rho_1 = min2; rho_1 <= max2 ; rho_1+= delta2){

		// for loop first variable parameter
		for (rho_0 = min1; rho_0 <= max1 ; rho_0 += delta1){

			// number of system that are simulated parallel.
			// Every system gets another seed
			for(int j = 0; j < systems; j++){

				i++;
				file2 << i << endl;

				seed = gsl_rng_uniform_int(r,1000);

				out = "Parameters";
				out.append(do_to_str(i));
				out.append(".txt");

				file.open(out.c_str());

				// system
				file << ensemble 			<< endl;		//  0
				file << sites				<< endl;		//  1
				file << boundary_cond		<< endl;		//  2
				file << obc_option			<< endl;		//  3
				file << empty_sys			<< endl;		//  4
				file << rho_0 				<< endl;		//  5
				file << rho_1				<< endl;		//  6

				// particles
				file << r_right				<< endl;		//  7
				file << r_left				<< endl;		//  8
				file << r_push 				<< endl;		//  9
				file << r_repel 			<< endl;		// 10
				file << r_on 				<< endl;		// 11
				file << r_off 				<< endl;		// 12
				file << r_in_left 			<< endl;		// 13
				file << r_out_left 			<< endl;		// 14
				file << r_in_right 			<< endl;		// 15
				file << r_out_right 		<< endl;		// 16

				// seed
				file << seed 				<< endl;		// 17

				// measurement
				file << t_max				<< endl;		// 18
				file << delta_t				<< endl;		// 19
				file << rec_density			<< endl;		// 20
				file << rec_current			<< endl;		// 21


				file.close();
			}

		}

	}

	file2.close();

	return 0;
}
