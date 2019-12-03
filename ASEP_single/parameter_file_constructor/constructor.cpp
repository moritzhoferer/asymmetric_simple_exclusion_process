//============================================================================
// Name        : constructor_short.cpp
// Author      : Moritz Hoferer
// Date		   : 2015/08/26
// Version     : 1.0
// Copyright   : n/a
// Description : Constructs a set of parameter files
//============================================================================

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <fstream>
#include <math.h>
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
	gsl_rng_set(r, 13);

	int i, sites, end_after_tmax, rec_current, rec_density, boundary_cond, empty_sys, systems, steps1, steps2;
	long int iteration_max, seed;
	double t_max, t_min, r_right, r_left, r_push, r_repel, r_on, r_off, r_in_left, r_out_left, r_in_right, r_out_right, rho_0, min1, min2, max1, max2, delta1, delta2;

	ofstream file,file2;
	string out,out2;

	// system configuration
	sites 			= 100000;
	boundary_cond	= 0; 			// pbc (1) OR opc (0)
	empty_sys		= 1;			// empty system (1) OR system with density rho_0 (0)
	rho_0			= 0.0;
	systems			= 1;			// number of system, that are simulated parallel


	// particle properties
	r_right			= 1.0;
	r_left			= 0.0;
	r_push			= 1.0;
	r_repel			= 2.0;
	r_on			= 0.0;
	r_off			= 0.0;
	r_in_left		= 0.0;
	r_out_left		= 0.0;
	r_in_right		= 0.0;
	r_out_right		= 0.0;

	// simulation adjustment
	// yes (1) OR no (0)
	end_after_tmax	= 1;		// ending after t_max (1) OR ending after iteration_max (0)
	rec_current		= 1;		// record on (1) OR record off (0)
	rec_density		= 1;		// record on (1) OR record off (0)
	iteration_max	= 10000;
	t_max			= 41000;
	t_min			= 1000;

	i = 0;

	// first variable parameter
	// r_repel
	min1 	= 0.0;
	max1 	= 1.0;
	steps1 	= 20;
	delta1 	= (max1 - min1) / steps1;

	// second variable parameter
	// r_in_left
	min2 	= 0.0;
	max2	= 1.0;
	steps2	= 20;
	delta2 = (max2 -min2) / steps2;

	out2 = "option.txt";
	file2.open(out2.c_str());

	double c_concentr[3] = {10,50,250};
	double c_repel[3] = {13 , 20, 50};

	//for loop first variable parameter
	for (r_out_right = min1; r_out_right < max1 + delta1; r_out_right += delta1){
	//for(int m = 0; m<3; m++){

		// for loop second variable parameter
		for (r_in_left = min2; r_in_left < max2 + delta2 ; r_in_left += delta2){
		//for(int n = 0; n<3; n++){

			r_out_left = 1-r_in_left;

			if(r_out_left < 0){
				r_out_left = 0.0;
			}

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

				file << sites				<< endl;
				file << t_max				<< endl;
				file << t_min				<< endl;
				file << iteration_max		<< endl;
				file << end_after_tmax		<< endl;
				file << rec_current			<< endl;
				file << rec_density			<< endl;
				file << boundary_cond		<< endl;
				file << empty_sys			<< endl;
				file << r_right				<< endl;
				file << r_left				<< endl;
				file << r_push 				<< endl;
				//file << c_repel[m] 			<< endl;
				//file << r_on * c_concentr[n]				<< endl;
				file << r_repel 			<< endl;
				file << r_on 				<< endl;
				file << r_off 				<< endl;
				file << r_in_left 			<< endl;
				file << r_out_left 			<< endl;
				file << r_in_right 			<< endl;
				file << r_out_right 		<< endl;
				file << rho_0 				<< endl;
				file << seed 				<< endl;

				file.close();
			}

		}

	}

	file2.close();

	return 0;
}
