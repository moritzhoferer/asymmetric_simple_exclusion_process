// Author: Moritz Hoferer
// Date: 2015/08/07
// Further informations in the main.cpp

#include <cmath>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <deque>

class Lattice{
private:

	// True VS False
	// boundary_cond:	periodic boundary conditions (pbc) VS open boundary conditions (obc)
	// empty_sys:		start with empty system VS finite density in the system
	bool boundary_cond, empty_sys;

	// end_after_tmax:		maximum simulation time
	// action_error:		True if there is no more possible action in the system for a next Gillespie step
	// time_on:				True if t_now bigger than t_min
	// record the current:	True if current is recorded
	// record the density:	True if average density of every single site is recorded
	bool end_after_tmax, action_error , time_on, rec_current, rec_density;

	// sites:		number of sites in the system
	// position:	is used for selecting the right site out of the possibility arrays
	// seed:		seed for random function
	int sites, position, seed;

	// step_count:		count iteration steps
	// iteration_max:	maximum number of iteration steps
	long int step_count, iteration_max;

	// RATES:
	// r_ right:	rate to jump to the right
	// r_left:		rate to jump to the left
	// r_on:		rate to attach on one site
	// r_off:		rate to detach from one site
	// r_repel:		being repelled rate
	// r_push:		being pushed rate
	// r_in_right:	rate to enter the system from the right (alpha) just for the open system
	// r_out_right:	rate to leave the system from the last site (beta) just for the open system
	// tot_rate:	sum over the rates of all possible next actions
	// rho_0: 		initial density
	// rho_final: 	to check if particle number is conserved with pbc and without attachment and detachment
	//				to check evolution of density with obc and/or attachment and detachment
	double r_right, r_left, r_on, r_off, r_repel, r_push, r_in_left, r_in_right, r_out_left, r_out_right, tot_rate, rho_0, rho_final;

	// TIMES:
	// t_now: 		simulation time
	// t_before: 	the Gillespie step
	// t_min: 		time to start recording
	// t_max: 		maximum simulation time
	// t_first: 	first time value that is bigger t_min
	double t_now, t_before, t_min, t_max, t_first ;


	//gsl random number generator
	gsl_rng * r;

	// ARRAYS
	// state: 		occupied (True) or not occupied (False) site
	// urn:			used to fill the system with the right initial density rho_0
	// density:		record the occupation time to divide it in the end by the record time -> average density of every site
	// current:		count the transitions between neighboring sites. Divide by the record time -> average current between all sites
	std::deque<bool> state;
	std::deque<int> urn;
	std::deque<double> density;
	std::deque<double> current;

	// POSSIBILITIES
	// p_rigth: 	sites with particle, that can hop to the right
	// p_left:		sites with particle, that can hop to the left
	// p_push:		sites with particle, that can be jostled by a particle from behind
	// p_ repel:	sites with particle, that can be repelled by the a particle from the front
	// p_on:		not occupied sites for the possibility, that a particle attaches
	// p_off		occupied sites for the possibility, that a particle detaches
	std::deque<unsigned int> p_right;
	std::deque<unsigned int> p_left;
	std::deque<unsigned int> p_push;
	std::deque<unsigned int> p_repel;
	std::deque<unsigned int> p_on;
	std::deque<unsigned int> p_off;
	std::deque<double> p_rates;

	// possible actions
	void hop_right();
	void hop_left();
	void be_pushed();
	void be_repelled();
	void attach();
	void detach();
	void in_left();
	void out_left();
	void in_right();
	void out_right();


public:

	// initialize Lattice
	// insert number of the relevant parameters
	Lattice(const double Parameters[21]);

	//destructor
	~Lattice();

	// fill the system with the "exact" initial spatial density rho_0
	void fill_sys();

	// fill the array with the numbers of the sites which can make a action
	void list_actions();

	// build the array with the needed rates for different actions
	void make_rates();


	// simulate the next Gillespie step
	void simulate();

	// adds up the occupation times of the system
	void record_density();

	// print the current state of the system with length and time
	void print_state() const;

	// print the density
	void print_density() const;

	// print the current
	void print_current() const;

	// print the final lattice density
	void get_final_rho();

	// print the number of simulated Gillespie steps
	long int get_count() const;

	// print current simulation time
	double get_time() const;

	// check if number of simulation steps and/or time is not over the maximal values. True none of the maxima is exceeded
	bool not_finished() const;


	// True if the density is recorded
	bool recording_density() const;

	// True if the current is recorded
	bool recording_current() const;

	// True if system is a closed ring
	bool using_pbc() const;

		// writes state of the lattice and time in a file
	void write_state(std::ofstream &write) const;

	// norm the density profile and write it into write
	void write_density(std::ofstream &write) const;

	// divide current entries by the time
	// particle moving to the right is defined as positive current.
	void write_current(std::ofstream &write) const;

	// writes the used parameters in a file
	void write_parameter(std::ofstream &write) const;

	};


// used to convert a number (int, double, etc.) into a string
std::string do_to_str(const double &d);

// used to read in the configuration parameters out of the parameters[i].txt-file
void read_in_parameters(const std::string &str, double Params[]);
