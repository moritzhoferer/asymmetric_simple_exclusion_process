// Author: Moritz Hoferer
// Date: 2015/09/08
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

	// ensemble:	number of systems
	// sites:		number of sites in the system
	// position:	is used for selecting the right site out of the possibility arrays
	// obc_option:	(1) normal obc w/o kicking out at the edge
	//				(2) normal obc w/  kicking out at the edge
	//				(3) pseudo pbc
	// seed:		seed for random function
	int ensemble, sites, position, obc_option, seed;

	// step_count:		count iteration steps
	// iteration_max:	maximum number of iteration steps
	long int step_count, iteration_max;

	// RATES:
	// r_ right:		rate to jump to the right
	// r_left:			rate to jump to the left
	// r_on:			rate to attach on one site
	// r_off:			rate to detach from one site
	// r_repel:			being repelled rate
	// r_push:			being pushed rate
	// r_in_left
	// r_out_left
	// r_repel_out_lef:
	// r_in_right:		rate to enter the system from the right (alpha) just for the open system
	// r_out_right:		rate to leave the system from the last site (beta) just for the open system
	// r_push_out_right:
	// tot_rate:		sum over the rates of all possible next actions
	// rho_0: 			initial density in the left half of the lattice
	// rho_1:			initial density in the right half of the lattice
	// rho_final: 		to check if particle number is conserved with pbc and without attachment and detachment
	//					to check evolution of density with obc and/or attachment and detachment
	double r_right, r_left, r_on, r_off, r_repel, r_push, r_in_left, r_in_right, r_out_left, r_repel_out_left , r_out_right, r_push_out_right, tot_rate, rho_0,rho_1, rho_final;

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
	Lattice(const double Parameters[22]);

	//destructor
	~Lattice();

	// fill the system with the "exact" initial spatial density rho_0 in the left half and with rho_1 in the right half of the lattice
	void fill_sys();

	// fill the array with the numbers of the sites which can make a action
	void list_actions();

	// build the array with the needed rates for different actions
	void make_rates();

	// function that resets the lattice for a new run
	void reset();

	// simulates the time progress of the Gillespie step
	void simulate_time();

	// simulates the movement of the Gillespie step
	void simulate_movement();

	// print the current state of the system with length and time
	void print_state() const;

	// TRANSFER THIS TO MEASUREMENT CLASS
	// print the final lattice density
	void get_final_rho();

	// print the number of simulated Gillespie steps
	long int get_count() const;

	// print current simulation time
	double get_time() const;

	// give 1 if site d is occupied and gives 0 if site d is vacant
	bool get_single_state(int &i);

	// True if system is a closed ring
	bool using_pbc() const;

		// writes state of the lattice and time in a file
	void write_state(std::ofstream &write) const;

	// writes the used parameters of the system in a file
	void write_parameter(std::ofstream &write) const;

	};

// include this further class in a second step
class Measurement{
private:
	// boundary condition to check the length of the current deque
	// boundary_cond: 	pbc (1) OR obc (0)
	// rec_current:		yes(1) OR no(0)
	// rec_density: 	yes(1) OR no(0)
	bool boundary_cond, rec_current, rec_density ;

	// emsemble:	number of lattices in the ensemble
	// sites:		number of sites of the lattice
	int ensemble, sites, obc_option;

	// t_max:		maximum simulation time
	// delta_t:		time between the recording time points
	double t_max, delta_t;

	// time:		record time points
	// density:		recorded density profile at the time points
	// current:		recorded current profile at the time points
	std::deque<double> record_times;
	std::deque<double> density;
	std::deque<double> current;



public:
	// constructor
	Measurement(const double Parameters[22]);

	// destructor
	~Measurement();

	// gives the next time point to record
	double get_time_point(int &d) const;

	// gives number of simulated systems
	int get_ensemble() const;

	// gives number of simulated time point plus 1 for the initial state
	int get_numb_time_points() const;

	// check if number of simulation steps and/or time is not over the maximal values. True none of the maxima is exceeded
	bool not_finished(class Lattice &cl) const;

	// True if the density is recorded
	bool recording_density() const;

	// True if the current is recorded
	bool recording_current() const;

	// checks what should be recorded and records
	void record(class Lattice &cl, int &d);

	// writes the used parameters of the measurement in a file
	void write_parameter(std::ofstream &write) const;

	// checks if the density was recorded and writes it in the output file
	void write_density(std::ofstream &write,int &d) const;

	// checks if the current was record and writes it in the out file
	void write_current(std::ofstream &write,int &d) const;

};

// used to convert a number (int, double, etc.) into a string
std::string do_to_str(const double &d);

// used to read in the configuration parameters out of the parameters[i].txt-file
void read_in_parameters(const std::string &str, double Params[]);
