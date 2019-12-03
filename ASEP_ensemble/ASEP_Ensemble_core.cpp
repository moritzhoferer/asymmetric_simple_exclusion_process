// Author: Moritz Hoferer
// Date: 2015/09/08
// Further informations in the main.cpp

#include "ASEP_Ensemble_core.h"

using namespace std;

Lattice::Lattice(const double Parameters[22]){
	step_count = 0;
	t_first = 0.;
	t_now = 0.;
	time_on = false;
	action_error = false;

	//set random number generator
	r = gsl_rng_alloc(gsl_rng_mt19937);

	// set all parameters
	// system parameters
	ensemble			= Parameters[0];
	sites				= Parameters[1];
	boundary_cond		= bool(Parameters[2]);
	obc_option			= Parameters[3];
	empty_sys			= bool(Parameters[4]);
	rho_0				= Parameters[5];
	rho_1				= Parameters[6];

	// particle parameters
	r_right				= Parameters[7];
	r_left				= Parameters[8];
	r_push				= Parameters[9];
	r_repel				= Parameters[10];
	r_on				= Parameters[11];
	r_off				= Parameters[12];
	r_in_left			= Parameters[13];
	r_out_left			= Parameters[14];
	r_in_right			= Parameters[15];
	r_out_right			= Parameters[16];

	// seed
	seed				= Parameters[17];

	// measurement
	t_max				= Parameters[18];
	t_min				= Parameters[19];
	rec_density			= bool(Parameters[20]);
	rec_current			= bool(Parameters[21]);



	gsl_rng_set(r, seed);

	//introduce lattice of given size with initial density rho_0
	state.clear();
	state.resize(sites,false);

	// fill the lattice with "exact rounded" spatial initial density
	if(!empty_sys){
		fill_sys();
	}
	else{
		rho_0 = 0.0;
	}

	// normal obc w/ kicking out at the edge
	if((!boundary_cond) && (obc_option == 2)){
		r_repel_out_left = r_repel;
		r_push_out_right = r_push;
	}


	// pseudo pbc
	if((!boundary_cond) && (obc_option == 3)){
		r_in_left = (1 - rho_0) * rho_0 * r_right + rho_0 * rho_0 * r_push;
		r_out_left = (1- rho_0) * r_right;
		r_repel_out_left = (1 - rho_0) * r_repel;
		r_in_right = (1 - rho_1) * rho_1 * r_left + rho_1 * rho_1 * r_repel;
		r_out_right = (1 - rho_1) * r_left;
		r_push_out_right = (1 - rho_1) * r_push;
	}
}

Lattice::~Lattice(){
	//gsl_rng_free(r);

	//clear all arrays
	state.clear();
	p_right.clear();
	p_left.clear();
	p_push.clear();
	p_repel.clear();
	p_on.clear();
	p_off.clear();
	p_rates.clear();
	density.clear();
	current.clear();

}

void Lattice::reset(){
	// reset the system in a state that satisfies the initial conditions
	state.clear();
	state.resize(sites,false);
	if(!empty_sys){	fill_sys();}


	// reset step_count, t_first, t_now, time_on, action_error
	step_count = 0;
	t_first = 0.;
	t_now = 0.;
	time_on = false;
	action_error = false;

}

// fill the system with the "exact" initial spatial density rho_0
void Lattice::fill_sys(){
	urn.clear();
	urn.resize(sites/2 , 0);
	for(unsigned int i=0;i < sites/2;i++){urn[i] = i;}

	// fill_stop:	number of particle to fill into the system
	int fill_stop_0 = trunc(rho_0 * sites/2);

	for(unsigned int i = 0; i < fill_stop_0; i++ ){
		// select one site out of the urn
		position = gsl_rng_uniform_int(r,urn.size());
		// set the selected site to "True"
		state[urn[position]] = true;
		// remove the selected site out of the urn
		urn.erase(urn.begin() + position);
	}

	urn.resize(sites/2 , 0);
	for(unsigned int i = 0; i < sites/2; i++){urn[i] = sites/2  + i;}

	int fill_stop_1 = trunc(rho_1 * sites/2);

	for(unsigned int i = 0; i < fill_stop_1; i++ ){
			// select one site out of the urn
			position = gsl_rng_uniform_int(r,urn.size());
			// set the selected site to "True"
			state[urn[position]] = true;
			// remove the selected site out of the urn
			urn.erase(urn.begin() + position);
		}
}

// Includes all possible action with pbc and obc
void Lattice::list_actions(){
	// Reset all arrays respectively set the boundary values to false (not possible)
	p_right.clear();
	p_left.clear();
	p_push.clear();
	p_repel.clear();
	p_on.clear();
	p_off.clear();

	// periodic boundary conditions
	if(boundary_cond){
		// first site is occupied
		if(state[0]){
			// if second site vacant
			if(!state[1]){
				// if last site vacant -> hopping right
				if(!state[sites-1]){p_right.push_back(0);}
				// if last site occupied -> being pushed
				else{p_push.push_back(0);}
			}
			// if last site is vacant
			if(!state[sites-1]){
				// if second site is vacant -> hopping left
				if(!state[1]){p_left.push_back(0);}
				// if second site is occupied -> being repelled
				else{p_repel.push_back(0);}
			}
		}

		// last site is occupied
		if(state[sites-1]){
			// if first site is vacant
			if(!state[0]){
				// if second last site is vacant -> hopping right
				if(!state[sites -2]){p_right.push_back(sites-1);}
				// if second last site occupied -> being pushed
				else{p_push.push_back(sites-1);}
			}
			// if second last site is vacant
			if(!state[sites-2]){
				// if second last site is empty -> hopping left
				if(!state[0]){p_left.push_back(sites-1);}
				// if fist site is occupied -> being repelled
				else{p_repel.push_back(sites-1);}
			}
		}
	}
	// open boundary conditions
	else{
		// if first site is occupied and second not -> hopping right
		if((state[0]) && (!state[1])){p_right.push_back(0);}


		// if  last sites is occupied and second last not -> hopping left
		if((state[sites -1])&& (!state[sites-2])){p_left.push_back(sites-1);}
	}

	//second to second last site: hopping left and right (independent of the boundary conditions)
	for(unsigned int i = 1; i < sites-1;i++){
		// occupied site
		if(state[i]){
			// right neighbor is vacant
			if(!state[i+1]){
				// left neighbor is vacant -> hopping right
				if(!state[i-1]){p_right.push_back(i);}
				// left neighbor is occupied -> being pushed
				else{p_push.push_back(i);}
			}
			// left neighbor is vacant
			if(!state[i-1]){
				// right neighbor is vacant -> hopping left
				if(!state[i+1]){p_left.push_back(i);}
				// right neighbor is empty -> being repelled
				else{p_repel.push_back(i);}
			}
		}
	}

	// first to last site: attachment and detachment (independent of the boundary conditions)
	for(unsigned int i = 0; i < sites; i++){
		// occupied site -> detach
		if(state[i]){p_off.push_back(i);}
		// vacant site -> attach
		else{p_on.push_back(i);}
	}
}

void Lattice::make_rates(){
	p_rates.clear();
	// content order of p_rates: ('p_right',		'p_left',		'p_push',		'p_repel',		'p_on',
	//							  'p_off',			'p_in_left',	'p_out_left',	'p_in_right',	'p_out_right')
	p_rates.push_back(p_right.size() * r_right);
	p_rates.push_back(p_rates.back() + p_left.size() * r_left);
	p_rates.push_back(p_rates.back() + p_push.size() * r_push);
	p_rates.push_back(p_rates.back() + p_repel.size() * r_repel);
	p_rates.push_back(p_rates.back() + p_on.size() * r_on);
	p_rates.push_back(p_rates.back() + p_off.size() * r_off);
	if(!boundary_cond){
		// obc w/o kicking out at the edge
		if(obc_option == 1){
		p_rates.push_back(p_rates.back() + (!state[0]) * r_in_left);
		p_rates.push_back(p_rates.back() + (state[0])  * r_out_left);
		p_rates.push_back(p_rates.back() + (!state[sites-1]) * r_in_right);
		p_rates.push_back(p_rates.back() + (state[sites-1]) * r_out_right);
		}
		// obc w/ kicking out at the edge
		else if(obc_option == 2){
			p_rates.push_back(p_rates.back() + (!state[0]) * r_in_left);
			p_rates.push_back(p_rates.back() + (state[0]) * ( r_out_left +  (state[1]) * r_repel_out_left) );
			p_rates.push_back(p_rates.back() + (!state[sites-1]) * r_in_right);
			p_rates.push_back(p_rates.back() + (state[sites-1]) * ( r_out_right + (state[sites -2]) * r_push_out_right));
		}
		// pseudo pbc
		else if(obc_option == 3){
			p_rates.push_back(p_rates.back() + (!state[0]) * r_in_left  );
			p_rates.push_back(p_rates.back() + (state[0])  * ( (!state[1]) * r_out_left + (state[1]) * r_repel_out_left));
			p_rates.push_back(p_rates.back() + (!state[sites-1]) * r_in_right);
			p_rates.push_back(p_rates.back() + (state[sites-1]) * ((!state[sites -2]) * r_out_right + (state[sites-2]) * r_push_out_right ));
		}
	}
}



// This should work if the list_actions-function distinguishes right between pbc and obc
// Also included is the measurement of the current
void Lattice::hop_right(){
	// choose with uniform distribution  the one to move
	position = gsl_rng_uniform_int(r,p_right.size());

	// move the chosen one one site to the right
	if(p_right[position] != sites-1){
		state[p_right[position]] = false;state[p_right[position] + 1] = true;
		// check if current should be recorded and if t_before > t_min
	//	if((rec_current) && (time_on)){current[p_right[position]+1] += 1;}
	}
	else{
		state[sites-1] = false;state[0] = true;
	//	if((rec_current) && (time_on)){current[0] += 1;}
	}
}

void Lattice::hop_left(){
	// choose randomly one to move
	position = gsl_rng_uniform_int(r,p_left.size());

	if(p_left[position] != 0){
		// move the chosen one one site to the left
		state[p_left[position]]=false; state[p_left[position] -1] = true;
		// record movement for the current measurement
	//	if((rec_current) && (time_on)){current[p_left[position]] -= 1;}
	}
	else{state[0] = false; state[sites-1] = true;
	//	if((rec_current) && (time_on)){current[sites-1] -= 1;}
	}
}

void Lattice::be_pushed(){
	position = gsl_rng_uniform_int(r,p_push.size());

	if(p_push[position] != (sites-1)){
		state[p_push[position]] = false;	state[p_push[position] + 1] = true;
	//	if((rec_current) && (time_on)){current[p_push[position]+1] += 1;}
	}
	else{
		state[sites-1] = false;state[0] = true;
	//	if((rec_current) && (time_on)){current[0] += 1;}
	}
}

void Lattice::be_repelled(){
	position = gsl_rng_uniform_int(r,p_repel.size());

	if(p_repel[position] != 0){
		state[p_repel[position]]=false; state[p_repel[position] -1] = true;
	//	if((rec_current) && (time_on)){current[p_repel[position]] -= 1;}
	}
	else{
		state[0] = false; state[sites-1] = true;
	//	if((rec_current) && (time_on)){current[p_repel[position]] -= 1;}
	}
}

void Lattice::attach(){
	position = gsl_rng_uniform_int(r,p_on.size());
	state[p_on[position]] = true;
}

void Lattice::detach(){
	position = gsl_rng_uniform_int(r,p_off.size());
	state[p_off[position]] = false;
}

void Lattice::in_left(){
	state[0] = true;
	// particle is entering the lattice on the left end, means particle is moving to the right
	//if((rec_current) && (time_on)){current[0] += 1;}
}

void Lattice::out_left(){
	state[0] = false;
	// particle is leaving the lattice on the left end, means particle is moving to the left
	//if((rec_current) && (time_on)){current[0] -= 1;}
}

void Lattice::in_right(){
	state[sites-1] = true;
	// particle is entering the lattice on the right end, means particle is moving left
	//if((rec_current) && (time_on)){current[sites] -= 1;}
}

void Lattice::out_right(){
	state[sites-1] = false;
	// particle is leaving the lattice on the left end, means particle is moving right
	//if((rec_current) && (time_on)){current[sites] += 1;}
}

void Lattice::simulate_time(){
	list_actions();
	make_rates();

	if(p_rates.back() != 0.0){

		// update time
		t_before = t_now;
		t_now -= log(gsl_rng_uniform_pos(r))/p_rates.back();

		step_count++;

		// if run the fist time this loop save time as t_first
		if((t_before > t_min) && (!time_on)){
			t_first = t_before;
			time_on = true;
		}

	}

	else{
		// activate the action error
		action_error = true;
		cout << "Simulation has been truncated after t = "<< t_now <<",  because there are no more possible actions."<< endl;
	}

}

void Lattice::simulate_movement(){

		// record density if rec_density and time_on are True
		//if((rec_density) && (time_on)){record_density();}

		// select the action to execute
		// THINK ABOUT WHICH IS THE RIGTH RANDOM COMAND (w/ or w/o 0. and 1.?!)
		double select = gsl_rng_uniform(r) * p_rates.back();

		// content order of p_rates: ('p_right','p_left','p_jostle','p_repel','p_on','p_off',
		// 'p_in_left',p_out_left','p_in_right','p_out_right')
		// Alternatively use the BinarySearch algorithm
		if(select <= p_rates[0]){hop_right();}
		else if(select <= p_rates[1]){hop_left();}
		else if(select <= p_rates[2]){be_pushed();}
		else if(select <= p_rates[3]){be_repelled();}
		else if(select <= p_rates[4]){attach();}
		else if(select <= p_rates[5]){detach();}
		else if(select <= p_rates[6]){in_left();}
		else if(select <= p_rates[7]){out_left();}
		else if(select <= p_rates[8]){in_right();}
		else if(select <= p_rates[9]){out_right();}


}

void Lattice::print_state() const{
	cout << "Time:" << t_now <<"\t Current state:";
	for(unsigned int i = 0; i < sites; i++ ){
		if(state[i]){cout << "1";}
		else{cout << "0";}
	}
	cout << endl;
}

void Lattice::get_final_rho(){
	double c = 0.0;
	// count the particles in the lattice
	for(unsigned int i; i< sites; i++){c += 1.0 * (state[i]);}
	rho_final = c / sites;
	//cout << "rho_final: " << rho_final << endl;

}

long int Lattice::get_count() const{
	return step_count;
}

double Lattice::get_time() const{
	return t_now;
}

// give 1 if site d is occupied and gives 0 if site d is vacant
bool Lattice::get_single_state(int &i){
	return state[i];
}

bool Lattice::using_pbc() const{
	return boundary_cond;
}

void Lattice::write_parameter(std::ofstream &write) const{
	write << "# Simulation of a system with the following properties and parameters:" << endl;
	write << "#"<< endl;
	write << "# System configuration" << endl;
	write << "# Number of systems:			" << ensemble << endl;
	write << "# System size:					" << sites << endl;
	write << "# Boundary conditions:			";
	if(boundary_cond){write << "periodic" << endl;}else{write <<"open"<< endl;};
	write << "# Start with empty system:		";
	if(empty_sys){write << "yes" << endl;}
	else{write << "no" << endl;}
	write << "# Initial spatial density:		" << rho_0 << endl ;
	write << "# Final spatial density:		" << rho_final << endl;
	write << "#" << endl;
	write << "# Particle properties"<< endl;
	write << "# Rates of the different actions:" << endl;
	write << "# Hop right:					" << r_right << endl;
	write << "# Hop left:						" << r_left << endl;
	write << "# Been pushed:					" << r_push << endl;
	write << "# Been repelled:				" << r_repel << endl;
	write << "# Attach:						" << r_on << endl;
	write << "# Detach:						" << r_off << endl;
	if(!boundary_cond){
	write << "# Enter from the left:			" << r_in_left << endl;
	write << "# Leave at the left:			" << r_out_left << endl;
	write << "# Enter from the right:			" << r_in_right << endl;
	write << "# Leave at the right:			" << r_out_right << endl;
	}
	write << "#"<< endl;
	write << "# Seed:							" << seed << endl;
	write << "#"<< endl;
}

void Lattice::write_state(std::ofstream &write) const{
	write << "State of the lattice at time(" << t_now << "): \n";
	for(unsigned int i = 0; i < state.size(); i++){
		if(state[i]){write << "1";} else{write << "0";}
	}
	write << endl;
}




Measurement::Measurement(const double Parameters[22]){
	// system
	ensemble		= Parameters[0];
	sites		  	= Parameters[1];
	boundary_cond 	= bool(Parameters[2]);
	obc_option		= Parameters[3];

	// particles
	t_max			= Parameters[18];
	delta_t			= Parameters[19];
	rec_current		= bool(Parameters[21]);
	rec_density		= bool(Parameters[20]);



	// this is used to fill the deque with the time points
	double t = 0.;

	while(t <= t_max ){
		record_times.push_back(t);
		t += delta_t;
	}

	// construct the current deque depending on the boundary conditions
	if(boundary_cond){
		current.resize(sites*record_times.size(),0.);
	}
	else{
		current.resize((sites+1)*record_times.size(),0.);
	}

	density.resize(sites*record_times.size(),0);

}

Measurement::~Measurement(){
	record_times.clear();
	density.clear();
	current.clear();
}

double Measurement::get_time_point(int &d) const{
	return record_times[d];
}

int Measurement::get_ensemble() const{
	return ensemble;
}

int Measurement::get_numb_time_points() const{
	return floor(t_max/delta_t);
}


void Measurement::write_parameter(std::ofstream &write) const{
	write << "# Simulation time t =			" << t_max << endl;
	write << "# Time between snaps shots dt=	" << delta_t << endl;
	/*write << "# Number of simulated steps:	" << cl.get_count() << endl;
	write << "# Simulation ended because:		";
	if(action_error){write << "no more possible actions" << endl;}
	else if(end_after_tmax){write << "t_max was reached"<< endl;}
	else{write << "iteration_max was reached"<< endl;}*/
	write << "# Recorded quantities:			";
	if((rec_density) && (rec_current)){write << "density and current "<<endl;}
	else if(rec_density){write << "density" << endl;}
	else if(rec_current){write << "current" << endl;}
	else {write << "Nothing was recorded!!!" << endl;}
	write << "#" << endl;
}

bool Measurement::not_finished(class Lattice &cl) const{
	if((t_max > cl.get_time()) /*&& (!g)*/){return true;}
	else{/*cout <<"end now!"<< endl;*/ return false;}
}

bool Measurement::recording_density() const{
	return rec_density;
}

bool Measurement::recording_current() const{
	return rec_current;
}

void Measurement::record(class Lattice &cl, int &d){
	if(rec_density){
		for(int i = 0; i < sites; i++){
			density[d* sites + i] +=  1.0/ensemble * (cl.get_single_state(i));
		};
	}
	// introduce in class Lattice a variable signs the position of the last movement with direction
	if(rec_current){
		cout << d;
	}
}



void Measurement::write_density(std::ofstream &write, int &d) const {
	write << "# Density profile in the system at time(" << record_times[d] << "):" << endl;
	for(unsigned int i = 0; i < sites; i++){
		write << density[d*sites + i] << endl;
	}
}

void Measurement::write_current(std::ofstream &write, int &d) const {

	write << "# Current in the system at time(" << record_times[d] << "):" << endl;
	// pbc
	if(boundary_cond){
		for(signed int i = 0; i < sites; i++){
			/*write << current[d*sites + i] << endl*/;
		}
	}
	// obc
	else{
		for(signed int i = 0; i < sites + 1; i++){
				/*write << current[d*(sites+1) + i] << endl*/;
		}
	}
}

// used to convert a number (int, double, etc.) into a string
string do_to_str(const double &d){
	stringstream ss;
	ss << d;
	return(ss.str());
}

// used to read in the configuration parameters out of the parameters[i].txt-file
void read_in_parameters(const string &str, double Params[]){
	ifstream in(str.c_str());
	double data;
	int i = 0;


	while(in >> data){
		Params[i]=data;
		i++;
	}
	in.close();
}
