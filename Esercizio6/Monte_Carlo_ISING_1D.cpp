/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main(){
 
	Input(); //Inizialization
	
	if (equilibration) { // This only checks the energy to see equilibration
		cout << "Checking the system's energy for equilibration" << endl;
		Equilibration();	
	}
	else {
		for (int i = 0; i < ntemp; i++) { // Cicle over different temperatures
			for (int j = 0; j < step_equi; j++) 
				Move(metro); // Reach equilibration
			
			for (int iblk = 1; iblk <= nblk; iblk++) { //Simulation
				Reset(iblk);   //Reset block averages
				for(int istep = 1; istep <= nstep; istep++){
					Move(metro);
					Measure();
					Accumulate(); //Update block averages
				}
				Averages(iblk);   //Add to sum and sum squared
			}
			Results();	//Print result with error
			cout << "Finished temperature T = " << 0.5+1.5*double(i)/double(ntemp-1) << endl;
			beta = 1. / (0.5 + 1.5 * double(i+1) / double(ntemp-1)); // Change temperature
		}
	}
	ConfFinal(); //Write final configuration

	return 0;
}


void Input(void){

	ifstream ReadInput;

	cout << "Classic 1D Ising model             " << endl;
	cout << "Monte Carlo simulation             " << endl << endl;
	cout << "Nearest neighbour interaction      " << endl << endl;
	cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
	cout << "The program uses k_B=1 and mu_B=1 units " << endl;

	//Read seed for random numbers
	int p1, p2;
	ifstream Primes("Primes");
	Primes >> p1 >> p2 ;
	Primes.close();

	ifstream input("seed.in");
	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
	rnd.SetRandom(seed,p1,p2);
	input.close();
  
	//Read input informations
	ReadInput.open("input.dat");

	ReadInput >> temp;
	beta = 1.0/temp;
	cout << "Temperature = " << temp << endl;

	ReadInput >> nspin;
	cout << "Number of spins = " << nspin << endl;

	ReadInput >> J;
	cout << "Exchange interaction = " << J << endl;

	ReadInput >> h;
	cout << "External field = " << h << endl << endl;
    
	ReadInput >> metro; // if=1 Metropolis else Gibbs

	ReadInput >> nblk;

	ReadInput >> nstep;
	
	ReadInput >> equilibration; // if = 1 only perform equilibration
	
	ReadInput >> conf_final; // if = 1 read previous configuration
	
	ReadInput >> ntemp; // number of temperatures to sample in [0.5, 2]

	if(metro==1) cout << "The program perform Metropolis moves" << endl;
	else cout << "The program perform Gibbs moves" << endl;
	
	cout << "Number of blocks = " << nblk << endl;
	cout << "Number of steps in one block = " << nstep << endl << endl;
	ReadInput.close();

	//Prepare arrays for measurements
	iu = 0; //Energy
	ic = 1; //Heat capacity
	im = 2; //Magnetization
	ix = 3; //Magnetic susceptibility
 
	n_props = 4; //Number of observables

	//Generate the initial configuration
	if (conf_final) {
		cout << "Reading previous configuration" << endl << endl;
		ifstream config;
		config.open("config.final");
		for (int i=0; i<nspin; i++){
			config >> s[i];
		}
		config.close();
	}
	else {
		cout << "Generating a new configuration" << endl << endl;
		for (int i=0; i<nspin; i++){
			if(rnd.Rannyu() >= 0.5) s[i] = 1;
			else s[i] = -1;
		}
	} 
}

void Move(int metro){

	int o;
	double prob;

	for(int i=0; i<nspin; i++){
	
		//Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
		o = (int)(rnd.Rannyu()*nspin);
		attempted ++;
		if(metro){ //Metropolis
			prob = Boltzmann(s[o], o);
			if (prob > rnd.Rannyu()){
				s[o] *= -1;
				accepted ++;
			}
		}
		else{ //Gibbs sampling
			prob = 1. / (1 + Boltzmann(1., o));
			accepted ++;
			if (prob > rnd.Rannyu())
				s[o] = 1;
			else
				s[o] = -1;
		}
	}
}

double Boltzmann(int sm, int ip){

	//Calculate the probability for the Metropolis algorithm
	double ene = J * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) + h;
	return exp(- 2. * beta * ene * sm);
}

void Measure(){

	double u = 0.0, m = 0.0;

	//Cycle over spins
	for (int i=0; i<nspin; i++){
		u += - s[i] * ( J * s[Pbc(i+1)] + h);
		m += s[i];
	}
	
	//Calculate the current values for the energy and others
	walker[iu] = u;
	walker[ic] = u*u;
	walker[im] = m;
	walker[ix] = m*m;
}

void Reset(int iblk){ //Reset block averages
   
	if(iblk == 1){
		for(int i=0; i<n_props; i++){
			glob_av[i] = 0;
			glob_av2[i] = 0;
		}
	}

	for(int i=0; i<n_props; i++)
		blk_av[i] = 0;

	blk_norm = 0;
	attempted = 0;
	accepted = 0;
}


void Accumulate(void){ //Update block averages

	for(int i=0; i<n_props; i++)	
		blk_av[i] = blk_av[i] + walker[i];
		
	blk_norm = blk_norm + 1.0;
}


void Averages(int iblk){ //Print results for current block
    
	stima[iu] = blk_av[iu] / (blk_norm * (double)nspin); //Energy
	stima[im] = blk_av[im] / (blk_norm * (double)nspin); //Magnetization
	stima[ix] = blk_av[ix] / (blk_norm * (double)nspin) * beta; //Susceptibility
	stima[ic] = (blk_av[ic] - blk_av[iu] * blk_av[iu] / blk_norm) * beta * beta / (blk_norm * (double)nspin); //Specific heat
	
	for(int i=0; i<n_props; i++){
		glob_av[i] += stima[i];
		glob_av2[i] += stima[i]*stima[i];
	}
}


void ConfFinal(void){

	ofstream WriteConf;

	cout << endl << "Print final configuration to file config.final " << endl << endl;
	WriteConf.open("config.final");
	for (int i=0; i<nspin; i++){
	
	WriteConf << s[i] << endl;
	}
	WriteConf.close();

	rnd.SaveSeed();
}

int Pbc(int i){  //Algorithm for periodic boundary conditions

	if(i >= nspin) i = i - nspin;
	else if(i < 0) i = i + nspin;
	return i;
}

void Results(void){

	ofstream energia, calore, magne, susce;
	energia.open("output.energia" + to_string(metro),ios::app);
	calore.open("output.calore" + to_string(metro),ios::app);
	magne.open("output.magne" + to_string(metro),ios::app);
	susce.open("output.susce" + to_string(metro),ios::app);

	//Calculate the errors
	for(int i=0; i < n_props; i++)
		error[i] = sqrt((glob_av2[i] / (double)nblk - pow(glob_av[i] / (double)nblk, 2)) / (double)(nblk - 1));
	
	//Print the results
	energia << glob_av[iu] / (double)nblk << " " << error[iu] << endl;
	calore << glob_av[ic] / (double)nblk << " " << error[ic] << endl;
	magne << glob_av[im] / (double)nblk << " " << error[im] << endl;
	susce << glob_av[ix] / (double)nblk << " " << error[ix] << endl;
	
	energia.close();
	calore.close();
	magne.close();
	susce.close();

}

void Equilibration(void){

	double u;
	ofstream equi;
	equi.open("output.equi");
	
	//Change configuration and only measure the energy
	for (int i = 0; i < nstep; i++){
		Move(metro);
		u = 0;
		for (int j = 0; j < nspin; j++)
			u += - s[j] * ( J * s[Pbc(j+1)] + h);
		equi << i + 1 << " " << u << endl;	
	}
	
	equi.close();
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
