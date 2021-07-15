/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include "MolDyn_NVE.h"

using namespace std;

int main(){ 
	Input();	//Inizialization
	int nconf = 1;
	for (int i = 0; i < nblocks; i++){
		Reset();
		for (int istep=1; istep <= nstep/nblocks; ++istep){
			Move();	//Move particles with Verlet algorithm
			if (istep%10 == 0){
				Measure();	//Properties measurement
				//ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 				
				nconf += 1;
			}
		}
		GetSums();
		Uncertainty(i+1);
		cout << "Finished block " << i+1 << "/" << nblocks << endl;
	}	
	ConfFinal();	//Write final configuration to restart
	ConfPrevious();	//Write previous configuration to restart
	
	return 0;
}


void Input(){	//Prepare all stuff for the simulation
	ifstream ReadInput,ReadConf;
	double ep, ek, pr, et, vir;

	cout << "Classic Lennard-Jones fluid        " << endl;
	cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
	cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
	cout << "The program uses Lennard-Jones units " << endl;

	seed = 1;	//Set seed for random numbers
	srand(seed);	//Initialize random number generator
  
	ReadInput.open("input.dat");	//Read input

	ReadInput >> temp;

	ReadInput >> npart;
	cout << "Number of particles = " << npart << endl;

	ReadInput >> rho;
	cout << "Density of particles = " << rho << endl;
	vol = (double)npart/rho;
	cout << "Volume of the simulation box = " << vol << endl;
	box = pow(vol,1.0/3.0);
	cout << "Edge of the simulation box = " << box << endl;

	ReadInput >> rcut;
	ReadInput >> delta;
	ReadInput >> nstep;
	ReadInput >> nblocks;

	cout << "The program integrates Newton equations with the Verlet method " << endl;
	cout << "Time step = " << delta << endl;
	cout << "Number of steps = " << nstep << endl << endl;
	cout << "Number of blocks = " << nblocks << endl << endl;
	
	//Initialize the variables for the sum and the sum squared
	sum_pot = 0;
	sum_kin = 0;
	sum_etot= 0;
	sum_temp = 0;
	sum2_pot = 0;
	sum2_kin = 0;
	sum2_etot = 0;
	sum2_temp = 0;

	//Read initial configuration
	cout << "Read initial configuration from file config.0 " << endl << endl;
	ReadConf.open("config.0");
	for (int i=0; i<npart; ++i){
		ReadConf >> x[i] >> y[i] >> z[i];
		x[i] = x[i] * box;
		y[i] = y[i] * box;
		z[i] = z[i] * box;
	}
	ReadConf.close();

	//Read if previous configuration was given
	ReadInput >> previous;
	ReadInput.close();
	
	//Prepare initial velocities if previous configuration was not given
	if (!previous){
		cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
		double sumv[3] = {0.0, 0.0, 0.0};
		for (int i=0; i<npart; ++i){
			vx[i] = rand()/double(RAND_MAX) - 0.5;
			vy[i] = rand()/double(RAND_MAX) - 0.5;
			vz[i] = rand()/double(RAND_MAX) - 0.5;

			sumv[0] += vx[i];
			sumv[1] += vy[i];
			sumv[2] += vz[i];
		}
		for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
		double sumv2 = 0.0, fs;
		for (int i=0; i<npart; ++i){
			vx[i] = vx[i] - sumv[0];
			vy[i] = vy[i] - sumv[1];
			vz[i] = vz[i] - sumv[2];

			sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
		}
		sumv2 /= (double)npart;

		fs = sqrt(3 * temp / sumv2);	// fs = velocity scale factor 
		for (int i=0; i<npart; ++i){
			vx[i] *= fs;
			vy[i] *= fs;
			vz[i] *= fs;

			xold[i] = Pbc(x[i] - vx[i] * delta);
			yold[i] = Pbc(y[i] - vy[i] * delta);
			zold[i] = Pbc(z[i] - vz[i] * delta);
		}
	}
	
	//Read previous configuration if given
	else{
		cout << "Read previous configuration from file config.-1 " << endl << endl;
		ReadConf.open("config.-1");
		for (int i=0; i<npart; ++i){
			ReadConf >> xold[i] >> yold[i] >> zold[i];
			xold[i] = xold[i] * box;
			yold[i] = yold[i] * box;
			zold[i] = zold[i] * box;
		}
		ReadConf.close();
		
		double sumv2 = 0.0, fs;
		for (int i=0; i<npart; ++i){
			//Get v(t - dt/2)
			vx[i] = Pbc(x[i] - xold[i])/delta;
			vy[i] = Pbc(y[i] - yold[i])/delta;
			vz[i] = Pbc(z[i] - zold[i])/delta;

			sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
		}
		sumv2 /= (double)npart;

		fs = sqrt(3 * temp / sumv2);	// fs = velocity scale factor 
		for (int i=0; i<npart; ++i){
			vx[i] *= fs;
			vy[i] *= fs;
			vz[i] *= fs;

			x[i] = Pbc(xold[i] + vx[i] * delta);
			y[i] = Pbc(yold[i] + vy[i] * delta);
			z[i] = Pbc(zold[i] + vz[i] * delta);
		}
	}
	return;
}


void Move(void){	//Move particles with Verlet algorithm
	double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

	for(int i=0; i<npart; ++i){	//Force acting on particle i
		fx[i] = Force(i,0);
		fy[i] = Force(i,1);
		fz[i] = Force(i,2);
	}

	for(int i=0; i<npart; ++i){	//Verlet integration scheme

		xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
		ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
		znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

		vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
		vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
		vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

		xold[i] = x[i];
		yold[i] = y[i];
		zold[i] = z[i];

		x[i] = xnew;
		y[i] = ynew;
		z[i] = znew;
	}
	return;
}

double Force(int ip, int idir){	//Compute forces as -Grad_ip V(r)
	double f=0.0;
	double dvec[3], dr;

	for (int i=0; i<npart; ++i){
		if(i != ip){
			dvec[0] = Pbc( x[ip] - x[i] );	// distance ip-i in pbc
			dvec[1] = Pbc( y[ip] - y[i] );
			dvec[2] = Pbc( z[ip] - z[i] );

			dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
			dr = sqrt(dr);

			if(dr < rcut){
				f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8));	// -Grad_ip V(r)
      		}
		}
	}
  
	return f;
}

void Measure(){	//Properties measurement
	int bin;
	double v, t, vij;
	double dx, dy, dz, dr;
	ofstream Epot, Ekin, Etot, Temp;

	Epot.open("output_epot.dat",ios::app);
	Ekin.open("output_ekin.dat",ios::app);
	Temp.open("output_temp.dat",ios::app);
	Etot.open("output_etot.dat",ios::app);

	v = 0.0;	//reset observables
	t = 0.0;
	
	//reset the hystogram of g(r)
	for (int i=0; i < nbins; i++) gofr[i] = 0;

	//cycle over pairs of particles
	for (int i=0; i<npart-1; ++i){
		for (int j=i+1; j<npart; ++j){

			dx = Pbc( xold[i] - xold[j] );	// here I use old configurations [old = r(t)]
			dy = Pbc( yold[i] - yold[j] );	// to be compatible with EKin which uses v(t)
			dz = Pbc( zold[i] - zold[j] );	// => EPot should be computed with r(t)

			dr = dx*dx + dy*dy + dz*dz;
			dr = sqrt(dr);

			//update of the histogram of g(r)
			if(dr < box/2)
				gofr[int(dr * nbins * 2. / box)] += 2;

			//Potential energy
			if(dr < rcut)
				v += 4.0/pow(dr,12) - 4.0/pow(dr,6);
		}          
	}
	
	//Kinetic energy
	for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);

	// Add to the block calculation
	blocco_pot += v;
	blocco_kin += t;
	blocco_etot += v + t;
	blocco_temp += (2./3.) * t;
	for (int i = 0; i < nbins; i++)
 		blocco_gofr[i] += gofr[i];
   
	stima_pot = v/(double)npart;	//Potential energy per particle
	stima_kin = t/(double)npart;	//Kinetic energy per particle
	stima_temp = (2.0 / 3.0) * t/(double)npart;	//Temperature
	stima_etot = (t+v)/(double)npart;	//Total energy per particle

	Epot << stima_pot  << endl;
	Ekin << stima_kin  << endl;
	Temp << stima_temp << endl;
	Etot << stima_etot << endl;

	Epot.close();
	Ekin.close();
	Temp.close();
	Etot.close();

	return;
}


void ConfFinal(void){	//Write final configuration
	ofstream WriteConf;

	cout << "Print final configuration to file config.final " << endl << endl;
	WriteConf.open("config.final");

	for (int i=0; i<npart; ++i){
		WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
	}
	WriteConf.close();
	return;
}

void ConfPrevious(void){	//Write previous configuration
	ofstream WriteConf;

	cout << "Print previous configuration to file config.previous " << endl << endl;
	WriteConf.open("config.previous");

	for (int i=0; i<npart; ++i){
		WriteConf << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
	}
	WriteConf.close();
	return;
}

void ConfXYZ(int nconf){	//Write configuration in .xyz format
	ofstream WriteXYZ;

	WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
	WriteXYZ << npart << endl;
	WriteXYZ << "This is only a comment!" << endl;
	for (int i=0; i<npart; ++i){
		WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
	}
	WriteXYZ.close();
}

double Pbc(double r){	//Algorithm for periodic boundary conditions with side L=box
	return r - box * rint(r/box);
}


void Reset(void){

	// Reset the variables for a block
	blocco_pot = 0;
	blocco_kin = 0;
	blocco_etot = 0;
	blocco_temp = 0;
	for (int i = 0; i < nbins; i++)
 		blocco_gofr[i] = 0;

}

void GetSums(){

	double deltaV;

	// Calculate the value of the block
	blocco_pot = 10 * blocco_pot * (double)nblocks / ((double)nstep * (double)npart);
	blocco_kin = 10 * blocco_kin * (double)nblocks / ((double)nstep * (double)npart);
	blocco_etot = 10 * blocco_etot * (double)nblocks / ((double)nstep * (double)npart);
	blocco_temp = 10 * blocco_temp * (double)nblocks / ((double)nstep * (double)npart);

	for (int i = 0; i < nbins; i++){
		deltaV = M_PI/24. * pow(box/nbins,3) * (1. + 12.*pow(i+1,2));
 		blocco_gofr[i] = 10 * blocco_gofr[i] * (double)nblocks / ((double)nstep * rho * (double)npart * deltaV);
	}

	// Add to mean and mean squared
	sum_pot += blocco_pot;
	sum_kin += blocco_kin;
	sum_etot += blocco_etot;
	sum_temp += blocco_temp;
	sum2_pot += blocco_pot*blocco_pot;
	sum2_kin += blocco_kin*blocco_kin;
	sum2_etot += blocco_etot*blocco_etot;
	sum2_temp += blocco_temp*blocco_temp;

	for (int i = 0; i < nbins; i++){
 		sum_gofr[i] += blocco_gofr[i];
 		sum2_gofr[i] += blocco_gofr[i] * blocco_gofr[i];
	}
	
}

void Uncertainty(int num_blocks){

	ofstream Epot, Ekin, Etot, Temp, Gave;

	Epot.open("ave_epot.dat",ios::app);
	Ekin.open("ave_ekin.dat",ios::app);
	Temp.open("ave_temp.dat",ios::app);
	Etot.open("ave_etot.dat",ios::app);	

	// Calculate the error
	if (num_blocks == 1){
		error_pot = 0;
		error_kin = 0;
		error_etot = 0;	
		error_temp = 0;
	}
	else {
		error_pot = sqrt((sum2_pot/(double)num_blocks-pow(sum_pot/(double)num_blocks,2))/((double)num_blocks-1));
		error_kin = sqrt((sum2_kin/(double)num_blocks-pow(sum_kin/(double)num_blocks,2))/((double)num_blocks-1));
		error_etot = sqrt((sum2_etot/(double)num_blocks-pow(sum_etot/(double)num_blocks,2))/((double)num_blocks-1));
		error_temp = sqrt((sum2_temp/(double)num_blocks-pow(sum_temp/(double)num_blocks,2))/((double)num_blocks-1));
	}

	Epot << sum_pot/(double)num_blocks << " " << error_pot << endl;
	Ekin << sum_kin/(double)num_blocks << " " << error_kin << endl;
	Temp << sum_temp/(double)num_blocks << " " << error_temp << endl;
	Etot << sum_etot/(double)num_blocks << " " << error_etot << endl;
	
	//g(r)		
	if(num_blocks == nblocks){
		Gave.open("ave.gofr.dat");
		for (int i = 0; i < nbins; i++){
			error_gofr = sqrt((sum2_gofr[i]/(double)num_blocks-pow(sum_gofr[i]/(double)num_blocks,2))/((double)num_blocks-1));
			Gave << sum_gofr[i]/(double)num_blocks << " " << error_gofr << endl;
		}
		Gave.close();
	}

	Epot.close();
	Ekin.close();
	Temp.close();
	Etot.close();	
		
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
