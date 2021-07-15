#include <iostream>

#include "medieblocchi.h"

using namespace std;

int main(){

	//////////////////////////////////////////
	// Calcolo dell'energia dato mu e sigma //
	//////////////////////////////////////////

	cout << endl << "Inizio calcolo dell'energia" << endl << endl;

	double posizione = 1;	// Posizione iniziale
	int length = 100000;
	int n_blocchi = 100;
	Random rand = generatore();
	bool griglia = false;	// True calcola il valore dell'energia su una griglia
	
	if (!griglia){
		double mu = 0.80;
		double sigma = 0.62;
		double step = 2.5;
		Metropolis metro(length, step, posizione, mu, sigma);
		
		for(int i = 0; i < n_blocchi; i++){
			metro.Reset();
			metro.CalcolaBlocco(rand);
			metro.GetSomme();
			metro.Results("energia.txt", i+1);
			
			if ((i+1) % 10 == 0)
				cout << "Finito blocco " << i+1 << "/" << n_blocchi << endl;
		}
	}
	else {
		double griglia [] = {0.7, 0.55, 0.9, 0.7}; // Estremi della griglia
		int punti [] = {21, 16}; // Numero righe e colonne
		CalcolaGriglia(griglia, punti, length, posizione, rand);
	}
	
	return 0;
}
