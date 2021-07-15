#include <iostream>

#include "medieblocchi.h"

using namespace std;

int main(){

	/////////////////////////////////////////////////////////
	// Calcolo del valore medio con distribuzione uniforme //
	/////////////////////////////////////////////////////////

	cout << endl << "Inizio calcolo del valore medio con distribuzione uniforme" << endl;

	double posizione[3] = {100., 100., 100.};	//posizione iniziale
	int length = 1e4;	
	int n_blocchi = 100;	
	double step1 = 1.2;
	double step2 = 3;
	Uniforme100 unif1(length, step1, posizione);
	Uniforme210 unif2(length, step2, posizione);
	Random rand = generatore();
	
	for(int i = 0; i < n_blocchi; i++){
		// Resetta le variabili per il prossimo blocco
		unif1.Reset();
		unif2.Reset();
		
		// Calcola il valore del blocco
		unif1.CalcolaBlocco(rand);
		unif2.CalcolaBlocco(rand);
		
		// Aggiungi le somme e le somme quadratiche
		unif1.GetSomme();
		unif2.GetSomme();
		
		// Calcola l'errore e stampa su file
		unif1.Results("uniforme_100.txt", i+1);
		unif2.Results("uniforme_210.txt", i+1);
	}



	//////////////////////////////////////////////////////////
	// Calcolo del valore medio con distribuzione gaussiana //
	//////////////////////////////////////////////////////////

	cout << endl << "Inizio calcolo del valore medio con distribuzione gaussiana" << endl << endl;
	
	step1 = 0.8;
	step2 = 1.8;
	Gauss100 gauss1(length, step1, posizione);
	Gauss210 gauss2(length, step2, posizione);
		
	for(int i = 0; i < n_blocchi; i++){
		gauss1.Reset();
		gauss2.Reset();
		
		gauss1.CalcolaBlocco(rand);
		gauss2.CalcolaBlocco(rand);
		
		gauss1.GetSomme();
		gauss2.GetSomme();
		
		gauss1.Results("gaussiana_100.txt", i+1);
		gauss2.Results("gaussiana_210.txt", i+1);
	}
	
	return 0;
}
