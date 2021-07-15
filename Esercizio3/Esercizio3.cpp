#include <iostream>
#include <fstream>

#include "medieblocchi.h"

using namespace std;

int main(){
	
	// Dati iniziali
	double S0 = 100.;
	double K = 100.;
	double T = 1.;
	double r = 0.1;
	double sigma = 0.25;
	int steps = 100;
	
	// Contenitori dati e altro
	Random rand = generatore();	
	int length = 1e4;
	int n_blocchi = 100;
	CallDiretto callfin(length, r, sigma, T, S0, K);
	PutDiretto putfin(length, r, sigma, T, S0, K);
	CallDiscreto calldis(length, r, sigma, T, S0, K, steps);
	PutDiscreto putdis(length, r, sigma, T, S0, K, steps);
	
	for(int i = 0; i < n_blocchi; i++){
		// Resetta le variabili per il prossimo blocco
		callfin.Reset();
		putfin.Reset();
		calldis.Reset();
		putdis.Reset();
		
		// Calcola il valore di un blocco
		callfin.CalcolaBlocco(rand);
		putfin.CalcolaBlocco(rand);
		calldis.CalcolaBlocco(rand);
		putdis.CalcolaBlocco(rand);

		// Aggiungi somme e somme quadratiche
		callfin.GetSomme();
		putfin.GetSomme();
		calldis.GetSomme();
		putdis.GetSomme();

		// Calcola l'errore e stampa i risultati
		callfin.Results("call_diretto.txt", i+1);
		putfin.Results("put_diretto.txt", i+1);
		calldis.Results("call_discreto.txt", i+1);
		putdis.Results("put_discreto.txt", i+1);
		
		if ((i+1)%5 == 0)
			cout << "Finito blocco " << i+1 << " su " << n_blocchi << endl;	
	}

	return 0;
}
