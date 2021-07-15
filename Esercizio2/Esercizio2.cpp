#include <cmath>
#include <fstream>
#include <iostream>

#include "medieblocchi.h"

using namespace std;

int main(){
	
	////////////////////////////
	// Calcolo dell'integrale //
	////////////////////////////	
	
	cout << endl << "Inizio calcolo dell'integrale con le due distribuzioni" << endl;
	
	int length = 1e5;
	int n_blocchi = 100;
	Uniforme uniform(length);
	Lineare linear(length);
	Random rand = generatore();
	
	for (int i = 0; i < n_blocchi; i++){
		// Resetta le variabili per il prossimo blocco
		uniform.Reset();
		linear.Reset();
		
		// Calcola il valore di un blocco
		uniform.CalcolaBlocco(rand);
		linear.CalcolaBlocco(rand);
		
		// Aggiungi somme e somme quadratiche
		uniform.GetSomme();
		linear.GetSomme();
		
		// Calcula l'errore e stampa su file
		uniform.Results("distribuzione_uniforme.txt", i+1);
		linear.Results("distribuzione_lineare.txt", i+1);
	}
	
	
	
	////////////////////////////////////////////////////
	// Calcolo della distanza quadratica media con RW //
	////////////////////////////////////////////////////
	
	cout << endl << "Inizio simulazioni dei random walk" << endl << endl;
	
	length = 1e4;
	n_blocchi = 100;
	int size = 100; //	Lunghezza di un random walk
	double step = 2.;
	PassoReticolo reticolo(length, size, step);
	PassoIsotropo isotropo(length, size, step);
	
	for (int i = 0; i < n_blocchi; i++){
		reticolo.Reset();
		isotropo.Reset();
		reticolo.CalcolaBlocco(rand);
		isotropo.CalcolaBlocco(rand);
		reticolo.GetSomme();
		isotropo.GetSomme();
		if ((i+1)%5 == 0)
			cout << "Finito blocco " << i+1 << " su " << n_blocchi << endl;	
	}
	
	reticolo.Results("reticolo.txt", n_blocchi);
	isotropo.Results("isotropo.txt", n_blocchi);
	
	return 0;
}
