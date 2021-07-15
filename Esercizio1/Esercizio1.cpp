#include <fstream>
#include <iostream>
#include <string>

#include "medieblocchi.h"

using namespace std;

int main(){

	////////////////////////////////////////
	// Calcolo di valore medio e varianza //
	////////////////////////////////////////

	cout << endl << "Inizio calcolo del valore medio e della varianza" << endl << endl;
	
	int n_blocchi = 100;
	int length = 1e4;
	MediaUnif media(length);
	VarUnif varianza(length);
	Random rand = generatore();
	
	for(int i = 0; i< n_blocchi; i++){
	
		// Resetta le variabili per il prossimo blocco
		media.Reset();	
		varianza.Reset();
		
		// Calcola il valore di un blocco
		media.CalcolaBlocco(rand);
		varianza.CalcolaBlocco(rand);
		
		// Aggiungo somme e somme quadratiche
		media.GetSomme();
		varianza.GetSomme();
		
		// Calcola l'errore e stampa su file
		media.Results("valore_medio.txt",i+1);
		varianza.Results("varianza.txt",i+1);
	}	
	
	
	
	////////////////////////////
	// Calcolo del chi quadro //
	////////////////////////////
	
	cout << "Inizio calcolo del chi quadro" << endl << endl;
	
	n_blocchi = 1000; 
	length = 1e4;
	int M = 100;	// Numero sottointervalli di [0,1]
	ChiQuadro chi(length, M);
	
	for (int i = 0; i < n_blocchi; i++){
		chi.Reset();
		chi.CalcolaBlocco(rand);
		chi.PrintBlocco("chi_quadro.txt");
	}
	
	
	
	/////////////////////////////////////
	// Istrogrammi delle distribuzioni //
	/////////////////////////////////////
	
	cout << "Inizio calcolo per istrogrammi delle distribuzioni" << endl << endl;
	
	n_blocchi = 1e4;
	int N[4] = {1,2,10,100};	// Numero elementi su cui fare la media
	
	for (int i = 0; i < 4; i++){
		MediaUnif unif(N[i]);
		Esponenziale expo(N[i]);
		Lorentziana cauchy(N[i]);
		
		for (int j = 0; j < n_blocchi; j++){
			unif.Reset();
			expo.Reset();
			cauchy.Reset();
		
			unif.CalcolaBlocco(rand);
			expo.CalcolaBlocco(rand);
			cauchy.CalcolaBlocco(rand);
		
			unif.PrintBlocco("uniforme_" + to_string(i) + ".txt");
			expo.PrintBlocco("esponenziale_" + to_string(i) + ".txt");
			cauchy.PrintBlocco("cauchy_" + to_string(i) + ".txt");
		}	
	}		
	
	
	
	///////////////////////////
	// Esperimento di Buffon //
	///////////////////////////
	
	cout << "Inizio calcolo per esperimento di Buffon" << endl << endl;
	
	double lunghezza = 1.;
	double distanza = 2*lunghezza;
	length = 1e6;
	n_blocchi = 100;
	Buffon esperimento(length, lunghezza, distanza);
	
	for(int i = 0; i<n_blocchi; i++){
		esperimento.Reset();	
		esperimento.CalcolaBlocco(rand);
		esperimento.GetSomme();
		esperimento.Results("pi_greco.txt",i+1);
	}
	
	return 0;
}
