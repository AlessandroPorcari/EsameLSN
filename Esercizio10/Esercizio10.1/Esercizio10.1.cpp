#include <iostream>
#include <fstream>
#include <string>

#include "cromosomi.h"

using namespace std;

int main(){
	
	//////////////////////////////////////////////////
	// Commesso viaggiatore con simulated annealing //
	//////////////////////////////////////////////////
	
	// Dati iniziali
	Random rand = generatore();
	int n_città = 32;
	int n_step = 100000; 	// Numero di step metropolis ogni temperatura
	bool circonferenza = true; // Se true genera città sulla circonferenza
	double temp_in = 5;
	double temp_fin = 0.01;
	int n_temp = 20;
	
	// Vettori per posizioni, cromosoma iniziale e cromosoma migliore
	vector <vector <double>> pos_città (GeneraPosizioni(n_città, circonferenza, rand));
	Cromosoma cromo(n_città, rand);
	cromo.Cost(pos_città);
	Cromosoma migliore(cromo);
	
	// Salvo il valore del cromosoma iniziale
	migliore.PrintCost(0);
	
	cout << "Comincio a utilizzare Metropolis" << endl;
	
	// Dopo ogni temperatura divido questa per ks
	double k = pow(temp_in / temp_fin, 1 / double(n_temp - 1));
	double temp = temp_in;
	
	// Ciclo sulle diverse temperature
	for (int i = 0; i < n_temp; i++){	
		for (int j = 0; j < n_step; j++){	
		
			// Prova a fare un passo Metropolis		
			cromo.Metropolis(temp, pos_città, rand);		
			if (cromo.GetCost() < migliore.GetCost())
				migliore = cromo;
			if ((j+1)%100 == 0)
				migliore.PrintCost(i * n_step + j + 1);
		}	
		
		// Diminuisco la temperatura
		temp /= k;
		cout << "Finito la temperatura " << i + 1 << "/" << n_temp << endl; 
	}
	
	// Stampo su file l'ordine finale con le posizioni
	migliore.Configurazione(pos_città);
	
	return 0;
}
