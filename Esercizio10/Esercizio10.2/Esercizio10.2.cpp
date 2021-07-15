#include <iostream>
#include <fstream>
#include <string>

#include "mpi.h"
#include "cromosomi.h"

using namespace std;

int main(int argc, char *argv[]){

	////////////////////////////////////////////////////////////////
	// Commesso viaggiatore con algoritmo genetico e librerie MPI //
	////////////////////////////////////////////////////////////////

	// Cose per MPI
	int size, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// Controlla di avere 4 processi
	if (size != 4){
		cout << "Servono 4 processi, termino il programma" << endl;
		return -1;
	}

	// Dati iniziali
	Random rand = generatore(rank);
	int n_città = 32;
	int n_popolazione = 1000; // Nota deve essere pari
	int n_generazioni = 200;
	bool circonferenza = false; // Se true genera città sulla circonferenza
	int scambio = 30; // Ogni quante generazioni devo scambiare i dati
	int quanti = 1; // Numero di cromosomi da inviare ogni scambio
	
	// Probabilità di utilizzare gli operatori
	double prob_crossover = 0.7;
	double prob_mutazione = 0.1;
	vector <double> probs (GetProbs(n_popolazione)); // Per la selezione
	
	// Vettori per posizioni delle città
	vector <vector <double>> pos_città(n_città, vector <double> (3));
	if (rank == 0)
		pos_città = GeneraPosizioni(n_città, circonferenza, rand);
	for (int i = 0; i < n_città; i++)
		MPI_Bcast(&pos_città[i][0], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	// Vettori per i cromosomi
	vector <Cromosoma> cromo (GeneraCromosomi(n_popolazione, pos_città, rand));
	vector <Cromosoma> new_gen;
	
	// Ordino in base al costo e salvo il valore della migliore
	Ordina(cromo);
	PrintCost(cromo, 0, rank);
	
	// Ciclo su ogni generazione
	for (int i = 0; i < n_generazioni; i++){			
		// Ciclo per ottenere una nuova generazione
		for (int j = 0; j < n_popolazione/2; j++){	

			// Seleziono due individui e provo a fare il crossover
			vector <int> index (Selezione(probs, rand));	
			vector <Cromosoma> figli (Crossover(cromo[index[0]], cromo[index[1]], prob_crossover, rand));	
				
			// Provo ad aggiungere mutazioni
			figli[0].Mutazioni(prob_mutazione, rand);
			figli[1].Mutazioni(prob_mutazione, rand);	
				
			// Calcolo il costo
			figli[0].Cost(pos_città);
			figli[1].Cost(pos_città);	
					
			// Aggiungo alla nuova generazione
			new_gen.push_back(figli[0]);
			new_gen.push_back(figli[1]);
		}		
		
		// Copio la nuova generazione sulla vecchia e svuoto la nuova
		cromo = new_gen;
		new_gen.clear();
		
		// Ordino in base al costo e salvo il valore della migliore
		Ordina(cromo);
		PrintCost(cromo, i + 1, rank);
		
		// Se sono passate un tot di generazioni faccio lo scambio
		if ((i+1) % scambio == 0)
			Scambia(quanti, rank, cromo, rand);
		
		if (((i+1) % 10 == 0) and rank == 0)
			cout << "Finito la generazione " << i+1 << "/" << n_generazioni << endl; 
	}	
	
	// Stampo su file l'ordine finale con le posizioni
	cromo[0].Configurazione(pos_città, rank);
	
	MPI_Finalize();
	
	return 0;
}
