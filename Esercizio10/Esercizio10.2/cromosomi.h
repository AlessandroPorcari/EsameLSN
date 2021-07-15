#ifndef __cromosomi_h__
#define __cromosomi_h__

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include "mpi.h"
#include "../../Random/random.h"

using namespace std;

class Cromosoma{

	public:
		
		// Costruttori
		Cromosoma(int n_città, Random &rand);
		Cromosoma(vector <int> geni);
		
		double GetCost() const;
		vector <int> GetGeni() const;
		void Cost(const vector <vector <double>> & città);
		void Mutazioni(double prob, Random &rand);
		void Configurazione(const vector <vector<double>> &città, int rank);
		
	private:
	
		vector <int> m_geni;
		double m_cost;
		
		int PBC(int index);
		void Swap(int index1, int index2);
		void Mutazione1(Random &rand);
		void Mutazione2(Random &rand);
		void Mutazione3(Random &rand);
		void Mutazione4(Random &rand);
};

// Funzione per il generatore casuale
Random generatore(int rank);

// Funzioni per operare sui cromosomi
void Ordina(vector <Cromosoma>  &cromo);
void PrintCost(const vector <Cromosoma> &cromosomi, int ngen, int rank);
vector <double> GetProbs(int n_popolazione);
vector <int> Selezione(const vector <double> &probs, Random &rand);
vector <vector <double>> GeneraPosizioni(int n_città, bool circonferenza, Random &rand);
vector <Cromosoma> GeneraCromosomi(int popolazione, const vector <vector <double>> &città, Random &rand);
vector <Cromosoma> Crossover(const Cromosoma &padre, const Cromosoma &madre, double prob, Random &rand);

// Funzione per scambio informazioni
void Scambia(int quanti, int rank, vector <Cromosoma> &cromo, Random &rand);

#endif
