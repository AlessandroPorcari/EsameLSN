#ifndef __cromosomi_h__
#define __cromosomi_h__

#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>
#include "../../Random/random.h"

using namespace std;

class Cromosoma{

	public:
		
		// Costruttori
		Cromosoma(int n_città, Random &rand);
		Cromosoma(vector <int> geni);
		
		double GetCost() const;
		vector <int> GetGeni() const;
		void PrintCost(int passi);
		void Mutazioni(Random &rand);
		void Cost(const vector <vector <double>> & città);
		void Configurazione(const vector <vector<double>> &città);
		void Metropolis(double temp, const vector <vector <double>> &pos_città, Random &rand);
		
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

Random generatore();
vector <vector <double>> GeneraPosizioni(int n_città, bool circonferenza, Random &rand);

#endif
