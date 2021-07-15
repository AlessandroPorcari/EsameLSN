#include "cromosomi.h"

Cromosoma::Cromosoma(int n_città, Random &rand){

	// Genera un vettore contenente numeri da 1 fino a n_città
	for (int i = 0; i < n_città; i++)
		m_geni.push_back(i + 1);
	
	// Rimescola a caso tutti gli elementi del vettore tranne il primo
	for (int i = 1; i < n_città; i++){
		int m = int(rand.Rannyu(1,n_città));
		Swap(i,m);
	}
	
}

Cromosoma::Cromosoma(vector <int> geni){

	m_geni = geni;

}

double Cromosoma::GetCost() const{

	return m_cost;

}

vector <int> Cromosoma::GetGeni() const{

	return m_geni;
}

void Cromosoma::Swap(int index1, int index2){

	int k = m_geni[index1];
	m_geni[index1] = m_geni[index2];
	m_geni[index2] = k;

}

void Cromosoma::Cost(const vector <vector <double>> & città){

	// Come distanza uso L1 perchè così uso la normale nozione di distanza
	m_cost = 0;
	double dist;
	
	for (unsigned i = 0; i < m_geni.size() - 1; i++){
		dist = 0;
		for (unsigned j = 0; j < città[0].size(); j++)
			dist += pow(città[m_geni[i]-1][j] - città[m_geni[i+1]-1][j], 2);
		m_cost += sqrt(dist);
	}
	// Collega l'ultimo punto con il primo
	dist = 0;
	for (unsigned i = 0; i < città[0].size(); i++)
		dist += pow(città[m_geni[m_geni.size()-1]-1][i] - città[0][i], 2);
	m_cost += sqrt(dist);

}

int Cromosoma::PBC(int index){

	// Queste periodic boundary conditions sono pensate per lavorare con le mutazioni, quindi
	// Quando vengono applicate non tengono conto del primo elemento di un cromosoma che deve
	//sempre rimanere uguale a 1

	if (index < int(m_geni.size()))
		return index;
	else
		return index - m_geni.size() + 1;
}

void Cromosoma::Mutazioni(Random &rand){

	// Applica soltanto una mutazione
	double prob = rand.Rannyu();
	
	if (prob < 0.25)
		Mutazione1(rand);
	else if (0.25 <= prob and prob < 0.5)
		Mutazione2(rand);
	else if (0.5 <= prob and prob < 0.75)
		Mutazione3(rand);
	else
		Mutazione4(rand);
		
}

void Cromosoma::Mutazione1(Random &rand){

	// Questa mutazione prende due geni (tranne il primo) e li scambia
		
	int index1;
	int index2;
	
	index1 = int(rand.Rannyu(1,m_geni.size()));
	do{ index2 = int(rand.Rannyu(1,m_geni.size()));
	} while (index1 == index2);
	
	Swap(index1, index2);

}

void Cromosoma::Mutazione2(Random &rand){

	// Questa mutazione prende una serie di geni contigui e li sposta avanti di una certa
	// quantità, i geni che vengono incontrati nel mentre che si sposta vengono spostati
	// indietro invece che avanti

	int primo_blocco = int(rand.Rannyu(1,m_geni.size()));
	int n_contiguous = int(rand.Rannyu(1,m_geni.size()));
	int n_spostamenti;
	vector <int> newGeni;
	
	// Se devo spostare tutti i geni di una quantità
	if (n_contiguous == int(m_geni.size()) - 1){
		n_spostamenti = int(rand.Rannyu(1,m_geni.size() - 1));
		newGeni = m_geni;
		for (unsigned i = 1; i < m_geni.size(); i++)
			m_geni[PBC(i + n_spostamenti)] = newGeni[i];
	}	
	// Se non devo spostare tutti i geni
	else{
		n_spostamenti = int(rand.Rannyu(1,m_geni.size() - n_contiguous));
		for (int i = 0; i < n_contiguous + n_spostamenti; i++)
			newGeni.push_back(m_geni[PBC(primo_blocco + i)]);
		for (int i = 0; i < n_spostamenti; i++)
			m_geni[PBC(primo_blocco + i)] = newGeni[n_contiguous + i];
		for (int i = 0; i < n_contiguous; i++)
			m_geni[PBC(primo_blocco + n_spostamenti + i)] = newGeni[i];
	}

}

void Cromosoma::Mutazione3(Random &rand){

	// Questa mutazione scambia due blocchi della stessa dimensione
	
	int n_contiguous = int(rand.Rannyu(1,m_geni.size()/2)); // blocchi da spostare	
	int primo_blocco = int(rand.Rannyu(1,m_geni.size())); // primo elemento primo blocco
	
	int secondo_blocco; // primo elemento secondo blocco
	secondo_blocco = int(rand.Rannyu(primo_blocco + n_contiguous, m_geni.size() + primo_blocco - n_contiguous));
	
	for (int i = 0; i < n_contiguous; i++)
			Swap(PBC(primo_blocco + i), PBC(secondo_blocco + i));
			
}

void Cromosoma::Mutazione4(Random &rand){

	// Questa mutazione prende una serie di geni e li inverte
	
	int n_contiguous = int(rand.Rannyu(2,m_geni.size()));
	int primo_blocco = int(rand.Rannyu(1,m_geni.size()));
	
	// Inverto i geni scelti
	for (int i = 0; i < n_contiguous/2; i++)
			Swap(PBC(primo_blocco + i), PBC(primo_blocco + n_contiguous - 1 - i));

}

void Cromosoma::Configurazione(const vector <vector<double>> &città){

	ofstream output;
	output.open("posizioni_città.txt");
	for (unsigned i = 0; i < m_geni.size(); i++){
		for (unsigned j = 0; j < città[0].size(); j++)
			output << città[m_geni[i]-1][j] << " "; 
		output << endl;
	}
	
	// Stampa la posizione iniziale per averla anche come finale
	for (unsigned i = 0; i < città[0].size(); i++)
			output << città[m_geni[0]-1][i] << " ";
	output.close();

}

void Cromosoma::PrintCost(int passi){

	// Stampa i risultati
	ofstream output;
	output.open("costo.txt", ios::app);
	output << passi << " " << m_cost << endl;
	output.close();
}

void Cromosoma::Metropolis(double temp, const vector <vector <double>> &pos_città, Random &rand){

	// Propongo una mossa	
	Cromosoma new_geni(m_geni);
	new_geni.Mutazioni(rand);
	new_geni.Cost(pos_città);

	// Accetto la mossa secondo metropolis
	double prob = exp((m_cost - new_geni.GetCost())/temp);
	if (rand.Rannyu() < prob){
		m_geni = new_geni.GetGeni();
		m_cost = new_geni.GetCost();
	}
}

Random generatore(){ 	// Prepare the random number generator

	Random rnd;
	int seed[4];
	int p1, p2;
	
	ifstream Primes("../../Random/Primes");	
	if (Primes.is_open())
		Primes >> p1 >> p2;
	else
		cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("../../Random/seed.in");
	string property;
	if (input.is_open()){
		while (!input.eof()){
			input >> property;
			if(property == "RANDOMSEED"){
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed,p1,p2);
			}
		}
		input.close();
	}
	else
		cerr << "PROBLEM: Unable to open seed.in" << endl;
   
	return rnd;
}

vector <vector <double>> GeneraPosizioni(int n_città, bool circonferenza, Random &rand){

	vector <vector <double>> posizioni;
	vector <double> pos;
	
	for (int i = 0; i < n_città; i++){
		pos.clear();
		
		// Città sulla circonferenza
		if(circonferenza){
			double theta = rand.Rannyu(0,2*M_PI);
			pos.push_back(cos(theta));
			pos.push_back(sin(theta));
		}
		// Città nel quadrato
		else{
			pos.push_back(rand.Rannyu());
			pos.push_back(rand.Rannyu());
		}
		posizioni.push_back(pos);
	}
	
	return posizioni;
}
