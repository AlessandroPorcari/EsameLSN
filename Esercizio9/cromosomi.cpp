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
	// sempre rimanere uguale a 1

	if (index < int(m_geni.size()))
		return index;
	else
		return index - m_geni.size() + 1;
}

void Cromosoma::Mutazioni(double prob, Random &rand){

	// Prova ad applicare le quattro possibili mutazioni
	if (prob > rand.Rannyu())
		Mutazione1(rand);
	if (prob > rand.Rannyu())
		Mutazione2(rand);
	if (prob > rand.Rannyu())
		Mutazione3(rand);
	if (prob > rand.Rannyu())
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

Random generatore(){ 	// Prepare the random number generator

	Random rnd;
	int seed[4];
	int p1, p2;
	
	ifstream Primes("../Random/Primes");	
	if (Primes.is_open())
		Primes >> p1 >> p2;
	else
		cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("../Random/seed.in");
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

void Ordina(vector <Cromosoma> &cromo){

	// Ordina i cromosomi in base al costo, mette prima quelli con il costo più basso
	
	double costo_min, costo;
	unsigned index;
	
	for (unsigned i = 0; i < cromo.size() - 1; i++){
		index = i;
		costo_min = cromo[i].GetCost();
		for (unsigned j = i + 1; j < cromo.size(); j++){
			costo = cromo[j].GetCost();
			if (costo < costo_min){
				costo_min = costo;
				index = j;
			}
		}
		if (index != i){
			Cromosoma ordinato = cromo[i];
			cromo[i] = cromo[index];
			cromo[index] = ordinato;
		}			
	}
}

void PrintCost(const vector <Cromosoma> &cromosomi, int ngen){

	// Calcola la media sulla metà migliore della popolazione
	double media = 0;
	for (unsigned i = 0; i < cromosomi.size()/2; i++)
		media += cromosomi[i].GetCost();
	media = media * 2 / cromosomi.size();

	// Stampa i risultati
	ofstream output;
	output.open("generazioni.txt", ios::app);
	output << ngen << " " << cromosomi[0].GetCost() << " " << media << endl;
	output.close();
}

vector <double> GetProbs(int n_popolazione){

	// Preparo un vettore con le probabilità per la scelta di un certo elemento
	// La probabilità è proporzionale a 1/n
	
	double norma = 0;
	vector <double> probabilities;
	
	for (int i = 0; i < n_popolazione; i++){
		norma += 1/double(i+1);	
		probabilities.push_back(norma);
	}
	for (int i = 0; i < n_popolazione; i++)
		probabilities[i] /= norma;

	return probabilities;

}

vector <int> Selezione(const vector <double> &probs, Random &rand){
	
	vector <int> indici;
	int index1, index2;
	double prob;
	
	// Selezione del primo indice
	index1 = 0;
	prob = rand.Rannyu();
	while (prob > probs[index1]) index1++;
	indici.push_back(index1);
	
	// Selezione del secondo indice
	do {
		index2 = 0;
		prob = rand.Rannyu();
		while (prob > probs[index2]) index2++;
	} while (index1 == index2);	
	indici.push_back(index2);
	
	return indici;

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

vector <Cromosoma> GeneraCromosomi(int popolazione, const vector <vector <double>> &città, Random &rand){

	vector <Cromosoma> cromosomi;
	
	// Genero i cromosomi iniziali e calcolo il costo
	for (int i = 0; i < popolazione; i++){
		cromosomi.push_back(Cromosoma(città.size(), rand));
		cromosomi[i].Cost(città);
	}	
	
	return cromosomi;
}

vector <Cromosoma> Crossover(const Cromosoma &padre, const Cromosoma &madre, double prob, Random &rand){

	vector <Cromosoma> figli;

	// Se faccio il crossover
	if (rand.Rannyu() < prob){
		// Prendo i cromosomi del padre e della madre e preparo i figli
		vector <int> cromoPadre (padre.GetGeni());
		vector <int> cromoMadre (madre.GetGeni());
		vector <int> son1, son2;

		// cut è il primo indice dell'array che verrà tagliato, deve essere almeno 2
		// altrimenti scambio soltanto madre con padre, e non deve essere l'ultimo
		// altrimenti non cambio le sequenze
		int cut = int(rand.Rannyu(2,cromoPadre.size()));
		
		// Copio la parte da tenere
		for (int i = 0; i < cut; i++){
			son1.push_back(cromoPadre[i]);
			son2.push_back(cromoMadre[i]);	
		}

		// Completo il primo figlio
		for (unsigned i = 1; i < cromoMadre.size(); i++){
			bool presente = false;
			// Controlla se il numero è già presente nel figlio
			for (int j = 1; j < cut; j++){
				if (cromoMadre[i] == son1[j])
					presente = true;		
			}
			// Se non è presente inserisce il numero
			if (! presente)
				son1.push_back(cromoMadre[i]);	
		}
		
		// Completo il secondo figlio
		for (unsigned i = 1; i < cromoPadre.size(); i++){
			bool presente = false;
			for (int j = 1; j < cut; j++){
				if (cromoPadre[i] == son2[j])
					presente = true;		
			}
			if (! presente)
				son2.push_back(cromoPadre[i]);	
		}
		
		figli.push_back(Cromosoma(son1));
		figli.push_back(Cromosoma(son2));
		
		return figli;
	}
	else{
		figli.push_back(padre);
		figli.push_back(madre);
		
		return figli;
	}
}
