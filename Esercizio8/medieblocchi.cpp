#include "medieblocchi.h"

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

MedieBlocchi::MedieBlocchi(int length){

	m_length = length;
	m_somma = 0;
	m_somma2 = 0;
	m_blocco = 0;
	
}

void MedieBlocchi::GetSomme(){

	m_somma += m_blocco;
	m_somma2 += m_blocco * m_blocco;
	
}

void MedieBlocchi::Results(string nome, int num){

	// Calcola l'errore
	if (num == 1) m_errore = 0;
	else m_errore = sqrt((m_somma2/(double)num-pow(m_somma/(double)num,2))/((double)num-1));

	// Stampa a video
	ofstream output;
	output.open(nome,ios::app);
	output << num * m_length << " " << m_somma/(double)num << " " << m_errore << endl;	
	output.close();
}

Metropolis::Metropolis(int length, double step, double position, double mu, double sigma) : MedieBlocchi(length){

	m_step = step;
	m_mu = mu;
	m_sigma = sigma;
	m_position = position;
	m_prob = GetProbability(m_position);
	// Salvo la probabilità del punto in cui sono per non stare a ricalcolarla se non mi sposto

}

void Metropolis::CalcolaBlocco(Random &rand){

	double new_pos;
	double prob;
	int accepted = 0;
	
	for (int i = 0; i < m_length; i++){
		// Calcola una nuova posizione
		new_pos = m_position + m_step * rand.Rannyu(-1, 1);
		
		// Usa l'algoritmo di Metropolis
		prob = GetProbability(new_pos);
		if ((prob / m_prob) > rand.Rannyu()){
			accepted ++;
			m_position = new_pos;
			m_prob = prob;
		} 
		m_blocco += Integranda(m_position);
		
		// Togliere commmento per stampare su file la posizione ogni 10 passi
		//if ((i+1)%10 == 0)
		//	PrintPosizione("posizioni.txt");
	}
		
	m_blocco /= (double)m_length;
	// Togliere il commento se voglio controllare accettazione
	//cout << "Acceptance rate = " << (double)accepted / (double)m_length << endl;
}

void Metropolis::PrintPosizione(string nome){

	// Stampa su file la posizione corrente
	ofstream output;
	output.open(nome,ios::app);
	output << m_position << endl;	
	output.close();

}

double Metropolis::GetProbability(double pos){

	// Calcola la probabilità data dal modulo quadro
	double prob = exp(-pow(pos / m_sigma, 2));
	return prob * (1. + cosh(2. * pos * m_mu / pow(m_sigma, 2)));

}

double Metropolis::Integranda(double pos){

	double V;
	double T;

	// Calcola il termine (H psi) / psi	
	V = pow(pos, 4) - 2.5 * pow(pos, 2);
	T = pow(m_sigma, 2) - pow(pos + m_mu, 2);
	T += 4. * m_mu * pos / ( 1 + exp(-2 * m_mu * pos / pow(m_sigma, 2)));
	T /= 2. * pow(m_sigma, 4);
	
	return T + V;

}

void CalcolaGriglia(double griglia [4], int punti [2], int length, double posizione, Random & rand){

	// Calcola il valore dell'energia su una griglia facendo length passi per ogni punto
	double mu;
	double sigma;
	
	ofstream output;
	output.open("griglia.txt");	

	for (int i = 0; i < punti[0]; i++) {
		for (int j = 0; j < punti[1]; j++) {
			// Prendi un punto su una griglia
			mu = griglia[0] + (griglia[2] - griglia[0]) * double(i) / double(punti[0] - 1);
			sigma = griglia[1] + (griglia[3] - griglia[1]) * double(j) / double(punti[1] - 1);
			Metropolis metro(length, 0.5 + sigma + mu, posizione, mu, sigma);
			
			// Calcola il valore dell'energia e stampa su file
			metro.CalcolaBlocco(rand);
			output << mu << " " << sigma << " " << metro.Energia() << endl;
		}
		cout << "Finito " << i+1 << "/" << punti[0] << endl;
	}
	
	output.close();
}
