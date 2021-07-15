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
	
}

void MedieBlocchi::GetSomme(){

	m_somma += m_blocco;
	m_somma2 += m_blocco * m_blocco;
	
}

void MedieBlocchi::Results(string nome, int num){

	// Calcola l'errore
	if (num == 1) m_errore = 0;
	else m_errore = sqrt((m_somma2/(double)num-pow(m_somma/(double)num,2))/((double)num-1));

	// Stampa su file
	ofstream output;
	output.open(nome,ios::app);
	output << num * m_length << " " << m_somma/(double)num << " " << m_errore << endl;	
	output.close();

}

CallDiretto::CallDiretto(int length, double mu, double sigma, double T, double S0, double K) : MedieBlocchi(length){

	m_mu = mu;
	m_sigma = sigma;
	m_T = T;
	m_S0 = S0;
	m_K = K;

}

void CallDiretto::CalcolaBlocco(Random &rand){

	// Per il prezzo della Call devo generare una variabile gaussiana e poi
	// controllare se ho guadagnato rispetto a K. Infine devo sconto di exp(-mu*T)
	for (int i = 0; i < m_length; i++)
		m_blocco += max(0., m_S0 * exp((m_mu - m_sigma * m_sigma / 2) * m_T + m_sigma * rand.Gauss(0, sqrt(m_T))) - m_K);
	
	m_blocco = m_blocco * exp(-m_mu * m_T) / (double)m_length;
}

PutDiretto::PutDiretto(int length, double mu, double sigma, double T, double S0, double K) : MedieBlocchi(length){

	m_mu = mu;
	m_sigma = sigma;
	m_T = T;
	m_S0 = S0;
	m_K = K;

}

void PutDiretto::CalcolaBlocco(Random &rand){

	// Per il prezzo della Put devo generare una variabile gaussiana e poi
	// controllare se ho guadagnato rispetto a K. Infine devo sconto di exp(-mu*T)
	for (int i = 0; i < m_length; i++)
		m_blocco += max(0., m_K - m_S0 * exp((m_mu - m_sigma * m_sigma / 2) * m_T + m_sigma * rand.Gauss(0, sqrt(m_T))));
	
	m_blocco = m_blocco * exp(-m_mu * m_T) / (double)m_length;
}

CallDiscreto::CallDiscreto(int length, double mu, double sigma, double T, double S0, double K, int steps) : MedieBlocchi(length){

	m_mu = mu;
	m_sigma = sigma;
	m_T = T;
	m_deltaT = T / (double)steps;
	m_S0 = S0;
	m_K = K;
	m_steps = steps;

}

void CallDiscreto::CalcolaBlocco(Random &rand){
	
	// Per il prezzo della Call ogni passo devo generare una variabile gaussiana
	// Arrivato al tempo finale controllo se ho guadagnato rispetto a K e sconto di exp(-mu*T)
	for (int i = 0; i < m_length; i++){
		double S = m_S0;
		for (int j = 0; j < m_steps; j++){
			S = S * exp((m_mu - m_sigma * m_sigma / 2) * m_deltaT + m_sigma * rand.Gauss(0.,1.) * sqrt(m_deltaT));
		}
		m_blocco += max(0., S - m_K);
	}

	m_blocco = m_blocco * exp(-m_mu * m_T) / (double)m_length;
}

PutDiscreto::PutDiscreto(int length, double mu, double sigma, double T, double S0, double K, int steps) : MedieBlocchi(length){

	m_mu = mu;
	m_sigma = sigma;
	m_T = T;
	m_deltaT = T / (double)steps;
	m_S0 = S0;
	m_K = K;
	m_steps = steps;

}

void PutDiscreto::CalcolaBlocco(Random &rand){

	// Per il prezzo della Put ogni passo devo generare una variabile gaussiana
	// Arrivato al tempo finale controllo se ho guadagnato rispetto a K e sconto di exp(-mu*T)
	for (int i = 0; i < m_length; i++){
		double S = m_S0;
		for (int j = 0; j < m_steps; j++){
			S = S * exp((m_mu - m_sigma * m_sigma / 2) * m_deltaT + m_sigma * rand.Gauss(0.,1.) * sqrt(m_deltaT));
		}
		m_blocco += max(0., m_K - S);
	}

	m_blocco = m_blocco * exp(-m_mu * m_T) / (double)m_length;
}
