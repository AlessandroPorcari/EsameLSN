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

	// Stampa a video
	ofstream output;
	output.open(nome,ios::app);
	output << num * m_length << " " << m_somma/(double)num << " " << m_errore << endl;	
	output.close();

}

void MedieBlocchi::PrintBlocco(string nome){

	// Stampa a video
	ofstream output;
	output.open(nome,ios::app);
	output << m_blocco << endl;	
	output.close();

}

MediaUnif::MediaUnif(int length) : MedieBlocchi(length){

}

void MediaUnif::CalcolaBlocco(Random &rand){

	// Calcola valore medio distribuzione uniforme
	for (int i = 0; i < m_length; i++)
		m_blocco += rand.Rannyu();
		
	m_blocco /= (double)m_length;
}

VarUnif::VarUnif(int length) : MedieBlocchi(length){

}

void VarUnif::CalcolaBlocco(Random &rand){

	// Calcola varianza distribuzione uniforme
	for (int i = 0; i < m_length; i++)
		m_blocco += pow(rand.Rannyu()-0.5,2);
		
	m_blocco /= (double)m_length;
}

ChiQuadro::ChiQuadro(int length, int M) : MedieBlocchi(length){

	m_M = M;
	m_counter = new int[M];

}

ChiQuadro::~ChiQuadro(){

	delete [] m_counter;

}

void ChiQuadro::CalcolaBlocco(Random &rand){

	// Resetta il contatore dei sottointervalli
	for (int i = 0; i < m_M; i++) 
		m_counter[i] = 0;
		
	// Conta il numero di eventi in ogni sottointervallo
	for (int i = 0; i < m_length; i++) 
		m_counter[int(m_M * rand.Rannyu())] ++;
		
	// Calcola il chi quadro
	for (int i = 0; i < m_M; i++) 
		m_blocco += pow(m_counter[i] - (double)m_length/(double)m_M,2);
	
	m_blocco = m_blocco * (double)m_M/(double)m_length;
}

Esponenziale::Esponenziale(int length) : MedieBlocchi(length){

}

void Esponenziale::CalcolaBlocco(Random &rand){

	// Calcola il valore medio di variabili esponenziali
	for (int i = 0; i < m_length; i++)
		m_blocco += rand.Exponential(1.);
		
	m_blocco /= (double)m_length;
}

Lorentziana::Lorentziana(int length) : MedieBlocchi(length){

}

void Lorentziana::CalcolaBlocco(Random &rand){

	// Calcola il valore medio di variabili lorentziane
	for (int i = 0; i < m_length; i++)
		m_blocco += rand.Cauchy(0, 1.);
		
	m_blocco /= (double)m_length;
}

Buffon::Buffon(int length, double L, double d) : MedieBlocchi(length){

	m_L = L;
	m_d = d;
}

void Buffon::CalcolaBlocco(Random &rand){

	double x, y;
	int hit = 0;
	
	for (int i = 0; i < m_length; i++){
		// Genera un numero uniforme nel semicerchio superiore
		do{
			x = rand.Rannyu(-1,1);
			y = rand.Rannyu();
		} while(x*x + y*y > 1);

		// Prendi il coseno e genera un numero casuale in (0,1)
		x /= sqrt(x*x + y*y);
		y = rand.Rannyu();
		
		// Controlla se hai colpito una linea
		if (y == 0 or y == 1.) hit ++;
		y += x * m_L / m_d;
		if (y <= 0 or y >= 1.) hit ++;
	}
	
	m_blocco = 2 * m_L * double(m_length)/((double)hit * m_d);
}
