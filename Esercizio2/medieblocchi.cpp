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

Uniforme::Uniforme(int length) : MedieBlocchi(length){

}

void Uniforme::CalcolaBlocco(Random &rand){

	// Calcola l'integrale di f(x) = PI/2 * cos(x*PI/2) con distribuzione uniforme 
	for (int i = 0; i < m_length; i++)
		m_blocco += cos(M_PI * rand.Rannyu() / 2);
		
	m_blocco = m_blocco * M_PI / ((double)m_length * 2);
}

Lineare::Lineare(int length) : MedieBlocchi(length){

}

void Lineare::CalcolaBlocco(Random &rand){

	// Calcola l'integrale di f(x) = PI/2 * cos(x*PI/2) con distribuzione p(x) = 2(1-x) 
	double x;
	for (int i = 0; i < m_length; i++){
		x = 1 - sqrt(1 - rand.Rannyu());
		m_blocco += cos(M_PI * x / 2) / (1 - x);
	}
		
	m_blocco = m_blocco * M_PI / ((double)m_length * 4);
}

MedieBlocchiVett::MedieBlocchiVett(int length, int size){

	m_length = length;
	m_size = size;
	m_somma = new double[m_size];
	m_somma2 = new double[m_size];
	m_blocco = new double[m_size];
	m_errore = new double[m_size];
	
	for (int i = 0; i < m_size; i++){
		m_somma[i] = 0;
		m_somma2[i] = 0;
	}
	
}

MedieBlocchiVett::~MedieBlocchiVett(){

	delete [] m_somma;
	delete [] m_somma2;
	delete [] m_blocco;
	delete [] m_errore;

}

void MedieBlocchiVett::GetSomme(){

	for (int i = 0; i < m_size; i++){
		m_somma[i] += m_blocco[i];
		m_somma2[i] += m_blocco[i] * m_blocco[i];
	}
}

void MedieBlocchiVett::Reset(){

	for (int i = 0; i < m_size; i++)
		m_blocco[i] = 0;
}

void MedieBlocchiVett::Results(string nome, int num){

	// Calcola gli errori
	if (num == 1){
		for (int i = 0; i < m_size; i++)
			m_errore[i] = 0;
	}
	else {
		for (int i = 0; i < m_size; i++)
			m_errore[i] = sqrt((m_somma2[i]/(double)num-pow(m_somma[i]/(double)num,2))/((double)num-1));
	}

	// L'errore del random walk di dimensione 1 deve essere nullo perchè avendo fatto
	// un solo passo la distanza è sempre uguale. Se non metto la seguente riga ottengo
	// nan nel random walk continuo. Questo è dovuto alle approssimazioni fatte nel
	// calcolo del seno e del coseno che modificano leggermente la distanza. Da questo
	// viene che in certi casi l'errore nel random walk di lunghezza 1 non è più nullo
	// ma diventa una radice di un numero negativo molto piccolo (sqrt(10^-15))
	m_errore[0] = 0;

	// Stampa a video
	ofstream output;
	output.open(nome);
	for (int i = 0; i < m_size; i++)
		output << i+1 << " " << m_somma[i]/(double)num << " " << m_errore[i] << endl;	
	output.close();
}

PassoReticolo::PassoReticolo(int length, int size, double passo) : MedieBlocchiVett(length, size){

	m_passo = passo;
	
}

void PassoReticolo::CalcolaBlocco(Random &rand){
	
	int segno;
		
	for (int i = 0; i < m_length; i++){
		// Inizia dall'origine
		for (int i = 0; i < 3; i++)
			m_posizione[i] = 0;
		
		// Genera un random walk
		for (int j = 0; j < m_size; j++){
			// Fai un passo in una direzione
			if (rand.Rannyu() > 0.5) segno = 1;
			else segno = -1;
			m_posizione[int(3*rand.Rannyu())] += m_passo * segno;
			
			// Calcola distanza al quadrato
			for (int k = 0; k < 3; k++)
				m_blocco[j] += m_posizione[k] * m_posizione[k];
		}
	}
	
	// Calcola il valore della distanza quadratica media del blocco
	for (int i = 0; i < m_size; i++)
		m_blocco[i] = sqrt(m_blocco[i] / (double)m_length);
			
}

PassoIsotropo::PassoIsotropo(int length, int size, double passo) : MedieBlocchiVett(length, size){

	m_passo = passo;
	
}

void PassoIsotropo::CalcolaBlocco(Random &rand){
	
	double theta, phi;
		
	for (int i = 0; i < m_length; i++){
		// Inizia dall'origine
		for (int i = 0; i < 3; i++)
			m_posizione[i] = 0;
			
		// Genera un random walk
		for (int j = 0; j < m_size; j++){
			// Fai un passo in una direzione casuale
			theta = rand.Rannyu(0,M_PI);
			phi = rand.Rannyu(0,2*M_PI);
			m_posizione[0] += m_passo * sin(theta) * cos(phi);
			m_posizione[1] += m_passo * sin(theta) * sin(phi);
			m_posizione[2] += m_passo * cos(theta);

			// Calcola distanza al quadrato
			for (int k = 0; k < 3; k++)
				m_blocco[j] += m_posizione[k] * m_posizione[k];
		}
	}
	
	// Calcola il valore della distanza quadratica media del blocco
	for (int i = 0; i < m_size; i++)
		m_blocco[i] = sqrt(m_blocco[i] / (double)m_length);
			
}
