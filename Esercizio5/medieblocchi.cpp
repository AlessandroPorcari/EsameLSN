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

Uniforme100::Uniforme100(int length, double step, double position[3]) : MedieBlocchi(length){

	m_step = step;
	m_r = 0;
	for (int i = 0; i < 3; i++){
		m_position[i] = position[i];
		m_r += m_position[i]*m_position[i];
	}
	
	m_r = sqrt(m_r);

}

void Uniforme100::CalcolaBlocco(Random &rand){

	double new_pos[3];
	double r_new;
	double prob;
	int accepted = 0;
	
	for (int i = 0; i < m_length; i++){
		// Calcola una nuova posizione
		r_new = 0;
		for (int j = 0; j < 3; j++){
			new_pos[j] = m_position[j] + m_step * rand.Rannyu(-1, 1);
			r_new += new_pos[j]* new_pos[j];
		}
		r_new = sqrt(r_new);
		
		// Usa l'algoritmo di Metropolis
		prob = exp(2.*(m_r - r_new));
		if (prob > rand.Rannyu()){
			accepted ++;
			m_r = r_new;
			for (int j = 0; j < 3; j++)
				m_position[j] = new_pos[j];
		}
		m_blocco += m_r;
		
		// Togliere commento per stampare su file la posizione ogni passo
		if ((i+1)%10 == 0)
			PrintPosizione("pos_unif_100.txt");
	}
		
	m_blocco /= (double)m_length;
	// Togliere il commento se voglio controllare accettazione
	//cout << "Acceptance rate = " << (double)accepted / (double)m_length << endl;
}

void Uniforme100::PrintPosizione(string nome){

	ofstream output;
	output.open(nome,ios::app);
	output << m_position[0] << " " << m_position[1] << " " << m_position[2] << endl;	
	output.close();

}

Uniforme210::Uniforme210(int length, double step, double position[3]) : MedieBlocchi(length){

	m_step = step;
	m_r = 0;
	for (int i = 0; i < 3; i++){
		m_position[i] = position[i];
		m_r += m_position[i]*m_position[i];
	}

	m_r = sqrt(m_r);

}

void Uniforme210::CalcolaBlocco(Random &rand){

	double new_pos[3];
	double r_new;
	double prob;
	int accepted = 0;
	
	for (int i = 0; i < m_length; i++){
		// Calcola una nuova posizione
		r_new = 0;
		for (int j = 0; j < 3; j++){
			new_pos[j] = m_position[j] + m_step * rand.Rannyu(-1, 1);
			r_new += new_pos[j]* new_pos[j];
		}
		r_new = sqrt(r_new);
		
		// Usa l'algortimo di Metropolis
		prob = exp(m_r - r_new)*pow(new_pos[2]/m_position[2],2);
		if (prob > rand.Rannyu()){
			accepted ++;
			m_r = r_new;
			for (int j = 0; j < 3; j++)
				m_position[j] = new_pos[j];
		}
		m_blocco += m_r;
		
		// Togliere commento per stampare su file la posizione ogni passo
		if ((i+1)%10 == 0)
			PrintPosizione("pos_unif_210.txt");
	}
		
	m_blocco /= (double)m_length;
	// Togliere il commento se voglio controllare accettazione
	//cout << "Acceptance rate = " << (double)accepted / (double)m_length << endl;
}

void Uniforme210::PrintPosizione(string nome){

	ofstream output;
	output.open(nome,ios::app);
	output << m_position[0] << " " << m_position[1] << " " << m_position[2] << endl;	
	output.close();

}

Gauss100::Gauss100(int length, double step, double position[3]) : MedieBlocchi(length){

	m_step = step;
	m_r = 0;
	for (int i = 0; i < 3; i++){
		m_position[i] = position[i];
		m_r += m_position[i]*m_position[i];
	}
	
	m_r = sqrt(m_r);

}

void Gauss100::CalcolaBlocco(Random &rand){

	double new_pos[3];
	double r_new;
	double prob;
	int accepted = 0;
	
	for (int i = 0; i < m_length; i++){
		// Calcola una nuova posizione
		r_new = 0;
		for (int j = 0; j < 3; j++){
			new_pos[j] = m_position[j] + m_step * rand.Gauss(0, 1);
			r_new += new_pos[j]* new_pos[j];
		}
		r_new = sqrt(r_new);
		
		// Usa l'algortimo di Metropolis
		prob = exp(2.*(m_r - r_new));
		if (prob > rand.Rannyu()){
			accepted ++;
			m_r = r_new;
			for (int j = 0; j < 3; j++)
				m_position[j] = new_pos[j];
		}
		m_blocco += m_r;
		
		// Togliere commento per stampare su file la posizione ogni passo
		if ((i+1)%10 == 0)
			PrintPosizione("pos_gauss_100.txt");
	}
		
	m_blocco /= (double)m_length;
	// Togliere il commento se voglio controllare accettazione
	//cout << "Acceptance rate = " << (double)accepted / (double)m_length << endl;
}

void Gauss100::PrintPosizione(string nome){

	ofstream output;
	output.open(nome,ios::app);
	output << m_position[0] << " " << m_position[1] << " " << m_position[2] << endl;	
	output.close();

}

Gauss210::Gauss210(int length, double step, double position[3]) : MedieBlocchi(length){

	m_step = step;
	m_r = 0;
	for (int i = 0; i < 3; i++){
		m_position[i] = position[i];
		m_r += m_position[i]*m_position[i];
	}
	
	m_r = sqrt(m_r);

}

void Gauss210::CalcolaBlocco(Random &rand){

	double new_pos[3];
	double r_new;
	double prob;
	int accepted = 0;
	
	for (int i = 0; i < m_length; i++){
		// Calcola una nuova posizione
		r_new = 0;
		for (int j = 0; j < 3; j++){
			new_pos[j] = m_position[j] + m_step * rand.Gauss(0, 1);
			r_new += new_pos[j]* new_pos[j];
		}
		r_new = sqrt(r_new);
		
		// Usa l'algortimo di Metropolis
		prob = exp(m_r - r_new)*pow(new_pos[2]/m_position[2],2);
		if (prob > rand.Rannyu()){
			accepted ++;
			m_r = r_new;
			for (int j = 0; j < 3; j++)
				m_position[j] = new_pos[j];
		}
		m_blocco += m_r;
		
		// Togliere commento per stampare su file la posizione ogni passo
		if ((i+1)%10 == 0)
			PrintPosizione("pos_gauss_210.txt");
	}
		
	m_blocco /= (double)m_length;
	// Togliere il commento se voglio controllare accettazione
	//cout << "Acceptance rate = " << (double)accepted / (double)m_length << endl;
}

void Gauss210::PrintPosizione(string nome){

	ofstream output;
	output.open(nome,ios::app);
	output << m_position[0] << " " << m_position[1] << " " << m_position[2] << endl;	
	output.close();

}
