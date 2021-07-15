#ifndef __medieblocchi_h__
#define __medieblocchi_h__

#include <string>
#include <cmath>
#include <fstream>
#include <iostream>

#include "../Random/random.h"

using namespace std;

Random generatore();

class MedieBlocchi {

	public:
	
		MedieBlocchi(int length);
		
		void GetSomme();
		void Results(string nome, int num);
		void Reset() { m_blocco = 0; }
		
		virtual void CalcolaBlocco(Random &rand) = 0;		

	protected:
	
		int m_length;
		double m_somma, m_somma2;
		double m_errore, m_blocco;
};

class Metropolis : public MedieBlocchi {

	public:
		Metropolis(int length, double step, double position, double mu, double sigma);
		void CalcolaBlocco(Random &rand);
		void PrintPosizione(string nome);
		double Energia() {return m_blocco; }
	
	private:
		double m_step;
		double m_position;
		double m_mu;
		double m_sigma;
		double m_prob;
		
		double GetProbability(double pos);
		double Integranda(double pos);
};

void CalcolaGriglia(double griglia [4], int punti [2], int length, double posizione, Random & rand);

#endif
