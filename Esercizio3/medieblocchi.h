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

class CallDiretto : public MedieBlocchi {

	public:
		CallDiretto(int length, double mu, double sigma, double T, double S0, double K);
		void CalcolaBlocco(Random &rand);
	
	private:
		double m_mu, m_sigma, m_T, m_S0, m_K;
};

class PutDiretto : public MedieBlocchi {

	public:
		PutDiretto(int length, double mu, double sigma, double T, double S0, double K);
		void CalcolaBlocco(Random &rand);
	
	private:
		double m_mu, m_sigma, m_T, m_S0, m_K;
};

class CallDiscreto : public MedieBlocchi {

	public:
		CallDiscreto(int length, double mu, double sigma, double T, double S0, double K, int steps);
		void CalcolaBlocco(Random &rand);
	
	private:
		double m_mu, m_sigma, m_T, m_deltaT, m_S0, m_K;
		int m_steps;
};

class PutDiscreto : public MedieBlocchi {

	public:
		PutDiscreto(int length, double mu, double sigma, double T, double S0, double K, int steps);
		void CalcolaBlocco(Random &rand);
	
	private:
		double m_mu, m_sigma, m_T, m_deltaT, m_S0, m_K;
		int m_steps;
};

#endif
