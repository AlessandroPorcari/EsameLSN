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

class Uniforme : public MedieBlocchi {

	public:
		Uniforme(int length);
		void CalcolaBlocco(Random &rand);
	
};

class Lineare : public MedieBlocchi {

	public:
		Lineare(int length);
		void CalcolaBlocco(Random &rand);
	
};

class MedieBlocchiVett {

	public:
	
		MedieBlocchiVett(int length, int size);
		~MedieBlocchiVett();
		
		void GetSomme();
		void Results(string nome, int num);
		void Reset();
		
		virtual void CalcolaBlocco(Random &rand) = 0;		

	protected:
	
		int m_length;
		int m_size;
		double *m_somma, *m_somma2;
		double *m_errore, *m_blocco;
};

class PassoReticolo : public MedieBlocchiVett {

	public:
		PassoReticolo(int length, int size, double passo);
		void CalcolaBlocco(Random &rand);
	
	private:
		double m_passo;
		double m_posizione[3];
};

class PassoIsotropo : public MedieBlocchiVett {

	public:
		PassoIsotropo(int length, int size, double passo);
		void CalcolaBlocco(Random &rand);
	
	private:
		double m_passo;
		double m_posizione[3];
};

#endif
