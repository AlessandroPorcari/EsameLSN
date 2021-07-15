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

class Uniforme100 : public MedieBlocchi {

	public:
		Uniforme100(int length, double step, double position[3]);
		void CalcolaBlocco(Random &rand);
		void PrintPosizione(string nome);
	
	private:
		double m_step;
		double m_position[3];
		double m_r;
};

class Uniforme210 : public MedieBlocchi {

	public:
		Uniforme210(int length, double step, double position[3]);
		void CalcolaBlocco(Random &rand);
		void PrintPosizione(string nome);
	
	private:
		double m_step;
		double m_position[3];
		double m_r;
};

class Gauss100 : public MedieBlocchi {

	public:
		Gauss100(int length, double step, double position[3]);
		void CalcolaBlocco(Random &rand);
		void PrintPosizione(string nome);
	
	private:
		double m_step;
		double m_position[3];
		double m_r;
};

class Gauss210 : public MedieBlocchi {

	public:
		Gauss210(int length, double step, double position[3]);
		void CalcolaBlocco(Random &rand);
		void PrintPosizione(string nome);
	
	private:
		double m_step;
		double m_position[3];
		double m_r;
};

#endif
