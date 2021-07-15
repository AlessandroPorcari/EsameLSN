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
		void PrintBlocco(string nome);
		
		virtual void CalcolaBlocco(Random &rand) = 0;		

	protected:
	
		int m_length;
		double m_somma, m_somma2;
		double m_errore, m_blocco;
};

class MediaUnif : public MedieBlocchi {

	public:
		MediaUnif(int length);
		void CalcolaBlocco(Random &rand);
	
};

class VarUnif : public MedieBlocchi {

	public:
		VarUnif(int length);
		void CalcolaBlocco(Random &rand);

};

class ChiQuadro : public MedieBlocchi {

	public:
		ChiQuadro(int length, int M);
		~ChiQuadro();
		void CalcolaBlocco(Random &rand);
	
	private:
		int m_M; // Numero dei sottointervalli di [0,1]
		int *m_counter; // Eventi divisi nei sottointervalli
		
};

class Esponenziale : public MedieBlocchi {

	public:
		Esponenziale(int length);
		void CalcolaBlocco(Random &rand);

};

class Lorentziana : public MedieBlocchi {

	public:
		Lorentziana(int length);
		void CalcolaBlocco(Random &rand);
	
};


class Buffon : public MedieBlocchi {

	public:
		Buffon(int length, double L, double d);
		void CalcolaBlocco(Random &rand);
		
	private:
		double m_L; // Lunghezza dell'ago
		double m_d; // Distanza tra linee

};

#endif
