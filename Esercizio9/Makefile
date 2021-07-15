CFLAGS = -Wall -O3 --std=c++11

Esercizio9 : Esercizio9.o random.o cromosomi.o
	g++ random.o Esercizio9.o cromosomi.o -o Esercizio9
	
Esercizio9.o : Esercizio9.cpp
	g++ -c Esercizio9.cpp -o Esercizio9.o $(CFLAGS)
	
cromosomi.o : cromosomi.cpp cromosomi.h
	g++ -c cromosomi.cpp -o cromosomi.o $(CFLAGS)	
	
random.o : ../Random/random.cpp ../Random/random.h
	g++ -c ../Random/random.cpp -o random.o $(CFLAGS)
	
clean :
	rm *.o Esercizio9
	
cleantxt :
	rm *.txt
