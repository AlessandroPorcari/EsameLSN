CFLAGS = -Wall -O3 --std=c++11

Esercizio10.1 : Esercizio10.1.o random.o cromosomi.o
	g++ random.o Esercizio10.1.o cromosomi.o -o Esercizio10.1
	
Esercizio10.1.o : Esercizio10.1.cpp
	g++ -c Esercizio10.1.cpp -o Esercizio10.1.o $(CFLAGS)
	
cromosomi.o : cromosomi.cpp cromosomi.h
	g++ -c cromosomi.cpp -o cromosomi.o $(CFLAGS)	
	
random.o : ../../Random/random.cpp ../../Random/random.h
	g++ -c ../../Random/random.cpp -o random.o $(CFLAGS)
	
clean :
	rm *.o Esercizio10.1
	
cleantxt :
	rm *.txt
