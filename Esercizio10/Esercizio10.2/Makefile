CFLAGS = -Wall -O3 --std=c++11

Esercizio10.2 : Esercizio10.2.o random.o cromosomi.o
	mpicxx random.o Esercizio10.2.o cromosomi.o -o Esercizio10.2
	
Esercizio10.2.o : Esercizio10.2.cpp
	mpicxx -c Esercizio10.2.cpp -o Esercizio10.2.o $(CFLAGS)
	
cromosomi.o : cromosomi.cpp cromosomi.h
	mpicxx -c cromosomi.cpp -o cromosomi.o $(CFLAGS)	
	
random.o : ../../Random/random.cpp ../../Random/random.h
	g++ -c ../../Random/random.cpp -o random.o $(CFLAGS)
	
clean :
	rm *.o Esercizio10.2
	
cleantxt :
	rm *.txt
