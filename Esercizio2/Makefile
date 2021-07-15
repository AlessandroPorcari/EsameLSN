CFLAGS = -Wall -O3 --std=c++11

Esercizio2 : Esercizio2.o random.o medieblocchi.o
	g++ random.o Esercizio2.o medieblocchi.o -o Esercizio2
	
Esercizio2.o : Esercizio2.cpp
	g++ -c Esercizio2.cpp -o Esercizio2.o $(CFLAGS)
	
medieblocchi.o : medieblocchi.cpp medieblocchi.h
	g++ -c medieblocchi.cpp -o medieblocchi.o $(CFLAGS)
	
random.o : ../Random/random.cpp ../Random/random.h
	g++ -c ../Random/random.cpp -o random.o $(CFLAGS)
	
clean :
	rm *.o Esercizio2
	
cleantxt :
	rm *.txt
