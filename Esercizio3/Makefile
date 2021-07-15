CFLAGS = -Wall -O3 --std=c++11

Esercizio3 : Esercizio3.o random.o medieblocchi.o
	g++ random.o Esercizio3.o medieblocchi.o -o Esercizio3
	
Esercizio3.o : Esercizio3.cpp
	g++ -c Esercizio3.cpp -o Esercizio3.o $(CFLAGS)
	
medieblocchi.o : medieblocchi.cpp medieblocchi.h
	g++ -c medieblocchi.cpp -o medieblocchi.o $(CFLAGS)
	
random.o : ../Random/random.cpp ../Random/random.h
	g++ -c ../Random/random.cpp -o random.o $(CFLAGS)
	
clean :
	rm *.o Esercizio3
	
cleantxt :
	rm *.txt
