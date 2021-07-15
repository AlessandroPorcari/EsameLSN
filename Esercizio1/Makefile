CFLAGS = -Wall -O3 --std=c++11

Esercizio1 : Esercizio1.o random.o medieblocchi.o
	g++ random.o Esercizio1.o medieblocchi.o -o Esercizio1
	
Esercizio1.o : Esercizio1.cpp
	g++ -c Esercizio1.cpp -o Esercizio1.o $(CFLAGS)
	
medieblocchi.o : medieblocchi.cpp medieblocchi.h
	g++ -c medieblocchi.cpp -o medieblocchi.o $(CFLAGS)
	
random.o : ../Random/random.cpp ../Random/random.h
	g++ -c ../Random/random.cpp -o random.o $(CFLAGS)
	
clean :
	rm *.o Esercizio1
	
cleantxt :
	rm *.txt
