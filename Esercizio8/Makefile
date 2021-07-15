CFLAGS = -Wall -O3 --std=c++11

Esercizio8 : Esercizio8.o random.o medieblocchi.o
	g++ random.o Esercizio8.o medieblocchi.o -o Esercizio8
	
Esercizio8.o : Esercizio8.cpp
	g++ -c Esercizio8.cpp -o Esercizio8.o $(CFLAGS)
	
medieblocchi.o : medieblocchi.cpp medieblocchi.h
	g++ -c medieblocchi.cpp -o medieblocchi.o $(CFLAGS)
	
random.o : ../Random/random.cpp ../Random/random.h
	g++ -c ../Random/random.cpp -o random.o $(CFLAGS)
	
clean :
	rm *.o Esercizio8
	
cleantxt :
	rm *.txt
