CFLAGS = -Wall -O3 --std=c++11

Esercizio5 : Esercizio5.o random.o medieblocchi.o
	g++ random.o Esercizio5.o medieblocchi.o -o Esercizio5
	
Esercizio5.o : Esercizio5.cpp
	g++ -c Esercizio5.cpp -o Esercizio5.o $(CFLAGS)
	
medieblocchi.o : medieblocchi.cpp medieblocchi.h
	g++ -c medieblocchi.cpp -o medieblocchi.o $(CFLAGS)
	
random.o : ../Random/random.cpp ../Random/random.h
	g++ -c ../Random/random.cpp -o random.o $(CFLAGS)
	
clean :
	rm *.o Esercizio5
	
cleantxt :
	rm *.txt
