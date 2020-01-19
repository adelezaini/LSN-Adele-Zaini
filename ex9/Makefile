CC = g++
CFLAGS = -Wall -O3 --std=c++11

main.exe : main.o random.o tour.o population.o
	$(CC) random.o tour.o population.o main.o -o main.exe
main.o : main.cpp random.h TSP.h
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
tour.o : tour.cpp TSP.h random.h
	$(CC) -c tour.cpp -o tour.o $(CFLAGS)
population.o : population.cpp TSP.h random.h
	$(CC) -c population.cpp -o population.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
