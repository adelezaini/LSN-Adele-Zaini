CC = g++
CFLAGS = -Wall -O3 --std=c++11

MolDyn_NVE.exe : MolDyn_NVE.o random.o
	$(CC) random.o MolDyn_NVE.o -o MolDyn_NVE.exe
MolDyn_NVE.o : MolDyn_NVE.cpp MolDyn_NVE.h
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean.all:
	rm *.o MolDyn_NVE.exe seed.out old.0 old.final config.final
	bash clean.sh
clean:
	rm *.o seed.out old.0 old.final config.final
