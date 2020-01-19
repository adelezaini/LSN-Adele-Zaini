#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <vector>
#include "random.h"

using namespace std;

//*************************************************************************************************//
//MAIN.h
//*************************************************************************************************//

//random numbers
#include "random.h"
Random rnd;

//cities
string config;
int ncities;

//questo vettore contiene posizioni con corrispondenza biunivoca tra: indice nel vettore – posizione reale
//es: [4.5, 8.3, 9.4, 1.2] --> indice 3 <-> 9.4
vector<point> positions;

Tour x_tour;
Tour best_tour;
double best_lenght;

//simulation
double temp, t_rate, beta;
int tstep, nstep_final, nstep=0;
int tcount=0;
double alpha,e;
int accepted=0;
int attempted=0;

//output
string folder;
ofstream outL;

//functions
void Input(void);
void Check(Tour);

vector<point> Circonference(void);
vector<point> Square(void);

void Reset(void);
Tour Mutate(Tour);
void Metropolis(void);
void FindBest(void);
void NextTemperature(void);

void PrintPositions(void);
void PrintBestLenght2(int);
void PrintBest(void);
void PrintConsole();

void Clean(void);
