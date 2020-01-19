/*******************************************************************************************************************************
Laboratorio Simulazione Numerica
Excercise 9: Genetic Algorithm - Traveling Salesman Problem

 Adele Zaini
**********************************************************************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <vector>
#include "random.h"
#include "TSP.h"


using namespace std;

//NOTE: when I wrote this code I prefered to keep it as "clean" as possible to better understand the operational functions and the logical idea behind. At the same time Id liked to keep it open to improvements and future re-adaptions.

//*************************************************************************************************//
//MAIN.h
//*************************************************************************************************//
//input
string config;
int ncities, npop, ngeneration; // dim_elite;
double Pc, Pm;

//questo vettore contiene posizioni con corrispondenza biunivoca tra: indice nel vettore – posizione reale
//es: [4.5, 8.3, 9.4, 1.2] --> indice 3 <-> 9.4
vector<point> positions;
Population * population;

//output
string folder;
ofstream outL;

//functions
void Input(void);
vector<point> Circonference(void);
vector<point> Square(void);
void PrintPositions(void);
void PrintBestLenght2(int);
void PrintBest(void);
void Check(Population*, int);

Random rnd;

//*************************************************************************************************//
//MAIN.cpp
//*************************************************************************************************//
int main(){
  
  Input();
  
  //Genetic Algorithm Run: con condizione di stop
  for (int g=0; g<ngeneration; g++){//magari posso anche mettere altra condizione
    if (g%50==0) cout << g<< endl;
    Check(population, g);
    population->Offspring(Pc);
    population->Mutation(Pm);
    PrintBestLenght2(g); //unico file con 0=ngen, 1=best, 2=besthalf
  }
  
  PrintBest();
  
  cout << endl<<"Done √"<<endl;
  cout << "After "<<ngeneration<<" generations, the best lenght is = "<<population->BestLenght2()<<endl<<endl;
    cout <<endl;
  
  outL.close();
  rnd.SaveSeed();
  
  return 0;
}
  
//*************************************************************************************************//
//INITIALIZATION:

void Input(void){
  
  //Initialize random variable --> Tour and Population will enter data from 'seed.out'
  rnd.SetSeed();
  rnd.SaveSeed();
  
  ifstream ReadInput;

  cout << "---------------------------------------------------------" <<endl;
  cout << "          Travelling Salesman Problem     " << endl;
  cout << "              Genetic Algorithm                 "<<endl;
  cout << "---------------------------------------------------------" <<endl<<endl;
  
  ReadInput.open("input.dat");
  
//Set tour:
  ReadInput >> ncities;
  ReadInput >> config;
  
  cout << "The program evaluates the best tour (highest fitness)"<<endl;
  cout << "-> among " << ncities << " cities placed ";
  
  //Set positions of cities
  if (config=="Circonference"){
    positions=Circonference();
    folder="Circonference/";
    cout << "on a circonference of r=1" << endl;
  }else if (config=="Square"){
    positions=Square();
    folder="Square/";
    cout << "within a square of l=1" << endl;
  }else{
    cerr << "PROBLEM: this program works only with cities placed on a circonference or within a square."<<endl<<"Plese choose 'configuration input' between 'Circonference' or 'Square'."<< endl<<endl;
      exit(1);
    }
  //OSS: si possono aggiungere configurazioni qui
  PrintPositions();
  

//Set Genetic Algorithm parameters:
  cout <<endl<< "Simulation parameters:"<<endl;
  
  ReadInput >> npop;
  ReadInput >> ngeneration;
  
  ReadInput >> Pc;
  ReadInput >> Pm;
  
  cout << "– Number of individuals in each population = " << npop << endl;
  cout << "– Number of generations = " << ngeneration << endl;
  cout << "– Probability of Crossover= " << Pc << endl;
  cout << "– Probability of Mutation= " << Pm << endl;
  
    ReadInput.close();
  
//Create population of random tours
  population=new Population(npop,positions);
  
  cout <<endl<< "Population initialized √ -> Initial lenght = "<<population->BestLenght2()<<endl;
  
  outL.open(folder+"Best_Lenght.out");
}
//*************************************************************************************************//
//ACTUAL POSITIONS OF CITIES:

//Create vector of 'ncities' points located:
//a) on a circonference (r=1):
vector<point> Circonference(void){
  
  vector<point> cities_pos;
  for(int i=0;i<ncities;i++){
    double theta=rnd.Angle();
    double x=cos(theta);
    double y=sin(theta);
    point position={x,y};
    cities_pos.push_back(position);
  }
  return cities_pos;
  }

//b) within a square(2x2):
vector<point> Square(void){
  
  vector<point> cities_pos;
  for(int i=0;i<ncities;i++){
    double x=rnd.Rannyu(-0.5,0.5);
    double y=rnd.Rannyu(-0.5,0.5);
    point position={x,y};
    cities_pos.push_back(position);
  }
  return cities_pos;
}
//OBS: one can implement the possibility to choose the ray/side lenght, but now I prefer to keep it simple
//*************************************************************************************************//
//CHECK POPULATION:

void Check(Population* pop, int gen){
bool check=true;

  for (int i=0; i<npop; i++){
    Tour tour=pop->GetTour(i);
    check=tour.Check();
    if(!check){
      cerr <<"PROBLEM (Tour bonds) in generation n."<<gen<<endl; //già messaggio di errore in Tour
      exit(-1);
    }
  }
  
}
//*************************************************************************************************//
//PRINT:
void PrintPositions(void){
  ofstream Config(folder+"Positions.out");
  for (int i=0; i<ncities; i++){
    Config << i << " " << positions[i].x << " " << positions[i].y <<endl;
  }
  Config.close();
}

void PrintBestLenght2(int gen){
  outL << gen << " " << population->BestLenght2() << " " << population->BestHalfLenght2() <<endl;
}

void PrintBest(void){
  ofstream outB(folder+"Best_Tour.out");
  outB << ngeneration <<endl;
  
  Tour Best=population->Best();
  int size=Best.Size(); //sicuramente =ncities, ma per sicurezza lo tengo dipendente
  for (int i=0; i<size; i++){
    int city_index=Best.GetCity(i);
    outB <<city_index << " " << positions[city_index].x <<" " <<positions[city_index].y <<endl;
  }
  outB.close();
}

