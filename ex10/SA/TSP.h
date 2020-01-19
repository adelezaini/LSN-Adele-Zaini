
/*******************************************************************************************************************************
Laboratorio Simulazione Numerica
Excercise 9: Genetic Algorithm - Traveling Salesman Problem

 Adele Zaini
**********************************************************************************************************************************/

//NOTE: One can adapt to a generic Genetic Algorithm changing Tour->Chromosome/Individual & City->Gene

#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>
#include <vector>
#include "random.h"

using namespace std;

#ifndef _TSP_h_
#define _TSP_h_

struct point{
   double x,y;
};

class Tour{
  //Note: Tour uses positions vector to evaluate the fitness, but all operations are on the vector of cities indexes

public:
  
  Tour(){};
  Tour(vector<point>); //Create Tour with random city order
  Tour(vector<point>, vector<int>); //Create Tour with given city order
  //OBS: one can implement other constructors (e.g. Tour(), Tour(int ncities)...), but I prefer to keep the restriction on initializing position to avoid potential errors (i.e. evaluating distances...).
  
  ~Tour(){};
  
  //1) Methods working with city order(vector of integer):
  //Generals:
  void Initialize(); //Set random city order
  bool Check()const; //Check if bonds are satified (ALL city just ONCE) (=true)
  int Size()const{return ncities;}
  int Pbc(int); //x_{N+1}=x1 -> working on iterators, risk of going out of the boundaries
  int GetCity(int i)const{return cities[i];}
  //Tour& operator= (const Tour&);
  void Print()const; //Print city order
  void PrintSeed(){rnd->PrintSeed();}
  
  //Mutations of city order:
  void Permutation(); //Permute the whole vector
  void Swap(int,int); //Swap two elements
  void RandomSwap(); //=Pair permutation
  void ShiftTot(); //Shift of n positions
  void Shift(); //Shift of m elements of n positions
  void Inversion(); //Invert elements from i to n
  void InversionTot(); //Invert the whole vector
  
  
  //2) Methods specific of TSP:
  double GetLenght2()const; //give square lenght
  double Distance2(int, int)const; //calculate distance between two cities
  double Fitness()const; //evaluate the fitness (1/L2)
  vector<point> GetPositions()const{return positions;} //give the postion vector
  
private:
  vector<int> cities; //city order
  int ncities;
  
  vector<point> positions; //city positions (untouched by the methods-> Tour works on one-to-one correspondance index<->position --> it does not need changing positions, but it needs them to evaluate lenght and fitness
  Random* rnd;
  
};

//*************************************************************************************************//
/*
class Population{
public:
  
  Population(int, vector<point>); //Create population of tours given the npopulation dimension and the position vector
  ~Population(){};
  
  //Genetic Algorithm:
  void Offspring(double); //Create new generation: 1) elite? 2) crossover with probability Pc 3) no? select and copy a good one
  void Mutation(double); //Change each city order with probability Pm (mutation randomly chosen among RandomSwap, ShiftTot, Shift,Invertion);
  int Select()const; //use the Roussian Roulette method to select individuals with different probability weight (=fitness); //output: tour index
  Tour Crossover(int,int)const; //given the indexes of mum and dad, merge the first part of the mum with the missing elements from the dad in the same order //output: new Tour child
  
  Tour Best()const; //evaluates the fitness for each individual and gives the best tour (highiest fitness)
  double BestFitness()const; //fitness of the best
  double BestLenght2()const; //square lenght of the best
  double BestHalfLenght2(); //evaluates the mean value of fitness of the best half of individuals
  
  Tour GetTour(int i){ return tours[i];}
  
  void Sort(); //uses quicksort to sort the population on a fitness basis ('0'=highiest fitness,...,'npop-1'=lowest fitness)
  void quickSort(int, int); //supports Sort() (it's a recursive method)
  void Swap(int, int); //supports quickSort()
  
private:
  vector<Tour> tours;
  
  int npop;
  int ncities;
  
  double fitness_sum;
  int iter;
  
  Random* rnd;
  
};
*/
#endif
