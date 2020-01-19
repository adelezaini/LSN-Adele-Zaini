/*******************************************************************************************************************************
Laboratorio Simulazione Numerica
Excercise 9: Genetic Algorithm - Traveling Salesman Problem

 Adele Zaini
**********************************************************************************************************************************/

#include "TSP.h"

using namespace std;

/*------------------------------------------------------------------------//
                                TOUR
//------------------------------------------------------------------------*/

//*********************************************************//
//CREATORS:
Tour::Tour(vector<point> pos){ // pos_initialized=true;
  
  positions=pos;
  ncities=pos.size();
  
  Initialize();
  
  rnd=new Random();
  rnd ->SetSeed("seed.out");
    rnd->SaveSeed();
}

Tour::Tour(vector<point> pos, vector<int> cit){ //pos_initialized=true;
  
  cities=cit;
  positions=pos;
  ncities=pos.size();
  
  rnd=new Random();
  rnd ->SetSeed("seed.out");
    rnd->SaveSeed();
}

void Tour::Initialize(){
  for (int i=0; i<ncities; i++) cities.push_back(i); //0 1 2 3 4 5 6 7 8 9 ...

  Permutation(); //permutation of initial configuration
}

//*********************************************************//
//CHECK:
bool Tour::Check()const{ //check if ALL cities are visited ONLY ONCE
  //"ALL"-> for loop checking integer in [1,ncities]
  //"ONLY ONCE"->it's implicit, since the vector size is ncities and all have to be visited (if one city were visited twice then some other wouldnt be in the tour)
  bool check=true;
  for(int i=0; i<ncities; i++){
    auto it=find(cities.begin(),cities.end(),i);
    if(it==cities.end()){
      cerr << "PROBLEM: Tour does not satisfy bounds." <<endl; //Fill it with random sorting."
      check=false;
      return check;
    }
  }
  return check;
}
/*void Check(void){
for(int i=1; i<=ncities; i++){
  vector<int>::iterator it=find(cities.begin(),cities.end(),i);
  if(it==cities.end()){
    cout << "PROBLEM: Tour does not satisfy bounds. Fill it with random sorting." <<endl;
    Initialize();
    return;
  }
 }
}*/
//*********************************************************//
//FITNESS:
double Tour::GetLenght2()const{ //just if pos initialized
  double sum=0.;
  double n=ncities-1;
  
  for (int i=0; i<n; i++){
    int city1=cities[i];
    int city2=cities[i+1];
    sum+=Distance2(city1, city2); //(positions[i], positions[i+1]);
  }
  sum+=Distance2(cities[0], cities[ncities-1]); //(positions.front(), positions.back());
  
  return sum;
}

double Tour::Distance2(int a, int b)const{//Distance between i-th city and j-th city
  point p1=positions[a];//(point p1, point p2)const{
  point p2=positions[b];
  return (p1.x-p2.x)*(p1.x-p2.x)+(p1.y-p2.y)*(p1.y-p2.y);
}

double Tour::Fitness()const{
  double L2=GetLenght2();
  double fitness=1./L2;
  return fitness;
}

//*********************************************************//
//MUTATIONS:
void Tour::Swap(int i, int j){ //Swap of two cities
  swap(cities[Pbc(i)], cities[Pbc(j)]);
}

void Tour::Permutation(){//total permutation of elements
  random_shuffle(cities.begin(), cities.end());
}

void Tour::RandomSwap(){ //=pair permutation
  int i=rnd->Rannyu(0,ncities);
  int j=rnd->Rannyu(0,ncities);
  Swap(i,j);
}

void Tour::ShiftTot(){ //Shift of +n position of all vector
  int n=rnd->Rannyu(1,ncities-1);
  rotate(cities.begin(), cities.begin()+n, cities.end());
  //middle becomes the new first (middle-1 the new last)
}

void Tour::Shift(){ //Shift of +n position for m cities (i=index where to start counting m elemtn)
  //OBS: not exactly same proposed shift but analogous
  int n=rnd->Rannyu(1,ncities);
  int i=rnd->Rannyu(0,ncities);
  int m=rnd->Rannyu(2,ncities); //solo sottogruppi: [2,ncities-1]
  int en=n+m-1;
  for (int k=0; k<m; k++){
    for (int l=0; l<en; l++) Swap(Pbc(i+l), Pbc(i+l+1));
  }
  /*int b=Pbc(i+m); //otherwise segmentation fault
  int c=Pbc(n+i+m);
    //idea: a) elements before i are not influenced
  //b) for the moment also m elements not influenced and rotate RIGHT of n+m
  rotate(cities.begin()+b, cities.begin()+c, cities.end());
  //c) now rotate LEFT of n all the (i, ncities) elements
  rotate(cities.begin()+i, cities.begin() +cities.size()-n, cities.end());*/
}

void Tour::Inversion(){
  int i=rnd->Rannyu(0,ncities-1);
  int n=rnd->Rannyu(i,ncities);
  for (int k=i; k<n; k++) Swap(k,n-k);
}

void Tour::InversionTot(){
  int M=0.5*ncities;
  for (int k=0; k<M; k++) Swap(k,ncities-1-k);
}
  
//*********************************************************//
int Tour::Pbc(int i){
  while(i<0) i=i+ncities;
  while(i>=ncities) i=i-ncities;
  return i;
}

void Tour::Print()const{
  for (int i=0; i<ncities; ++i) cout << cities[i] << " ";
  cout << endl;
}

Tour& Tour::operator=(const Tour& tour2){
  cities=tour2.cities;
  positions=tour2.positions;

  return *this;
}
