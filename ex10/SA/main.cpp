/*******************************************************************************************************************************
Laboratorio Simulazione Numerica
Excercise 10: Simulated Annealing – Traveling Salesman Problem

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
#include "main.h"

using namespace std;


//*************************************************************************************************//
//MAIN

int main(){
  
  Input();

  for(int t=0; t<tstep; t++){
    Reset();
    for(int l=0; l<nstep; l++) {
      Metropolis();
    }
    PrintConsole();
    
    FindBest();
    PrintBestLenght2(t);
    NextTemperature();
  }
  
  PrintBest();
    
  Clean();

  return 0;
}

//*************************************************************************************************//
//SET SIMULATION PARAMETERS:

void Input(void){
    
    //Initialize random variable --> Tour will enter data from 'seed.out'
  //Firstly I check if there's a seed.out file to restart simulation and get different configurations
  ifstream seedout("seed.out");
  bool checkseed=seedout.is_open();
  seedout.close();
  
  if(checkseed) rnd.SetSeed("seed.out");
  else rnd.SetSeed();
  
    rnd.SaveSeed();
    
    ifstream ReadInput;
    
    cout << "---------------------------------------------------------" <<endl;
    cout << "          Travelling Salesman Problem     " << endl;
    cout << "              Simulated Annealing                 "<<endl;
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
  
  //Create initial configuration with city randomly sorted:
  Tour tour(positions);
  x_tour=tour;
    Check(x_tour);
  
    best_tour=x_tour;
    best_lenght=x_tour.GetLenght2();
    cout <<endl<< "Tour initialized √ -> Initial lenght = "<<best_lenght<<endl<<endl;
    
    
    //Set Simulated Annealing parameters:
    cout <<endl<< "Simulated Annealing parameters:"<<endl;
    
    ReadInput >> temp;
    ReadInput >> t_rate;
    ReadInput >> tstep;
    
    ReadInput >> nstep_final;
    
    cout << "– Initial temperature T = " << temp << endl;
    cout << "- Reduction rate of temperature α=" << t_rate <<endl;
    cout << "– Dimension of annealing schedule N = " << tstep << endl;
    
    
    cout << "At a given T, the program uses the Metropolis algorithm to sample configurations according to the Boltzmann distribution." << endl;
    cout << "– Number of Metropolis steps n = " << nstep_final << endl<<endl;
  
    nstep=2*nstep_final/temp;
  //nstep=nstep_final;
    
    ReadInput.close();
    
    outL.open(folder+"Best_Lenght.out");
  
  cout<<"...Starting simulation..."<<endl;
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
//CHECK:
void Check(Tour tour){
    
    bool check=tour.Check();
    if(!check){
      cerr <<"PROBLEM (Tour bonds) at temperature n."<<tcount-1<<", in the "<<attempted+1<< "step of Metropolis." <<endl<<"If n=-1, problem in initializing tour." <<endl; //già messaggio di errore in Tour
      exit(-1);
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

void PrintBestLenght2(int t){
  outL << t << " " << temp<< " " << best_lenght <<endl;
}

void PrintBest(void){
  ofstream outB(folder+"Best_Tour.out");
  //outB << tstep <<endl;
  
  for (int i=0; i<ncities; i++){
    int city_index=best_tour.GetCity(i);
    outB <<city_index << " " << positions[city_index].x <<" " <<positions[city_index].y <<endl;
  }
  outB.close();
}

void PrintConsole(){
if (tcount%10==0){
  cout <<"------------------------------------------------------------"<<endl;
  cout <<"  *  " <<tcount<< "-th temperature: " <<temp <<endl;
  cout <<"     n_"<<tcount<<" = "<<nstep<<endl;
  
  alpha=double(accepted)/double(attempted);
  e=(0.5-alpha)/0.5;
  cout <<"acceptance rate= " << alpha <<endl;
  cout <<"relative error= "<< abs(e) <<endl;
  }
}


//*************************************************************************************************//
//SIMULATED ANNEALING:
  void Reset(){
    beta=1./temp;
    tcount++;
    
    accepted = 0;
    attempted = 0;
  }

void FindBest(void){
  double lenght=x_tour.GetLenght2();
  if(lenght<best_lenght){
    best_tour=x_tour;
    best_lenght=lenght;
  }
}

void NextTemperature(void){
  temp*=t_rate;
  if (nstep<nstep_final) nstep+=100;
  else nstep=nstep_final;
}
//*************************************************************************************************//
//METROPOLIS:
void Metropolis(void){
  //Tour y_tour=x_tour;

  Tour y_tour=Mutate(x_tour); //---> y=T(x);
  Check(y_tour);
  
  //Probability to accept:
  double deltaL=y_tour.GetLenght2()-x_tour.GetLenght2();
  double q=exp(-beta*deltaL);
  double A=min(1.,q);
  
  //Accept?
  double r=rnd.Rannyu();
  if (r<=A){
    x_tour=y_tour;
    accepted++;
  }
  attempted++;
}

Tour Mutate(Tour a){
  Tour b=a;
  
  int type=rnd.Rannyu(0,3);
  
  if (type==0) b.RandomSwap(); //Swap 2 elements randomly choosen
  else if (type==1) b.Shift(); //Shift m elements of n positions
  else if (type==2) b.Inversion(); //Invert sub group of m elemnets

  return b;
}
//*************************************************************************************************//
void Clean(void){
  
  cout << endl<<"Done √"<<endl<<endl;
  
  alpha=double(accepted)/double(attempted);
  e=(0.5-alpha)/0.5;
  cout <<"Final acceptance rate= " << alpha <<endl;
  cout <<"Relative error= "<< abs(e) <<endl;
  
  cout << "After "<<tcount<<" temperatures, the final lenght is = "<<best_lenght<<endl<<endl;

  outL.close();
  rnd.SaveSeed();
}


//ho tolto i pointer
//ho disattivato uguale
//ho creato tour constructor
//potrei mettere il seed out anche in 9
//provare a vedere se cambia qualcosa se t->tcount in print best lenght
//togliere print seed
//devo pulire cartella e fare read.me
