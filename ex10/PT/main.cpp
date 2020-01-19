/*******************************************************************************************************************************
Laboratorio Simulazione Numerica
Excercise 10: Parallel Tempering – Traveling Salesman Problem

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
#include "mpi.h"

using namespace std;


//mpi
int Size, Rank;
Tour bestest_tour;
double bestest_lenght;
int bestest_rank;

void InitializeRank();
void FindBestest();
void PrintBestest();

//*************************************************************************************************//
//MAIN

int main (int argc, char **arg){ //*argv[]
  
  /*if (argc<2){ //controlla se giusto
    cerr << "Usage: " << argv[0] << " <n_ranks> ..."<<endl; //<output_file_name>" << endl;
      return -1;
  }else size=atoi(argv[1]);*/
  
  Input();
  
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &Size);
  MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
 
  InitializeRank();

  //Same for each rank:-------------------------------------------------------
  for(int t=0; t<tstep; t++){
    Reset();
    for(int l=0; l<nstep; l++) {
      Metropolis();
    }
    //PrintConsole();
    
    FindBest();
    PrintBestLenght2(t);
    NextTemperature();
  }
  PrintBest(); //I still want to find and print every single best tour
  Clean();
  //-------------------------------------------------------------------------
  
  //Compare the different process:
  FindBestest();
  //PrintBestest();
  
  MPI_Finalize();

  return 0;
}

//*************************************************************************************************//
//SET SIMULATION PARAMETERS:

void Input(void){
  
  
/*//******* Initialize vector of bestssss: *****************************************************************************
  //best_of_rank is a vector of structures (argv[1]=size)
  for(int i=0; i<size; i++){
    Tour tour_appo;
    best_of_rank.push_back(best());
    best_of_rank[i].tour=tour_appo;
    best_of_rank[i].lenght=0.;
  }
  //NB: I dont insert "push_back" directly during MPI_processes, but I split it in two steps:
  //1. create vector of right dimension
  //2. load it with initial tour in Input() –> best_tours[rank]=x_tour;
  //this is why, otherwise the one-to-one correspondance index<->rank could not be preserved (not sequential work, but parallel indeed...);*/
    
  
//******* Initialize random variable for configurations (same for all ranks): ******************************************
  rnd.SetSeed();
  rnd.SaveSeed();
  //NB: all the ranks refers to the same configuration
  
  
//******** Input *******************************************************************************************************
    ifstream ReadInput;
    
    cout << "------------------------------------------------------------" <<endl;
    cout << "          Travelling Salesman Problem     " << endl;
    cout << "              Parallel Tempering of Simulated Annealing simulations                 "<<endl;
    cout << "-------------------------------------------------------------------------" <<endl<<endl;
  
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
    
    //Set Simulated Annealing parameters:
    ReadInput >> temp;
    ReadInput >> t_rate;
    ReadInput >> tstep;
    
    ReadInput >> nstep_final;
    
    cout <<endl<< "Simulated Annealing parameters:"<<endl;
    cout << "– Initial temperature T = " << temp << endl;
    cout << "- Reduction rate of temperature α=" << t_rate <<endl;
    cout << "– Dimension of annealing schedule N = " << tstep << endl;
    
    cout << "At a given T, the program uses the Metropolis algorithm to sample configurations according to the Boltzmann distribution." << endl;
    cout << "– Number of Metropolis steps n = " << nstep_final << endl<<endl;
  
    nstep=2*nstep_final/temp;
  //nstep=nstep_final;
    
    ReadInput.close();

}
//*************************************************************************************************//
void InitializeRank(){
  
  if (Rank==0) cout << "Parallel Tempering:"<<endl<<"The programs works on "<<Size<< " nodes"<<endl<<endl;
  
  cout<<"...Starting simulation of "<<Rank+1<<"-th node..."<<endl;
  
//Create initial configuration with city randomly sorted and Random variable re-initialized (different for each rank):
  rnd.SetSeed(Rank);
  Tour tour(positions, &rnd);
  x_tour=tour;
  Check(x_tour);

  //li tengo di supporto per ogni processo:
  best_tour=x_tour;
  best_lenght=x_tour.GetLenght2();
  
  /*//per confronto e stampa finale:
  best_of_rank[rank].tour=x_tour;
  best_of_rank[rank].lenght=x_tour.GetLenght2();*/
  
  cout <<endl<< Rank+1<<") Tour initialized √ -> Initial lenght = "<<best_lenght<<endl<<endl;
  
  outL.open(folder+"Best_Lenght_"+to_string(Rank)+".out");
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
//SIMULATED ANNEALING and PARALLEL TEMPERING:
  void Reset(){
    beta=1./temp;
    tcount++;
    
    accepted = 0;
    attempted = 0;
  }

void FindBest(void){
  double xlenght=x_tour.GetLenght2();
      
  if(xlenght<best_lenght){
    
    best_tour=x_tour;
    best_lenght=xlenght;
    
    //best_of_rank[rank].tour=x_tour;
    //best_of_rank[rank].lenght=xlenght;
  }
}

void FindBestest(){
  struct lenght_rank{ //creo struttura per far lavorare in automatico MINLOC
      double val;
      int rank;
  }bestest_lenght_, best_lenght_;
  
  
  best_lenght_.val=best_lenght;
  best_lenght_.rank=Rank;
  
  //double bestest_lenght_rank[2]={0.}; //it is a vector: [0]=bestestlenght; [1]=rank; //vedi se puoi metterlo nel main
  //double best_lenght_rank[2]={best_lenght, double(Rank)};
  
  MPI_Reduce(&best_lenght_, &bestest_lenght_, 1, MPI_DOUBLE_INT,MPI_MINLOC,0,MPI_COMM_WORLD);
  
  if(Rank==0){
    bestest_lenght=bestest_lenght_.val;
    bestest_rank=bestest_lenght_.rank;
    
    cout << "The bestest tour is evaluated by the "<<bestest_rank+1<<"-th node, with the lenght: L2="<<bestest_lenght<<endl;
  }
  if (Rank==bestest_rank) PrintBestest(); //bestest_tour=best_tour;
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
  ofstream outB(folder+"Best_Tour_"+to_string(Rank)+".out");
  //outB << tstep <<endl;
  
  for (int i=0; i<ncities; i++){
    int city_index=best_tour.GetCity(i);
    outB <<city_index << " " << positions[city_index].x <<" " <<positions[city_index].y <<endl;
  }
  outB.close();
}

void PrintBestest(void){
  ofstream outB(folder+"Bestest_Tour.out");
  outB << bestest_rank <<endl;
  
  for (int i=0; i<ncities; i++){
    int city_index=best_tour.GetCity(i);
    outB <<city_index << " " << positions[city_index].x <<" " <<positions[city_index].y <<endl;
  }
  outB.close();
}

/*void PrintBestest(void){
  ofstream outB(folder+"Bestest_Tour.out");
  outB << bestest_rank <<endl;
  
  for (int i=0; i<ncities; i++){
    int city_index=bestest_tour.GetCity(i);
    outB <<city_index << " " << positions[city_index].x <<" " <<positions[city_index].y <<endl;
  }
  outB.close();
}*/

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
void Clean(void){
  
  cout << endl<<Rank<<") Done √"<<endl<<endl;
  
  /*alpha=double(accepted)/double(attempted);
  e=(0.5-alpha)/0.5;
  cout <<"Final acceptance rate= " << alpha <<endl;
  cout <<"Relative error= "<< abs(e) <<endl;*/
  
  cout << Rank<<") After "<<tcount<<" temperatures, the final lenght is = "<<best_lenght<<endl<<endl;

  outL.close();
  rnd.SaveSeed();
}

//#include "/Users/antonellatagliavini/mpi/include/mpi.h"
