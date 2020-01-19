/*******************************************************************************************************************************
Laboratorio Simulazione Numerica
Excercise 9: Genetic Algorithm - Traveling Salesman Problem

 Adele Zaini
**********************************************************************************************************************************/

#include "TSP.h"

using namespace std;

/*------------------------------------------------------------------------//
                             POPULATION
//------------------------------------------------------------------------*/
//*********************************************************//
//CREATOR
Population::Population(int dim, vector<point> pos){
  
  npop=dim;
  ncities=pos.size();
  
  for (int n=0; n<npop; n++){
    Tour tour(pos);
    tours.push_back(tour); //I use creator from positions -> random sorting of city index
  }
  
  iter=0;
  
  rnd=new Random();
  rnd->SetSeed("seed.out");
    rnd->SaveSeed();
  
  fitness_sum=0.;
}

//*********************************************************//
//NEW GENERATION:

void Population::Offspring(double Pc){
  
  vector<Tour> children;
  iter++;
  
//SORT on fitness
  Sort();//it does not change probability if then the vector is sorted (r is random), but important for elite

  
//ELITE:
  int nelite;
  if(iter>50){ //attivo elite solo dopo 50 iterazioni per essere sicura non finire in un minimo locale
    nelite=npop*0.05+1; //I choose ~5% of population; //if npop<10, there is at least one
    for (int c=0; c<nelite; c++) children.push_back(tours[c]);
  }else{
    nelite=0;
  }
  
//FITNESS TOT
  double sum=0.;
  for (int i=nelite; i<npop; i++) sum+=tours[i].Fitness();
  fitness_sum=sum;

//SELECT new population with possibility of crossover
  for (int c=nelite; c<npop; c++){
    
    if(rnd->Rannyu()<Pc){
      
      int i_mom=Select();
      int i_dad=Select();
      
      while(i_dad==i_mom) {i_dad=Select();}
      Tour child=Crossover(i_mom,i_dad);
      children.push_back(child);
      
    }else{
      int i_good=Select();
      children.push_back(tours[i_good]);
    }
  }
  if(tours.size()==children.size()) tours=children;
  else{
    cerr <<"PROBLEM: something went wrong during creating new population n. "<<iter<<"."<<endl<<endl;
    exit(-1);
  }
  
}
//*********************************************************//
//SELECTION
int Population::Select()const{//Roulette russa   //potrei modificarla!
    double r=rnd->Rannyu(0.,fitness_sum);
    double isum=0.;
    int j=0;
    do{
      isum+=tours[j].Fitness();
      j++;
    }while(isum<r and j<npop);
  return j-1;
}

//*********************************************************//
//CROSSOVER:
Tour Population::Crossover(int m, int d)const{
  
  int cut=rnd->Rannyu(1, ncities-1);
  vector<int> son;
  vector<int> part2;
  
  //Ricopio parte di mamma
  for(int i=0; i<cut;i++) son.push_back(tours[m].GetCity(i));
  
  //Scrivo rimanenti nell'ordine del papa
  for(int k=0; k<ncities; k++){
    int dad_appo=tours[d].GetCity(k);
    bool diff=true;
    for(int j=0; j<cut; j++){
      int mom_appo=tours[m].GetCity(j);
      if(mom_appo==dad_appo) diff=false;
    }
    if (diff!=0) son.push_back(dad_appo);
  }
  
    vector<point> pos=tours[m].GetPositions();
  if(pos.size()==son.size()){
    Tour child(pos, son); //deve essere oggetto di tipo tour --> deve essere creato con le stess  posizioni degli altri
    return child;
  }else{
    cerr << "PROBLEM: something went wrong during crossover of population n. "<<iter<<"."<<endl<<endl;
    exit(-1);//exit the program
  }
}

//*********************************************************//
//MUTATION:
void Population::Mutation(double Pm){
  
  for (int i=0; i<npop; i++){
    
    if(rnd->Rannyu()<Pm){
      
      int type=rnd->Rannyu(0,4);
      
      if (type==0) tours[i].RandomSwap();
      else if (type==1) tours[i].Shift();
      else if (type==2) tours[i].ShiftTot();
      else tours[i].InversionTot();
      
    }
  }
}
//*********************************************************//
//BEST
Tour Population::Best()const{
  double best_fitness=0.;
  int best_index;
  
  for (int i=0;i<npop;i++){
    double fitness=tours[i].Fitness();
    if(fitness>best_fitness){
      best_index=i;
      best_fitness=fitness;
    }
  }
  return tours[best_index];
}

double Population::BestFitness()const{
  Tour best=Best();
  return best.Fitness();
}
double Population::BestLenght2()const{
  Tour best=Best();
  return best.GetLenght2();
}

double Population::BestHalfLenght2(){
  Sort();
  double m=0.;
  int stop=npop*0.5;
  for (int i=0; i<stop; i++) m+=tours[i].GetLenght2();
  m/=double(stop);
  return m;
}
//*********************************************************//
//SORT:
void Population::Sort(){
  quickSort(0,npop-1);
}

void Population::quickSort(int first, int last){ //quicksort
  if(first>last){
    int temp=last;
    last=first;
    first=temp;
  }
  
  if (last-first<=1){
    if(tours[first].GetLenght2()>tours[last].GetLenght2()) Swap(first, last);
    return;
  }else{
    int middle=(first+last)*0.5;
    double pivot=tours[middle].GetLenght2();
    int low=first;
    int high=last;
    double fit_low, fit_high;
    
    while(low<high){
      do{
        fit_low=tours[low].GetLenght2();
        low++;
      }while(fit_low<pivot);
      do{
        fit_high=tours[high].GetLenght2();
        high--;
      }while(fit_high>pivot);
      
      low--; //sono andata un po' oltre
      high++;
      
      if(low<high){
        Swap(low, high);
        low++;
      }
    }
    quickSort(first,low-1);
    quickSort(low,last);
  }
}

void Population::Swap(int i, int j){
  swap(tours[i], tours[j]);
}

