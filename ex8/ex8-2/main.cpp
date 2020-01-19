/*******************************************************************************************************************************
Laboratorio Simulazione Numerica
Excercise 8: Variational Monte Carlo

 Adele Zaini
**********************************************************************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include "random.h"
#include "main.h"

//inserire nel jupyter calcolo di derivata seconda
//mettere read me con equilibration

//NOTE: ho preso main.cpp di ex5, cambiando:
//1. tutto da 3D a 1D (e.g. "Metropolis1D()")
//2. aggiunto funzioni Initialize() e Clean();
//3. cambiato ovviamente input.dat e quindi anche pdf to sample per Metropolis (unica wf2, no opzione di scelta)
//4. aggiunto Measure() invece che inglobarlo direttamente in Accumulate(), perchè logicamente due funzioni separate
//5. oss-> Measure() è già piu complicatina del solito
//6. aggiunto funzione Equilibrate() con ottimizzazione automatica di delta, magari simulazione meno efficiente e meno accurata, ma almeno tutto automatico visto che programma svolge più compiti rispetto a ex5.

using namespace std;

//*************************************************************************************************//
//MAIN

int main(){
  
  Input();
  Initialize();
  
  Equilibrate();

  for(int n=0; n<nblock; n++){
    Reset(n);
    for(int l=0; l<L; l++) {

      Metropolis1D();
      Measure();
      Accumulate();
    }
    Averages(n);
  }
  
  Clean();

  return 0;
}

//*************************************************************************************************//
//SET SIMULATION PARAMETERS:

void Input(void){
  
  ifstream ReadInput;

  cout << "-------------------------------------------------------------------------------" <<endl;
  cout << "          Ground state of a single quantum particle in 1D     " << endl;
  cout << "          Variational Monte Carlo simulation"<<endl;
  cout << "External potential V(x)=x^4-5/2*x^2" <<endl;
  cout << "All quantities are rescaled considering hbar=1 and m=1 " << endl;
  cout << "-------------------------------------------------------------------------------" <<endl<<endl;
  
  ReadInput.open("input.dat");

  ReadInput >> x_0;
  x=x_0;
  
  ReadInput >> delta; //mi serve comunque anche se poi lo cambio --> vedi Equilibrate per maggiori dettaglis
  ReadInput >> nblock;
  ReadInput >> L; //ATTENZIONE! Ora non è nstep tot, ma è nstep per blocco
  ReadInput >> tr;

  //L=nstep/nblock;
  
  //for equilibration:
  int expo=log10(x_0);
  nequi=pow(10,expo+1);

  cout << "The program uses the Metropolis algorithm to sample the probability density"<< endl;
  cout << "– Origin of the simulation = ("<<x_0<<")" <<endl;
  cout << "– Time step = " << delta << endl;
  cout << "– Number of blocks = " << nblock << endl;
  cout << "– Number of steps per block= " << L << endl;
  if (tr=="Uniform") {
    T=&Uniform;
    cout << "– Transition probability = Uniform distribution" <<endl<<endl;
  }else if(tr=="Gaussian"){
    T=&Gaussian;
    cout << "– Transition probability = Gaussian distribution" <<endl<<endl;
  }else{
    cerr << "PROBLEM: this program works only with uniform or gaussian distribution as transition probability"<<endl<<"Plese choose 'transition input' between 'Uniform' or 'Gaussian'."<< endl<<endl;
    exit(1);
  }
  //OSS: si possono aggiungerepdf anche qui
  
  ReadInput >> mu;
  ReadInput >> sigma;
  
  cout << "The trail wave function: wf_T=exp[-0.5*(x-μ)^2/σ^2]+exp[-0.5*(x+μ)^2/σ^2] depends on the following parameters:" <<endl;
  cout << "- μ = " <<mu<<endl;
  cout << "- σ = " <<sigma<<endl<<endl;
  
    ReadInput >> equi;
  
  if (equi) cout << "Automatic equilibration: ON" <<endl<<endl;
  else{cout << "Automatic equilibration: OFF." <<endl;
    cout <<"Manually equilibrate the time step looking at the final acceptance rate" <<endl<<endl;}
  
    ReadInput.close();
}

void Initialize(void){
  rnd.SetSeed();
  p=&wfgauss2;
  wfT=&wfgauss;
  wfT_d2=&wfgauss_d2;
  
  if (tr=="Uniform"){
    outE.open("EnergyT.0.out");
    outP.open("Positions.0.out");
  }else {
    outE.open("EnergyT_"+tr+".0.out");
    outP.open("Positions_"+tr+".0.out");
  }
}

//*************************************************************************************************//
void Equilibrate(void){
  //1. Approaching x_max:
  
  //OSS: dipende da delta, vero, ma meglio prima avvicinarsi e poi cambiare delta rispetto al contrario (scelta di delta maggiormente dipendente da posizione iniziale)
  for (int i=0; i<nequi; i++) Metropolis1D();
  
  //2. Optimizing delta:
  if(equi){
    double old_e; //old_delta, old_alpha;
    int count=0;
    bool min;
    
   do{
      //old_delta=delta;
     
     if (count!=0) old_e=e;
       //old_alpha=alpha;}
     
     accepted=0;
     attempted=0;
    for (int i=0; i<2000; i++) Metropolis1D();
     
     alpha=double(accepted)/double(attempted);
     e=(0.5-alpha)/0.5;
     
        //cout << "Accepted=" << accepted << "   Attempted= " << attempted << "    Alpha =" << alpha << "    e=" << e <<endl;
     
     if(count==0){
       //old_alpha=alpha;
       old_e=e;
       
       cout << "...Equilibration of the system: optimizing the time step 𝛿..." <<endl<<endl;
       //OSS: sfrutto che quando errore relativo su accettance rate è negativo, delta è più piccolo di quello ottimizzato:
       //ATT: andava bene per ex5, ma non sono sicura che qui uguale!
       if(e<0) min=true; // per capire se approccio valore accettazione da delta più grande e più piccolo, e capire dimensioni variazioni delta
       else min=false;
     }
     
    if (count==1000) break; //nel caso in cui si entri in un loop
      
      if(min){ //actualDelta ++++++ -> targetDelta (se lo supero ritorno un pochino indietro)
        if ((old_e*e)>=0.) delta+=0.003;
        else delta-=0.002;
      }else{   //targetDelta <- ------ actualDelta (se lo supero trono un pochino avanti)
        if ((old_e*e)>=0.) delta-=0.003;
        else delta+=0.002;
      }
     
         count++;
     
   } while(abs(e)>0.001 || count==1);
    
    /*delta=old_delta;
    alpha=old_alpha;*/
    cout << "Optimized time step 𝛿 = "<<delta<<", with trial acceptance rate = " << alpha<< " (reached after "<< count <<" iterations)." <<endl;

    //OSS: non mi preoccupo della posizione, anzi ora ancora più equilibrata
  }
}
  
  
//*************************************************************************************************//
void Metropolis1D(void){
  double y;
  double px,py, A, r;
  
  y=T(x);
  
  //Probability to accept:
  px=p(x);
  py=p(y);
  A=min(1.,py/px);
  
  //Accept?
  r=rnd.Rannyu();
  if (r<=A){
    x=y;
    accepted++;
  }
  attempted++;
  
  PrintPosition1D();
}

void PrintPosition1D(void){
  outP << x << endl;
}
//*************************************************************************************************//
void Measure(void){ //calculate H*wf/wf
  double Kx=-0.5*wfT_d2(x)/wfT(x);
  H=Kx+V(x);
}
//*************************************************************************************************//
//AVERAGE:
void Reset(int n){
    mb=0.;
    l_norm=0;
    if (n==0) {m=0.;m2=0.; err=0.;}
  accepted = 0;
  attempted = 0;
}

void Accumulate(void){
  mb+=H;
  l_norm++;
}

void Averages(int n){
    mb=mb/double(l_norm);
    m=(m*double(n)+mb)/double(n+1);
    m2=(m2*double(n)+mb*mb)/double(n+1);
    if (n!=0) err=sqrt((m2-m*m)/double(n));
  
  outE << (n+1) << " " << m << " " << err << endl;
}
//*************************************************************************************************//
void Clean(void){
  
  alpha=double(accepted)/double(attempted);
  e=(0.5-alpha)/0.5;
  cout << "Final acceptance rate= " << alpha <<endl;
  cout << "Relative error= "<< abs(e) <<endl<<endl;
  cout << "Final energy E_T= " << m <<endl<<endl;
  
  rnd.SaveSeed();
  outE.close();
  outP.close();
}

//*************************************************************************************************//

//Metropolis pdf:
//1.
double Uniform(double x){
  return x+rnd.Rannyu(-delta,delta);
}
//2.
double Gaussian(double x){
  return rnd.Gauss(x,delta);
}

//Trial wave function
double wfgauss(double x){
  double num1=-(x-mu)*(x-mu);
  double num2=-(x+mu)*(x+mu);
  double den=2.*sigma*sigma;
  return exp(num1/den)+exp(num2/den);
}

//Pdf to sample (square of wfT)
double wfgauss2(double x){
  return pow(wfgauss(x),2);
}

//Second derivative of wfT
double wfgauss_d2(double x){
  double m1=(x-mu)*(x-mu)/sigma/sigma;
  double m2=(x+mu)*(x+mu)/sigma/sigma;
  double expo1= exp(-0.5*m1);
  double expo2=exp(-0.5*m2);
  return (expo1*(m1-1.)+expo2*(m2-1.))/sigma/sigma;
}

//External potential
double V(double x){
  return pow(x,4)-5./2.*x*x;
}

double Error(double ave, double ave2,int n){
    if(n==0) return 0;
    else return sqrt((ave2-ave*ave)/n);
}
