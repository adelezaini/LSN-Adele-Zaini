/*******************************************************************************************************************************
Laboratorio Simulazione Numerica
Excercise 5: Metropolis

 Adele Zaini
**********************************************************************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include "random.h"
#include "main.h"

using namespace std;

//OSS: struttura analoga a ex4, ma diverso comando di output che non è "append" ma "overwrite" perchè inglobando fase di equilibrazione dentro al programma non ho bisogno di salvare valori di simulazioni precedenti, e anzi così non ho bisogno di mandare "bash clean.sh" ogni volta ma tutto automatico

//*************************************************************************************************//
//MAIN

int main(){
  
  Input();
  rnd.SetSeed();
  
  //Equilibration:
  for (int i=0; i<nequi; i++) Metropolis();
  
  for(int n=0; n<nblock; n++){
    Reset(n);
    for(int l=0; l<L; l++) {

      Metropolis();
      Accumulate(); //and Measure()
      
      /*int istep=n*L+l+1;
       if(istep%(nstep/10.) == 0) cout << "Number of time-steps: " << istep << endl;*/
    }
    Averages(n);
  }
  
  //Conclude:
  cout << "Acceptance rate= " << double(acc)/double(nstep) <<endl;
  cout << "Relative error= "<< (0.5-double(acc)/double(nstep))/0.5 <<endl<<endl;
  
  rnd.SaveSeed();
  outR.close();
  outP.close();


  return 0;
}

//*************************************************************************************************//
//SET SIMULATION PARAMETERS:

void Input(void){
  
  ifstream ReadInput;

  cout << "-------------------------------------------------------------------------------" <<endl;
  cout << "          Hydrogen atom       " << endl;
  cout << "Simulation of spatial distribution of electrons for 1s and 2p orbitals" << endl << endl;
  cout << "The program is based on the Metropolis algorithm to sample wave functions"<< endl;
  cout << "All quantities are rescaled in Bohr radius units " << endl;
  cout << "-------------------------------------------------------------------------------" <<endl<<endl;
  
  ReadInput.open("input.dat");

  ReadInput >> orbital; //string
  cout << "– The orbital chosen to be sampled is " << orbital << endl;
  
  //assegno la probabilità da campionare in base all'orbitale scelto:
  //OSS: possono essere implementate altre opzioni per altri orbitali (o in gen altre distribuzioni)
  if (orbital=="1s") {
    p=&wf100;
    cout << "The corresponding pdf is the square of the Hydrogen wave function with quantic numbers (n=1,l=0,m=0)" <<endl<<endl;
  }else if(orbital=="2p") {
    p=&wf210;
    cout << "The corresponding pdf is the square of the Hydrogen wave function with quantic numbers (n=2,l=1,m=0)" <<endl<<endl;
  }else {
    cerr << "PROBLEM: this program evaluates only 1s or 2p orbitals." <<endl<< "Please enter 'orbital input' again" <<endl<<endl;
    exit(1);
  }
  
  ReadInput >> x_0; x[0]=x_0;
  ReadInput >> y_0; x[1]=y_0;
  ReadInput >> z_0; x[2]=z_0;
  
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> nblock;
  ReadInput >> tr;

  L=nstep/nblock;
  int expo=max({log10(x_0),log10(y_0),log10(z_0)}); //è già un integer
  nequi=pow(10,expo+1);

  cout << "– Origin of the simulation = ("<<x_0<<", "<<y_0<<", "<<z_0<<")" <<endl;
  cout << "– Time step = " << delta << endl;
  cout << "– Number of steps = " << nstep << endl;
  cout << "– Number of blocks = " << nblock << endl;
  if (tr=="Uniform") {
    T=&Uniform;
    cout << "– Transition probability = Uniform distribution" <<endl<<endl;
  }else if(tr=="Gaussian"){
    T=&Gaussian;
    cout << "– Transition probability = Gaussian distribution" <<endl<<endl;
  }else{
    cerr << "PROBLEM: this program works only with uniform or gaussian distribution as transition probability"<<endl<<"Plese choose 'transition input' between 'uniform' or 'gaussian'."<< endl<<endl;
    exit(1);
  }
  //OSS: si possono aggiungerepdf anche qui
  
    ReadInput.close();
  
  outR.open(tr+"/Radius_"+orbital+".out");
  outP.open(tr+"/Positions_"+orbital+".out");
}

//*************************************************************************************************//
void Metropolis(void){ //evaluate next step
  double y[3];
  double px,py, A, r;
  
  for (int i=0; i<3; i++) y[i]=T(x[i]);
  
  //Probability to accept:
  px=p(x[0],x[1],x[2]);
  py=p(y[0],y[1],y[2]);
  A=min(1.,py/px);
  
  //Accept?
  r=rnd.Rannyu();
  if (r<=A){
  for(int i=0; i<3; i++) x[i]=y[i];
    acc++;
  }
  
  PrintPosition(); //la tengo separata per eventualmente disattivarla se necessario
}

void PrintPosition(void){
  outP << x[0] << " " << x[1] << " " << x[2] << endl;
}

//*************************************************************************************************//
//AVERAGE:
void Reset(int n){
    mb=0.;
    l_norm=0;
    if (n==0) {m=0.;m2=0.; err=0.;}
}

void Accumulate(void){ //and Measure() --> calcola raggio e allo stesso tempo accumula il valore per la media del blocco
  R=sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
  mb+=R;
  l_norm++;
}

void Averages(int n){
  
    mb=mb/double(l_norm);
    m=(m*double(n)+mb)/double(n+1);
    m2=(m2*double(n)+mb*mb)/double(n+1);
    if (n!=0) err=sqrt((m2-m*m)/double(n));
  
  //PrintRadius():
  outR << (n+1) << " " << m << " " << err << endl;
}
//*************************************************************************************************//

double Error(double ave, double ave2,int n){
    if(n==0) return 0;
    else return sqrt((ave2-ave*ave)/n);
}

double Uniform(double x){
  return x+rnd.Rannyu(-delta,delta);
}
double Gaussian(double x){
  return rnd.Gauss(x,delta);
}
double wf100(double x, double y, double z){ //ricordati che è funzione al quadrato!!
  double r=sqrt(x*x+y*y+z*z);
  return exp(-2.*r)/M_PI;
}

double wf210(double x, double y, double z){
  double r=sqrt(x*x+y*y+z*z);
  double cos_theta=z/r;
  return r*r*exp(-r)*cos_theta*cos_theta/32./M_PI; // se no torna qualcosa riscrivi bello questo
}


/*  do{
    //genero nuova posizione candidata:
    for (int i=0; i<3; i++) y[i]=T(x[i]);
    
    //Probability to accept:
    px=p(x[0],x[1],x[2]);
    py=p(y[0],y[1],y[2]);
    A=min(1.,py/px);
    
    //Accept?
    r=rnd.Rannyu();
    tot++;
  }while(r>A);
    
  for(int i=0; i<3; i++) x[i]=y[i];
 
 cout << "Acceptance rate=" << double(nstep)/double(tot) << endl<<endl;
 cout << 0.5-double(nstep)/double(tot) <<endl<<endl;
  
}*/

