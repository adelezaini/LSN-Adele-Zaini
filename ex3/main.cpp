/*******************************************************************************************************************************
Laboratorio Simulazione Numerica
Excercise 3: Econofisica

 Adele Zaini
**********************************************************************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include "random.h"

using namespace std;

double Error(double, double,int);


int main (int argc, char *argv[]){

//GENERATORE NUMERI CASUALI: lo inizializzo

   Random rnd;
    rnd.SetSeed();


    int M=100000;
    int N=100;
    int L=int(M/N);
    
    double const S0=100.;
    double const T=1.;
    double const K=100.;
    double const r=0.1;
    double const sigma=0.25;
    
    double  mC=0., mP=0., mC2=0., mP2=0., errC=0.,errP=0.;
  
// Parte a) Direct sampling
    
    
ofstream outC("direct-call.out");
    if(!outC.is_open()) cerr << "PROBLEM: Unable to open direct-call.out" << endl;
ofstream outP("direct-put.out");
    if(!outP.is_open()) cerr << "PROBLEM: Unable to open direct-put.out" << endl;
    
    for(int i=0;i<N;i++){
        double call=0., put=0.;
        for(int j=0;j<L;j++){
            double WT=rnd.Gauss(0.,T);
            double ST=S0*exp((r-0.5*sigma*sigma)*T+sigma*WT);
            call+=max(0.,ST-K)*exp(-r*T);
            put+=max(0.,K-ST)*exp(-r*T);
        }
        call/=double(L);
        put/=double(L);
        
        mC=(mC*double(i)+call)/double(i+1);
        mP=(mP*double(i)+put)/double(i+1);
        mC2=(mC2*double(i)+call*call)/double(i+1);
        mP2=(mP2*double(i)+put*put)/double(i+1);
        errC=Error(mC,mC2,i);
        errP=Error(mP,mP2,i);
        
        outC << (i+1)*L << " " << mC << " " << errC << endl;
        outP << (i+1)*L << " " << mP << " " << errP << endl;
    }

    outC.close();
    outC.clear();
    outP.close();
    outP.clear();
    
//Parte b) Discretized sampling
    int nsteps=100;
    
outC.open("discrete-call.out");
    if(!outC.is_open()) cerr << "PROBLEM: Unable to open discrete-call.out" << endl;
outP.open("discrete-put.out");
    if(!outP.is_open()) cerr << "PROBLEM: Unable to open discrete-put.out" << endl;
    
    for(int i=0;i<N;i++){
        double call=0., put=0.;
        for(int j=0;j<L;j++){
            double ST=S0;
            double dT=T/double(nsteps);
            for(int k=0; k<nsteps; k++){
                double Z=rnd.Gauss(0.,1.);
                ST=ST*exp((r-0.5*sigma*sigma)*dT+sigma*Z*sqrt(dT));
            }
            call+=max(0.,ST-K)*exp(-r*T);
            put+=max(0.,K-ST)*exp(-r*T);
        }
        call/=double(L);
        put/=double(L);
        
        mC=(mC*double(i)+call)/double(i+1);
        mP=(mP*double(i)+put)/double(i+1);
        mC2=(mC2*double(i)+call*call)/double(i+1);
        mP2=(mP2*double(i)+put*put)/double(i+1);
        errC=Error(mC,mC2,i);
        errP=Error(mP,mP2,i);
        
        outC << (i+1)*L << " " << mC << " " << errC << endl;
        outP << (i+1)*L << " " << mP << " " << errP << endl;
    }

    outC.close();
    outP.close();
    
//Finito!
    
   rnd.SaveSeed();
    
    
   return 0;
}

double Error(double ave, double ave2,int n){ //all'errore passo già "n-1"
    if(n==0) return 0;
    else return sqrt((ave2-ave*ave)/n);
}





