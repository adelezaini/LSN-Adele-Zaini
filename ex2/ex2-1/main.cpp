/*******************************************************************************************************************************
Laboratorio Simulazione Numerica
Excercise 2.1

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
double f(double);
double p(double);
double g(double);

int main (int argc, char *argv[]){

//GENERATORE NUMERI CASUALI: lo inizializzo

   Random rnd;
    rnd.SetSeed();


    int M=10000;
    int N=100;
    int L=int(M/N);
    
    double m=0.,m2=0., err=0.;
    
    
// Parte a) Uniform sampling
    
//Note: esercizio analogo a ex1-1 (lo compatto un po'...)
ofstream out("uniform.out");
    if(!out.is_open()) cerr << "PROBLEM: Unable to open uniform.out" << endl;

    
    for(int i=0;i<N;i++){
        double Int=0., Int2=0.;
        for(int j=0;j<L;j++){
            double x=rnd.Rannyu();
            Int+=f(x);
        }
        Int/=double(L);
        Int2=Int*Int;
        
        m=(m*double(i)+Int)/double(i+1);
        m2=(m2*double(i)+Int2)/double(i+1);
        err=Error(m,m2,i);
        
        out << (i+1)*L << " " << m << " " << err << endl;
    }

    out.close();
    out.clear();
    
//Parte b) Importance sampling
    
out.open("importance.out");
    if(!out.is_open()) cerr << "PROBLEM: Unable to open importance.out" << endl;
    
    for(int i=0;i<N;i++){
        double Int=0., Int2=0.;
        for(int j=0;j<L;j++){
            double x=rnd.Rannyu();
            double y=1.-sqrt(1.-x); //Uso metodo della inversione della cumulativa e scelgo soluzione con meno per avere giusto intervallo.
            Int+=g(y);
        }
        Int/=double(L);
        Int2=Int*Int;
        
        m=(m*double(i)+Int)/double(i+1);
        m2=(m2*double(i)+Int2)/double(i+1);
        err=Error(m,m2,i);
        
        out << (i+1)*L << " " << m << " " << err << endl;
    }

    out.close();
    
//Finito!
    
   rnd.SaveSeed();
    
    
   return 0;
}

double Error(double ave, double ave2,int n){
    if(n==0) return 0;
    else return sqrt((ave2-ave*ave)/n);
}

double f(double x){
    return M_PI/2.*cos(M_PI*x/2.);
}
double p(double x){
    return 2.*(1.-x);
}
double g(double x){
    return M_PI/2.*cos(M_PI*x/2.)/(2.*(1.-x));
}




