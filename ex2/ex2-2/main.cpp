/*******************************************************************************************************************************
Laboratorio Simulazione Numerica
Excercise 2.2: Random walk

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
double Distance(double[]);

int main (int argc, char *argv[]){

//GENERATORE NUMERI CASUALI: lo inizializzo

   Random rnd;
    rnd.SetSeed();


    int M=10000;
    int N=100;
    int L=int(M/N);
    
    int S=100;
    double pos[3]; //mi serve per immagazzinare la posizione prec
    double a=1.;
    
    double RMS[N][S]; // per il calcolo finale mi serve calcolare 1. la media sugli N blocchi 2. per ogni step s-esimo
    
    double m=0.,m2=0., err=0.;
    
    
// Parte a) Lattice
    
ofstream out("discrete.out");
    if(!out.is_open()) cerr << "PROBLEM: Unable to open discrete.out" << endl;
    
    for(int n=0;n<N;n++){
        for(int s=0; s<S; s++) RMS[n][s]=0.;
    }
    
    for(int n=0;n<N;n++){

        for(int j=0;j<L;j++){
    
            for (int dim=0; dim<3; dim++) pos[dim]=0;
            
            for(int s=0; s<S; s++){
                //step discreti:
                int dir=rnd.Rannyu(0.,3.); //scelgo x,y, o z
                int v=rnd.Rannyu(0., 2.); //scelgo avanti o indietro
                if (v==0) pos[dir]+=a;
                else pos[dir]-=a;
                
                double r=Distance(pos); //calcolo la distanza dall'origine dopo il nuovo passo
                RMS[n][s]+=r*r; //carico il vettore per le distanze quadre di ogni step fino a...
            //quando arrivo all'ultimo elemento del blocco n-esimo:
                if (j==(L-1)) RMS[n][s]/=double(L); //le faccio diventare le medie per ogni blocco
            }
        }
    }
    //A questo punto ho matrice caricata {r2_[n,s]} con r2 mediati su ogni blocco
    //Calcolo RMS per ogni step e stampo
    
    for(int s=0; s<S; s++){
      m=0; m2=0.; err=0.;
        for(int n=0;n<N;n++){
          m+=RMS[n][s];
          m2+=RMS[n][s]*RMS[n][s];
        }
        m/=double(N);
        m2/=double(N);

        err=Error(m,m2,N);

        out << (s+1) << " " << sqrt(m) << " " << err << endl;
    }
    
    out.close();
    out.clear();
    
// Parte b) Continuum
    
    out.open("continuous.out");
    if(!out.is_open()) cerr << "PROBLEM: Unable to open continuous.out" << endl;
    
    for(int n=0;n<N;n++){
        for(int s=0; s<S; s++) RMS[n][s]=0.;
    }
    
    for(int n=0;n<N;n++){
        for(int j=0;j<L;j++){
            for (int dim=0; dim<3; dim++) pos[dim]=0;
            
            for(int s=0; s<S; s++){
                
                //step continui:
                double theta=rnd.Theta3D();
                double phi=rnd.Angle();
                pos[0]+=a*sin(theta)*cos(phi);
                pos[1]+=a*sin(theta)*sin(phi);
                pos[2]+=a*cos(theta);
                
                double r=Distance(pos); //calcolo la distanza dall'origine dopo il nuovo passo
                RMS[n][s]+=r*r; //carico il vettore per le distanze quadre di ogni step fino a...
            //quando arrivo all'ultimo elemento del blocco n-esimo:
                if (j==(L-1)) RMS[n][s]/=double(L); //le faccio diventare le medie per ogni blocco
            }
        }
    }
    //A questo punto ho matrice caricata {r2_[n,s]} con r2 mediati su ogni blocco
    //Calcolo RMS per ogni step e stampo
    
    for(int s=0; s<S; s++){
        m=0; m2=0.; err=0.;
        for(int n=0;n<N;n++){
          m+=RMS[n][s];
          m2+=RMS[n][s]*RMS[n][s];
        }
        m/=double(N);
        m2/=double(N);

        err=Error(m,m2,N);

        out << (s+1) << " " << sqrt(m) << " " << err << endl;
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

double Distance(double p[]){
    return sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
}


