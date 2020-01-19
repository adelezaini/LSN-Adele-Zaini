/*******************************************************************************************************************************
Laboratorio Simulazione Numerica
Excercise 1.1

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

// Parte a) Calcolare media di una funzione continua tra [0,1)

    
	int M=10000; //Lanci totali
	int N=100; //N. blocchi
	int L=int(M/N); // N. lanci in ogni blocco


	double *m=new double[N]; //Array per medie
	double *m2=new double[N]; //Array per medie^2
    
    double n,s,s2,err; //Variabili per la statistica

ofstream out1("mean.out");
    if(!out1.is_open()) cerr << "PROBLEM: Unable to open mean.out" << endl;
	
	//Carico m[N] con medie sul singolo blocco:
	for(int i=0;i<N;i++){
		s=0.;
        for(int j=0;j<L;j++){
            n=rnd.Rannyu();
            s+=n;
        }
		m[i]=s/double(L);
		m2[i]=m[i]*m[i];
	}

	//Tendenza a convergere: calcolo somme cumulative e errori e stampo i risultati sul file output
	for(int i=0;i<N;i++){
        s=0.; s2=0.; err=0.;
		for(int j=0;j<(i+1);j++){
			s+=m[j];
			s2+=m2[j];
		}
		s/=double(i+1);
		s2/=double(i+1);
		err=Error(s,s2,i);
        
        out1 << (i+1)*L << " " << s << " " << err << endl;
	}

    out1.close();
    
    
// Parte b) Calcolare le varianze e i relativi errori
    
    M=100000;
    N=100;
    L=int(M/N);
    
ofstream out2("variance.out");
    if(!out2.is_open()) cerr << "PROBLEM: Unable to open variance.out" << endl;
    
    //Ricominciamo:
    for (int i=0;i<N;i++){
        m[i]=0.;
        m2[i]=0.;
    }

    
    //Tutto come prima ma cambiando formulina della funzione
    for(int i=0;i<N;i++){
        s=0.;
        for(int j=0;j<L;j++){
            n=rnd.Rannyu();
            s+=(n-0.5)*(n-0.5);
        }
        m[i]=s/double(L);
        m2[i]=m[i]*m[i];
    }

    for(int i=0;i<N;i++){
        s=0.; s2=0.; err=0.;
        
        for(int j=0;j<(i+1);j++){
            s+=m[j];
            s2+=m2[j];
        }
        
        s/=double(i+1);
        s2/=double(i+1);
        err=Error(s,s2,i);
        
        out2 << (i+1)*L << " " << s << " " << err << endl;
    }

    out2.close();
    
    
// Parte c) Test chi quadro
    
    M=100;
    n=10000;
    double ni=n/M;
    int* count=new int[100]; //array per contare throws in ciascun sub intervallo
    double somma=0.;
    
    ofstream out3("chi.out");
    if(!out3.is_open()) cerr << "PROBLEM: Unable to open chi.out" << endl;
    
    for (int k=0; k<100; k++){ //eseguo il test 100 volte
        
        for (int i=0;i<M;i++) count[i]=0;
        
        //1. Riempio il vettore count:
        for(int i=0; i<n; i++){
            double c=rnd.Rannyu();
            int j=c*M; //indice del subintervallo in cui finisce c
            count[j]++;
        }
        
        //2. Calcolo il chi quadro:
        double chi=0.;
        for (int i=0; i<M; i++){
            double xi=double(count[i]);
            chi+=(xi-ni)*(xi-ni)/ni;
        }
        somma+=chi;
        out3 << k+1 <<" "<< chi << endl;
    }
    //cout << somma/100 << endl;
    
    out3.close();

//Finito!
    
   rnd.SaveSeed();
    
    delete [] m;
    delete [] m2;
    delete [] count;
    
   return 0;
}

double Error(double ave, double ave2,int n){
	if(n==0) return 0;
	else return sqrt((ave2-ave*ave)/n);
}







