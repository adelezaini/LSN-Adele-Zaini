/*******************************************************************************************************************************
Laboratorio Simulazione Numerica
Excercise 1.3: Buffon's experiment

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
	int l=int(M/N);
    
    double d=1.;
    double L=0.8;
    int Nhit=0;
    
    double pi=0., pi2=0., m=0.,m2=0., err=0.;
    
ofstream out("pi.out");
    if(!out.is_open()) cerr << "PROBLEM: Unable to open pi.out" << endl;

    
	//Rispetto al solito compatto il codice:
    
	for(int i=0;i<N;i++){
        Nhit=0;
        
        //Ciclo per l'esperimento:
        for(int j=0;j<l;j++){
            //a) Lancio l'ago:
            double x=rnd.Rannyu(0.,d);
            //double theta=rnd.Angle(); //--> Nota: avevo implementato il metodo, ma per le condizioni che qui imposto (cos>0) posso considerare solo alcuni angoli, quindi preferisco tenere il metodo per casi generali e ora utilizzarne uno ad hoc
                double xt,yt, r;
                do{
                xt = rnd.Rannyu();
                yt = rnd.Rannyu();
                r = sqrt(xt*xt + yt*yt);
                }while(r>1.);
    
            //b) Interseca la linea?
            if ((x-L/2.*xt/r)<=0. || (x+L/2.*xt/r)>=d) Nhit++;

        }
        //cout << Nhit << endl;
		pi=2*L*double(l)/double(Nhit)/d;
		pi2=pi*pi;
        
        m=(m*double(i)+pi)/double(i+1);
        m2=(m2*double(i)+pi2)/double(i+1);
        err=Error(m,m2,i);
        
        out << (i+1)*l << " " << m << " " << err << endl;
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





//- in order to random generate the $\theta$ angle, I implemented a new method in the PRNG based on the inverse function: $\theta=arctan(y/x)$. Even though it is defined in $[-\pi/2,\pi/2]$, this is not a problem because of the problem simmetry.

