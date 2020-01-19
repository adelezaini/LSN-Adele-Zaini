/*******************************************************************************************************************************
Laboratorio Simulazione Numerica
Excercise 1.2
 
 Adele Zaini
**********************************************************************************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <iomanip>
#include "random.h"

using namespace std;


int main (int argc, char *argv[]){

//GENERATORE NUMERI CASUALI: lo inizializzo

   Random rnd;
	rnd.SetSeed();


//Creo cicli e sottocicli per avere 3x4 dataset di 10^4 elementi
    int n=10000;
    int N[4]={1,2,10,100};
    
    
    ofstream out1("dice.out");
    if(!out1.is_open()) cerr << "PROBLEM: Unable to open dice.out" << endl;
     ofstream out2("exp.out");
    if(!out2.is_open()) cerr << "PROBLEM: Unable to open exp.out" << endl;
     ofstream out3("lor.out");
    if(!out3.is_open()) cerr << "PROBLEM: Unable to open lor.out" << endl;
    
    // Creo ciclo generale per i vari N
    for (int k=0; k<4; k++){
        
    // Creo ciclo per il numero delle somme che voglio calcolare per ogni N fissato
        for(int i=0;i<n;i++){
            double dice=0.;
            double exp=0.;
            double lor=0.;
            
    // Creo ciclo per somme cumulative Sn, come es prec
            for(int j=0;j<N[k];j++){
                dice+=int(rnd.Rannyu(1.,7.));
                exp+=rnd.Exponential(1.);
                lor+=rnd.Lorentz(0.,1.);
            }
            dice/=double(N[k]);
            exp/=double(N[k]);
            lor/=double(N[k]);
            
            out1 << dice << endl;
            out2 << exp << endl;
            out3 << lor << endl;
        }
    }
    
    out1.close();
    out2.close();
    out3.close();
    

//Finito!

//Ricordati: sui file hai
    //(0,10000): N=1;
    //(10000,20000): N=2;
    //(20000,30000): N=10;
    //(30000,40000): N=100;
//--> impara come estrappolarli a blocchi con Python
    
   rnd.SaveSeed();
    
    
   return 0;
}






