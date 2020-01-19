using namespace std;

//random numbers
#include "random.h"
Random rnd;

// simulation
string orbital;
string tr;
int nstep, iprint, nblock, L;

//Metropolis parameters
double delta;
double (*p)(double, double, double);
double (*T)(double);
int nequi, acc=0;

// averages
int l_norm;
double m=0.,m2=0.,err=0.,mb=0.;

//observables
double x_0, y_0, z_0, R;
double x[3];

//output files
ofstream outR, outP;


//functions
void Input(void);

//void Equilibrate(void);
void Metropolis(void);
void PrintPosition(void);

void Reset(int);
void Accumulate(void);
void Averages(int);

double wf100(double, double, double);
double wf210(double, double, double);
double Uniform(double);
double Gaussian(double);

double Error(double, double,int);
