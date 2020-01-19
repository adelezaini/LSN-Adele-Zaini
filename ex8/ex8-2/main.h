using namespace std;

//random numbers
#include "random.h"
Random rnd;

// simulation
int L, nblock;
bool equi;

//Metropolis parameters
double delta;
double (*p)(double);
double (*T)(double);
string tr;
int nequi;
double alpha,e;
int accepted=0;
int attempted=0;

// averages
int l_norm;
double m=0.,m2=0.,err=0.,mb=0.;

//observables and parameters
double x_0, x, H, mu, sigma;
double (*wfT)(double);
double (*wfT_d2)(double);


//output files
ofstream outP, outE;


//functions
void Input(void);
void Initialize(void);
void Equilibrate(void);

void Metropolis1D(void);
void PrintPosition1D(void);
void Measure(void);

void Reset(int);
void Accumulate(void);
void Averages(int);

void Clean(void);

double Uniform(double);
double Gaussian(double);

double wfgauss(double);
double wfgauss2(double);
double wfgauss_d2(double);

double V(double);

double Error(double, double,int);
