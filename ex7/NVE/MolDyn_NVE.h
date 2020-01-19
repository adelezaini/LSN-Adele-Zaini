/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//random numbers
#include "random.h"
Random rnd;

//parameters, observables
//––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
const int m_props=1000;
int n_props, iv, ik, it, ie, ip, igofr;
double bin_size;
int nbins;
double walker[m_props];
const double pi=3.1415927;
//––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
/*const int n_props=5;
int iv=0,ik=1,ie=2,it=3, ip=4;*/
double stima_pot, stima_kin, stima_etot, stima_temp, stima_pres;

// averages
//double acc,att;
int l_norm;
double m[m_props]={0},m2[m_props]={0},err[m_props]={0},mb[m_props];
//double m[n_props]={0},m2[n_props]={0},err[n_props]={0},mb[n_props]; // mb2[m_props];

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double temp,vol,rho,box,rcut; //energy,

// simulation
int nstep, iprint, nblock, L; //seed;
double delta;
bool restart, rescale; //autoequi;

//functions
void Input(void);

void InitialConfig(void);
void RestartConfig(void);

void Move(void);

void RescaleVelocities(void);
void InitialVelocities(void);

void Measure(void);
double Force(int, int);
double Pbc(double);
//void Equilibrate(void);
void Reset(int);
void Accumulate(void);
void Averages(int);

void ConfXYZ(int);
void ConfOld(void);
void ConfFinal(void);

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/


