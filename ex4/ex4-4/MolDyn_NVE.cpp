/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

//Esercizi 4: Meccanica Statistica -> Dinamica Molecolare in microcanonico NVE

//4.1: implementare un'opzione di restart e rescale delle velocità (Equilibrazione!)
//4.2: aggiungere il calcolo delle medie con incertezze
//4.3: calcolare P
//4.4: applicarlo a casi di Argon e Kripton


#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <string>
#include "MolDyn_NVE.h"
#include "random.h"


using namespace std;

double Error(double, double, int);

//*************************************************************************************************//
//MAIN

int main(){
  
  Input();             //Inizialization
  int nconf = 1;
  
  for(int n=0; n<nblock; n++){
    
    Reset(n);
    
    for(int l=0; l<L; l++) {
          int istep=n*L+l+1;
          Move();           //Move particles with Verlet algorithm
           if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
           if(istep%10 == 0){
             Measure();     //Properties measurement
             Accumulate();
      //        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
              nconf += 1;
           }
      }
    Averages(n);
  }
  ConfOld();   //Per verlet mi serve anche posizione vecchia
  ConfFinal();  //Write final configuration to restart


  return 0;
}

//*************************************************************************************************//
//SET SIMULATION PARAMETERS:

void Input(void){ //Prepare all stuff for the simulation
  
  ifstream ReadInput;
  //double ep, ek, pr, et, vir;
  cout << "-----------------------------------------------------" << endl;
  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;
  cout << "-----------------------------------------------------" << endl<<endl;
  
//READ SIMULATION FILE **************************************************
  ReadInput.open("input.simulation"); //Read input
  
  ReadInput >> element; //settare a 0 se no real simul, altrimenti direttamente nome dell'elemento
  ReadInput >> phase;
  ReadInput >> restart;
  ReadInput >> rescale;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> nblock;
  ReadInput >> iprint;
  
    L=nstep/nblock;
  
  ReadInput.close();
  ReadInput.clear();
  
  /*string input_phase= "input."+ phase;
  string input_element= element+ "/input.element";*/
  
//READ ELEMENT FILE:**************************************************
  if (element!="0"){
    
    cout << "Simulation of the real element: " << element <<", in its "<< phase <<" phase." << endl <<endl;
    cout << "Results in SI units will be available in the element folder."<<endl<<endl;
    
    ReadInput.open(element+ "/input.element");
    
    if(!ReadInput.is_open()) {
      cerr << "PROBLEM: file "<< element <<"/input.element does not exit." << endl;
      exit(1);
    }
    
    ReadInput >> sigma;
    ReadInput >> epsilon_kb;
    ReadInput >> mass;
    
    epsilon=epsilon_kb*kb;
    
    cout << element <<" characteristic parameters in Lennard-Jones model are:" <<endl;
    cout << "– Lenght unit σ= " <<sigma <<" nm" <<endl;
    cout << "– Energy unit ε= " <<epsilon <<" J"<<endl;
    cout << "– Mass of particles m= "<< mass << " amu"<< endl<<endl;;
    
    sigma*=nm;
    mass*=amu;
    
    fp=epsilon/(pow(sigma,3));
    
    ReadInput.close();
    ReadInput.clear();
  }else{ cout << "Results are in LJ units" << endl<<endl;}
  
//READ PHASE FILE**************************************************
  ReadInput.open("input."+ phase);
  
  if(!ReadInput.is_open()) {
    cerr << "PROBLEM: file input."<< phase<< " does not exit." << endl;
    exit(1);
  }
  
  ReadInput >> temp;
  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;
  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;
  ReadInput >> rcut;

  ReadInput.close();
  
  cout <<endl<< "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl;
  cout << "Number of blocks = " << nblock <<endl <<endl;
  
//*******************************************************************
  
  if (!restart && rescale) restart=true; // Nel caso uno si dimentichi di settare il restart

  if (restart){
    cout << "Restart simulation";
    if (rescale) cout << " and rescale velocities" <<endl;
    cout << endl;
  }
  
  //First run? --> read from config.0 and random velocities
  //Restart simulation? --> read from old configurations
  //while rescaling velocities to equilibrate the system? --> use actual T and target T
  
  if(!restart){
    system("bash clean.sh");
    InitialConfig();
    InitialVelocities();
  }else{
    RestartConfig();
    if(rescale){
      RescaleVelocities();
    }else{
      InitialVelocities();
   }
  }

   return;
}

//*************************************************************************************************//
//READ CONFIGURATION:

void InitialConfig(void){//Read initial configuration
    ifstream ReadConf;
cout << "Read initial configuration from file config.0 " << endl << endl;
ReadConf.open("config.0");
for (int i=0; i<npart; ++i){
  ReadConf >> x[i] >> y[i] >> z[i];
  x[i] = x[i] * box;
  y[i] = y[i] * box;
  z[i] = z[i] * box;
}
ReadConf.close();
  return;
}

void RestartConfig(void){ //Read old configuration to restart simulation
ifstream ReadConf;
  
  ReadConf.open("old.0");
  
  if(!ReadConf.is_open()) {
    cerr << "PROBLEM: file old.0 does not exit." << endl << "You may need to run one first time the simulation. " << endl;
    exit(1);
  }
  
  for (int i=0; i<npart; ++i){
      ReadConf >> xold[i] >> yold[i] >> zold[i];
      xold[i] = xold[i] * box;
      yold[i] = yold[i] * box;
      zold[i] = zold[i] * box;
  }
  ReadConf.close();
  ReadConf.clear();

  ReadConf.open("old.final");
  
  if(!ReadConf.is_open()) {
    cerr << "PROBLEM: file old.final does not exit." << endl << " You may need to run one first time the simulation. " << endl;
    exit(1);
  }else {cout << "Read initial configuration from file old.final " << endl << endl;}
  
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }

ReadConf.close();
  
   return;
}

//*************************************************************************************************//
//SPATIAL MOTION OF PARTICLES (VERLET ALGORITHM)
//--> cambia config e ricalcola velocità, conta che forza cambia perchè cambia potenziale perchè cambia configuraione(più distanti)
//--> in pratica più la configurazione si sparpaglia più rallenta la diffusione

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}


//*************************************************************************************************//
//SET VELOCITIES:

void InitialVelocities(void){ //Prepare initial velocities
  
  rnd.SetSeed();

  cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
  double sumv[3] = {0.0, 0.0, 0.0};
  for (int i=0; i<npart; ++i){
    vx[i] = rnd.Rannyu(-0.5,0.5); //rnd.Rannyu() - 0.5;
    vy[i] = rnd.Rannyu(-0.5,0.5); //rnd.Rannyu() - 0.5;
    vz[i] = rnd.Rannyu(-0.5,0.5);

    sumv[0] += vx[i];
    sumv[1] += vy[i];
    sumv[2] += vz[i];
  }
  for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
  double sumv2 = 0.0, fs;
  for (int i=0; i<npart; ++i){
    vx[i] = vx[i] - sumv[0];
    vy[i] = vy[i] - sumv[1];
    vz[i] = vz[i] - sumv[2];

    sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
  }
  sumv2 /= (double)npart;

  fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
  for (int i=0; i<npart; ++i){
    vx[i] *= fs;
    vy[i] *= fs;
    vz[i] *= fs;

    xold[i] = Pbc(x[i] - vx[i] * delta);
    yold[i] = Pbc(y[i] - vy[i] * delta);
    zold[i] = Pbc(z[i] - vz[i] * delta);
  }
  
  rnd.SaveSeed();
  return;
}


void RescaleVelocities(void){
  
    cout << "Rescale velocities from old configuration" << endl << endl;
  
  double sumv2 = 0.0, fs;
  
  //1. After reading old configurations r(t-dt) and r(t), evaluate r(t+dt) with Verlet algorithm
  Move();
  //2. now 'rold'=r(t) and 'r'=r(t+dt), so evaluate new velocities from the two configurations

  for(int i=0; i<npart; ++i){
    sumv2+=vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
  }
  
  //3. actualT=<v2>/3 and fs=√(targetT/actualT)--> ma vado diretta al rescaling senza passare per calcolo esplicito T
  //aggiungere cout di verifica con Tactual
  
  sumv2 /= (double)npart;
  fs = sqrt(3 * temp / sumv2);
   
  //4. rescale velocities
  for (int i=0; i<npart; ++i){
    vx[i] *= fs;
    vy[i] *= fs;
    vz[i] *= fs;
      
  //5. evaluate new 'rold' to restart configuration:
    double x_i=xold[i];
    double y_i=yold[i];
    double z_i=zold[i];
    
    xold[i] = Pbc(x[i] - (vx[i])*2.*delta);
    yold[i] = Pbc(y[i] - (vy[i])*2.*delta);
    zold[i] = Pbc(z[i] - (vz[i])*2.*delta);
     
    x[i]=x_i;
    y[i]=y_i;
    z[i]=z_i;
  }
     return;
}

//*************************************************************************************************//
//MEASURE (INSTANT) THERMODINAMIC PROPERTIES:

void Measure(){ //Properties measurement
  //int bin;
  double v=0., t=0., w=0.;
  double vij, wij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp, Press;

  Epot.open(phase+"/"+"instant_epot.out",ios::app);
  Ekin.open(phase+"/"+"instant_ekin.out",ios::app);
  Etot.open(phase+"/"+"instant_etot.out",ios::app);
  Temp.open(phase+"/"+"instant_temp.out",ios::app);
  Press.open(phase+"/"+"instant_press.out",ios::app);

  /*v = 0.0; //reset observables t = 0.0;*/

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       wij = 48.0/pow(dr,12) - 24.0/pow(dr,6);

//Potential energy (extensive)
       v += vij;
       w += wij;
     }
    }
  }

//Kinetic energy (extensive)
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
  
//INSTANT INTENSIVE PROPERTIES:
  stima_pot = v/(double)npart; //Potential energy per particle
  stima_kin = t/(double)npart; //Kinetic energy per particle
  stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
  stima_etot = (t+v)/(double)npart; //Total energy per particle
  stima_press = rho*stima_temp + w/(vol*3.); //Pressure per particle
  
  if (element!="0"){
    
    Epot.close();  Epot.clear();
    Ekin.close();  Ekin.clear();
    Etot.close();  Etot.clear();
    Temp.close();  Temp.clear();
    Press.close(); Press.clear();
    
    string results=element+"/"+phase+"/";
    Epot.open(results+"instant_epot.out",ios::app);
    Ekin.open(results+"instant_ekin.out",ios::app);
    Etot.open(results+"instant_etot.out",ios::app);
    Temp.open(results+"instant_temp.out",ios::app);
    Press.open(results+"instant_press.out",ios::app);
    
    stima_pot *= epsilon;
    stima_kin *= epsilon;
    stima_etot *= epsilon;
    stima_temp *= epsilon_kb;
    stima_press *= fp;
    
  }
  
  Epot << stima_pot  << endl;
  Ekin << stima_kin  << endl;
  Etot << stima_etot << endl;
  Temp << stima_temp << endl;
  Press << stima_press << endl;
  
  Epot.close();
  Ekin.close();
  Etot.close();
  Temp.close();
  Press.close();
  
    return;
}
//*************************************************************************************************//
//AVERAGE:
void Reset(int n){
  for(int i=0; i<n_props; i++) {
    mb[i]=0;
    l_norm=0;
    if (n==0) {m[i]=0.;m2[i]=0.;err[i]=0.;}}
  
}

void Accumulate(void){
  mb[iv]+=stima_pot;
  mb[ik]+=stima_kin;
  mb[ie]+=stima_etot;
  mb[it]+=stima_temp;
  mb[ip]+=stima_press;
  l_norm++;
}

void Averages(int n){
  
    ofstream Epot, Ekin, Etot, Temp, Press;
  
  Epot.open(phase+"/"+"ave_epot.out",ios::app);
  Ekin.open(phase+"/"+"ave_ekin.out",ios::app);
  Etot.open(phase+"/"+"ave_etot.out",ios::app);
  Temp.open(phase+"/"+"ave_temp.out",ios::app);
  Press.open(phase+"/"+"ave_press.out",ios::app);
  
   // cout << "Block number " << n+1 << endl;
  
  for(int i=0; i<n_props;i++){
    mb[i]/=double(l_norm);
    m[i]=(m[i]*double(n)+mb[i])/double(n+1);
    m2[i]=(m2[i]*double(n)+mb[i]*mb[i])/double(n+1);
    err[i]=Error(m[i],m2[i], n);
  }
  
  if (element!="0"){
    
    Epot.close();  Epot.clear();
    Ekin.close();  Ekin.clear();
    Etot.close();  Etot.clear();
    Temp.close();  Temp.clear();
    Press.close(); Press.clear();
    
    string results=element+"/"+phase+"/";
    Epot.open(results+"ave_epot.out",ios::app);
    Ekin.open(results+"ave_ekin.out",ios::app);
    Etot.open(results+"ave_etot.out",ios::app);
    Temp.open(results+"ave_temp.out",ios::app);
    Press.open(results+"ave_press.out",ios::app);
    
    /*m[iv] *= epsilon;                     err[iv]*= epsilon;
    m[ik] *= epsilon;                     err[ik]*= epsilon;
    m[ie] *= epsilon;                     err[ie]*= epsilon;
    m[it] *= epsilon_kb;                  err[it]*= epsilon_kb;
    m[ip] *= fp;                          err[ip]*= fp;*/
    
  }
  
  Epot << (n+1)*l_norm << " " << m[iv] << " " << err[iv] << endl;
  Ekin << (n+1)*l_norm << " " << m[ik] << " " << err[ik] << endl;
  Etot << (n+1)*l_norm << " " << m[ie] << " " << err[ie] << endl;
  Temp << (n+1)*l_norm << " " << m[it] << " " << err[it] << endl;
  Press<< (n+1)*l_norm << " " << m[ip] << " " << err[ip] << endl;
  
  Epot.close();
  Ekin.close();
  Etot.close();
  Temp.close();
  Press.close();
}

//*************************************************************************************************//
//PRINT CONFIGURATION

void ConfOld(void){
  ofstream WriteOld0("old.0");
  ofstream WriteOldF("old.final");
  
  cout << endl <<"Print final configuration to file old.0 and old.final to eventually restart simulation" << endl << endl;
  
  for (int i=0; i<npart; ++i) WriteOld0 << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  for (int i=0; i<npart; ++i) WriteOldF << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  
  WriteOld0.close();
  WriteOldF.close();
  
     return;
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
  return;
}

//*************************************************************************************************//
double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

double Error(double ave, double ave2,int n){ //all'errore passo già "n-1"
    if(n==0) return 0.;
    else return sqrt((ave2-ave*ave)/n);
}

 
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
