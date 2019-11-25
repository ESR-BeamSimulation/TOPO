#include <iostream>
#include <fstream>
#include <time.h> 
#include "MethodMP.h"
#include "PIC.h"
#include "Pusher.h"
#include "LeapFrog.h"
#include "Global.h"
#include "iomanip"
#include <sys/time.h>   
#include <strstream>
# include <stdio.h>
# include <stdlib.h>
long getCurrentTime() 
{    
   struct timeval tv;    
   gettimeofday(&tv,NULL);    
   return tv.tv_sec * 1000 + tv.tv_usec / 1000;    
}    
  

//#include "RK4.h"

using std::cout;
using std::endl;

MethodMP::~MethodMP()
{}
MethodMP::MethodMP()
{}

void MethodMP::Push(Lattice &lattice,Mesh &mesh,Beam &beam)
{


  PIC pic(mesh,beam);

  Pusher *pusher;
  if(pusherType=="leapfrog")
  {
    pusher = new LeapFrog;
  }
  else if(pusherType == "rk4")
  {
    //pusher = new RK4;
  }
  
// print the the information of each step in output.dat  
  ofstream fout("output.dat");
  
    fout<<"  position(cm)  "
        <<"  EK(MeV)"
        <<"  xrms size(cm) "
        <<"  yrms size(cm) "
        <<"  zrms size(cm) "
        <<"  zrms size(deg) " 
        <<"  pxrms(mrad) " 
        <<"  pyrms(mrad) " 
        <<"  pzrms(mrad) "
        <<"  Nemitx(cm*mrad)"                
        <<"  Nemity(cm*mrad)"  
        <<"  Nemitz(cm*rad)"
        <<"  Nemitz(KeV*ns)  "
        <<"  betax(cm/rad)"
        <<"  betay(cm*rad)"
        <<"  betaz(cm*rad)"
        <<"  alphax       "
        <<"  alphay       "
        <<"  alphaz       "
        <<"  transm       "
        <<"  MBtransm     "
        <<endl;
        
    ofstream trajectory ("traj.dat");
    
    trajectory<<"partic index  "
              <<"transLoss     "
              <<"longiLoss     "
              <<"x  (cm)       "
              <<"px (mrad)     " 
              <<"y  (cm)       "
              <<"py (mrad)     "
              <<"z (cm)        "
              <<"pz (mrad)     "
              <<"SCEx (V/m)    "
              <<"SCEy (V/m)    "
              <<"SCEz (V/m)    "
              <<"ExtEx (V/m)   "
              <<"ExtEy (V/m)   "
              <<"ExtEz (V/m)   "
              <<"ExtBx (V/m)   "
              <<"ExtBy (V/m)   "
              <<"ExtBz (V/m)   "
              <<endl;
 
 

  double time1, time2;
  int nprint, printnumber; 
  ostringstream s1;
  string s2;


  for(int i = 0; i<stepMax; ++i)
  {
    nprint   = i % printStep;
    
    if (nprint==0)   // print the date at given time step
    {
        printnumber =   int (i / printStep); 
        s1 << "fort" << printnumber << ".dat";
        s2 = s1.str();
        s1.str("");
        s1.clear();
        ofstream output (s2);
        output<<"transLossFlag "
              <<"longiLoss     "
              <<"x  (cm)       "
              <<"px (mrad)     " 
              <<"y  (cm)       "
              <<"py (mrad)     "
              <<"z (cm)        "
              <<"pz (mrad)     "
              <<"SCEx (V/m)    "
              <<"SCEy (V/m)    "
              <<"SCEz (V/m)    "
              <<"ExtEx (V/m)   "
              <<"ExtEy (V/m)   "
              <<"ExtEz (V/m)   "
              <<"ExtBx (V/m)   "
              <<"ExtBy (V/m)   "
              <<"ExtBz (V/m)   "
              <<endl;
        
        for(int i=0;i<beam.particle.size();++i)
        {
            output<<setw(8) <<        beam.particle[i].transLossFlag                                               //1
                  <<setw(8) <<        beam.particle[i].longiLossFlag                                               //2  
                  <<setw(15)<<100 *  (beam.particle[i].xyz[0] - beam.average[0])                                   //3 
                  <<setw(15)<<1000*  (beam.particle[i].xyz[1] - beam.average[1]) / beam.particle[i].xyz[5]         //4 
                  <<setw(15)<<100 *  (beam.particle[i].xyz[2] - beam.average[2])                                   //5 
                  <<setw(15)<<1000*  (beam.particle[i].xyz[3] - beam.average[3]) / beam.particle[i].xyz[5]         //6
                  <<setw(15)<<100 *  (beam.particle[i].xyz[4] - beam.average[4])                                   //7
                  <<setw(15)<<1000*  (beam.particle[i].xyz[5] - beam.average[5]) /beam.particle[i].xyz[5]          //8
                  <<setw(15)<<        beam.particle[i].internalElecField[0]                                        //9
                  <<setw(15)<<        beam.particle[i].internalElecField[1]                                        //10
                  <<setw(15)<<        beam.particle[i].internalElecField[2]                                        //11
                  <<setw(15)<<        beam.particle[i].externalElecField[0]                                        //12
                  <<setw(15)<<        beam.particle[i].externalElecField[1]                                        //13
                  <<setw(15)<<        beam.particle[i].externalElecField[2]                                        //14
                  <<setw(15)<<        beam.particle[i].externalMegnField[0]                                        //15
                  <<setw(15)<<        beam.particle[i].externalMegnField[1]                                        //16
                  <<setw(15)<<        beam.particle[i].externalMegnField[2]                                        //17
                  <<setw(15)<<endl;
        }
        output.close(); 
    }
    


  
   for(int i=0;i<traj;++i)
   {
        trajectory<<setw(8)<<         i
                  <<setw(8)<<         beam.particle[i].transLossFlag
                  <<setw(8)<<         beam.particle[i].longiLossFlag
                  <<setw(15)<<100 *  (beam.particle[i].xyz[0]                  )
                  <<setw(15)<<1000*  (beam.particle[i].xyz[1]                  ) / beam.particle[i].xyz[5]    
                  <<setw(15)<<100 *  (beam.particle[i].xyz[2]                  )
                  <<setw(15)<<1000*  (beam.particle[i].xyz[3]                  ) / beam.particle[i].xyz[5] 
                  <<setw(15)<<100 *  (beam.particle[i].xyz[4]                  )
                  <<setw(15)<<1000*  (beam.particle[i].xyz[5] - beam.average[5]) / beam.particle[i].xyz[5]
                  <<setw(15)<<        beam.particle[i].internalElecField[0]
                  <<setw(15)<<        beam.particle[i].internalElecField[1]
                  <<setw(15)<<        beam.particle[i].internalElecField[2]
                  <<setw(15)<<        beam.particle[i].externalElecField[0]
                  <<setw(15)<<        beam.particle[i].externalElecField[1]
                  <<setw(15)<<        beam.particle[i].externalElecField[2]
                  <<setw(15)<<        beam.particle[i].externalMegnField[0]
                  <<setw(15)<<        beam.particle[i].externalMegnField[1]
                  <<setw(15)<<        beam.particle[i].externalMegnField[2] 
                  <<endl; 
    }


        
    time1=getCurrentTime();

    if(spaceChargeFlag==1)
    {
        pic.        GetInternalFiled      (mesh,  beam, lattice);     //InternalField.getfield();
    }
    lattice.      GetExternalField      (timeNow,  beam );          //ExternalField.getfield();

    beam.         SetField              (               );
    pusher->      Update                (timeStep, beam, lattice );                   //Pusher.update();
    lattice.      LostCount             (beam           );
    beam.         CaculateEmittance     (               );
    beam.         Print                 (fout           );
    timeNow     +=timeStep;
    
    
    time2=getCurrentTime();
 
//    win = (sqrt(1 + pow(beam.average[5],2) + pow(beam.average[3],2) + pow(beam.average[1],2)) -1) * BaseMassInMeV;
    win = (sqrt(1 + pow(beam.particle[0].xyz[5],2) + pow(beam.particle[0].xyz[3],2) + pow(beam.particle[0].xyz[1],2)) -1) * BaseMassInMeV;
    if(i%50==0)
    {
        cout<<endl;
        cout<<"  time       "
            <<"  posit(cm)  "
            <<"  EK(MeV)    " 
            <<"  rms x(cm)  "
            <<"  rms y (cm) "
            <<"  rms z (cm) "
            <<"  emitx (cm*mrad)  "
            <<"  emity (cm*mrad)  "
            <<"  emitz (cm*mrad)  "
            <<"  transm     "<<endl;
    }
    
    
    cout<<i<<"     "
        <<setw(15)<<100*beam.average[4]
        <<setw(15)<<win 
        <<setw(15)<<100*beam.sigma[0]  
        <<setw(15)<<100*beam.sigma[2]
        <<setw(15)<<100*beam.sigma[4]
        <<setw(15)<<10e5*beam.emittance[0]
        <<setw(15)<<10e5*beam.emittance[1]
        <<setw(15)<<10e5*beam.emittance[2]
//        <<setw(15)<<beam.emittance[2] * BaseMassInMeV *10e3 /C_light *10e9
        <<setw(15)<<beam.particleLiveNumber
        <<setw(15)<<beam.particleLiveInMesh
        <<setw(15)<<time2-time1
        <<setw(15)<<lattice.element[beam.particle[0].indexElem]->name
        <<endl;
    


    
    if(beam.average[4]>lattice.positionBack.back())   //beam are out of lattice
    {
      cout<<"lattice is ran out, and the total length is  "<<lattice.positionBack.back()<<"m"<<endl;
      break;
    }
    
    if(beam.particleLiveNumber==0)
    {
       cout<<"all particle are lost!!!  "<<endl;
       break ;
    }
  }
  
  fout.close();  
  trajectory.close();
  delete pusher;

}
