#include "LeapFrog.h"
#include <iostream>
#include <Global.h>

void LeapFrog::Update(double timeStep, Beam &beam, Lattice &lattice)
{
  int    particleNumber = beam.particle.size();
//#pragma omp for
    
    double pxTemp =0,pyTemp=0,pzTemp=0;
    double cc=0,megnField=0;
    double gamma=0; 
    double pBxTemp=0,pByTemp=0,pBzTemp=0;
    double tempB=0;
    double rgamma=0;


    for (int i=0;i<particleNumber;++i)
    {
        if (    (lattice.element[beam.particle[i].indexElem]->name) == "sol"
              ||(lattice.element[beam.particle[i].indexElem]->name) == "sol3d"
              ||(lattice.element[beam.particle[i].indexElem]->name) == "filedsb" )
        {
            beam.particle[i].xyz[1] = beam.particle[i].xyz[1] + beam.particle[i].externalMegnField[2] * C_light / (BaseMassInMeV *10e6) /2.0 * beam.particle[i].xyz[2];
            beam.particle[i].xyz[3] = beam.particle[i].xyz[3] - beam.particle[i].externalMegnField[2] * C_light / (BaseMassInMeV *10e6) /2.0 * beam.particle[i].xyz[0];
            beam.particle[i].xyz[5] = beam.particle[i].xyz[5];
        }
        
    }
    





//    if(beam.particle[2].totalMegnField[0]!=0)
//    {
//        for(int j=0;j<6;++j)
//        {
//            cout<<beam.particle[2].xyz[j]<<"\t";
//        }
//        cout<<endl;
//        for(int j=0; j<3; ++j)
//        {
//           cout<<beam.particle[2].totalMegnField[j]<<"\t";
//        }
//        cout<<endl;
//        
//        cout<<sqrt(1 + pow(beam.particle[2].xyz[1],2) +pow(beam.particle[2].xyz[3],2) +pow(beam.particle[2].xyz[5],2))<<endl;
//        
//        getchar();
//    }


  for(int i=0;i<particleNumber;++i)
  {
    if(beam.particle[i].transLossFlag)  //if lost transversely, continue otherwise do the loop
    {
      continue;
    }

    //1:half-step accleration in the electric field    
    gamma  =  sqrt(1 + pow(beam.particle[i].xyz[1],2) +pow(beam.particle[i].xyz[3],2) +pow(beam.particle[i].xyz[5],2));
    for(int j=0;j<3;++j)
    {
      beam.particle[i].xyz[2*j+1]+= beam.particle[i].qOverMass /gamma * beam.particle[i].totalElecField[j] /C_light *timeStep/2;        
    }



    //2:particle velocity rotate in magnetic field
    


    gamma   =   sqrt(1 + pow(beam.particle[i].xyz[1],2) +pow(beam.particle[i].xyz[3],2) +pow(beam.particle[i].xyz[5],2));

    pBxTemp =   beam.particle[i].totalMegnField[0] *beam.particle[i].qOverMass *timeStep/gamma;
    pByTemp =   beam.particle[i].totalMegnField[1] *beam.particle[i].qOverMass *timeStep/gamma;
    pBzTemp =   beam.particle[i].totalMegnField[2] *beam.particle[i].qOverMass *timeStep/gamma;
    tempB   =   1.0 / (1 + ( pow(pBxTemp,2)+ pow(pByTemp,2)+ pow(pBzTemp,2) )/4.0 );
    

      pxTemp=beam.particle[i].xyz[1];
      pyTemp=beam.particle[i].xyz[3];
      pzTemp=beam.particle[i].xyz[5];
      
      
      beam.particle[i].xyz[1] = ( 1   +   (pow(pBxTemp,2) - pow(pByTemp,2) -  pow(pBzTemp,2)   ) /4.0) * pxTemp * tempB
                                      +   (    pBzTemp    +     pBxTemp    *      pByTemp        /2.0) * pyTemp * tempB
                                      +   (   -pByTemp    +     pBxTemp    *      pBzTemp        /2.0) * pzTemp * tempB;  
                 
      beam.particle[i].xyz[3] =           (   -pBzTemp    +     pBxTemp    *      pByTemp        /2.0) * pxTemp * tempB 
                              + ( 1   +   (pow(pByTemp,2) - pow(pBzTemp,2) -  pow(pBxTemp,2)   ) /4.0) * pyTemp * tempB
                                      +   (    pBxTemp    +     pByTemp    *      pBzTemp        /2.0) * pzTemp * tempB; 

      beam.particle[i].xyz[5] =           (    pByTemp    +     pBzTemp    *      pBxTemp        /2.0) * pxTemp * tempB  
                                      +   (   -pBxTemp    +     pByTemp    *      pBzTemp        /2.0) * pyTemp * tempB
                              + ( 1   +   (pow(pBzTemp,2) - pow(pByTemp,2) -  pow(pBxTemp,2)   ) /4.0) * pzTemp * tempB;






    //3:curvilinear accleration
    //This is a linac disign,so the R is infinite,so this step won't make any change




    //4:half-step accleration in the electric field again
    gamma  =  sqrt(1 + pow(beam.particle[i].xyz[1],2) +pow(beam.particle[i].xyz[3],2) +pow(beam.particle[i].xyz[5],2));
    for(int j=0;j<3;++j)
    {
      beam.particle[i].xyz[2*j+1]+= beam.particle[i].qOverMass / gamma* beam.particle[i].totalElecField[j] /C_light *timeStep/2;  //unit:m/s
    }

    //5:advance the position
    gamma   =   sqrt(1 + pow(beam.particle[i].xyz[1],2) +pow(beam.particle[i].xyz[3],2) +pow(beam.particle[i].xyz[5],2));

    beam.particle[i].xyz[0] += C_light /gamma *beam.particle[i].xyz[1]*timeStep;                  //unit:m
    beam.particle[i].xyz[2] += C_light /gamma *beam.particle[i].xyz[3]*timeStep;
    beam.particle[i].xyz[4] += C_light /gamma *beam.particle[i].xyz[5]*timeStep;
  }
  
//    if(beam.particle[2].totalMegnField[0]!=0)
//    {
//        for(int j=0;j<3;++j)
//        {
//            cout<<beam.particle[2].xyz[2*j+1]<<"\t";
//        }
//        cout<<endl;
//        cout<<sqrt(1 + pow(beam.particle[2].xyz[1],2) +pow(beam.particle[2].xyz[3],2) +pow(beam.particle[2].xyz[5],2))<<endl;
//        getchar();
//      }
 

    for (int i=0;i<particleNumber;++i)
    {
        if (    (lattice.element[beam.particle[i].indexElem]->name) == "sol"
              ||(lattice.element[beam.particle[i].indexElem]->name) == "sol3s"
              ||(lattice.element[beam.particle[i].indexElem]->name) == "filedsb" )
        {
            beam.particle[i].xyz[1] = beam.particle[i].xyz[1] - beam.particle[i].externalMegnField[2] * C_light / (BaseMassInMeV *10e6) /2.0 * beam.particle[i].xyz[2];
            beam.particle[i].xyz[3] = beam.particle[i].xyz[3] + beam.particle[i].externalMegnField[2] * C_light / (BaseMassInMeV *10e6) /2.0 * beam.particle[i].xyz[0];
            beam.particle[i].xyz[5] = beam.particle[i].xyz[5];
        }
        
    }
 
 
 
 
 
 
  
      
}

