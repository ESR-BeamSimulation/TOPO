#include <iostream>
#include "Particle.h"
#include "Global.h"

Particle::Particle()
{
  xyz.                  resize(6);
  externalElecField.    resize(3);
  externalMegnField.    resize(3);
  internalElecField.    resize(3);
  internalMegnField.    resize(3);
  totalElecField.       resize(3);
  totalMegnField.       resize(3);
}

void Particle::GetTotalField()       
{

  for(int i = 0; i<3 ;++i)
  {    
    totalElecField[i] = (externalElecField[i] + internalElecField[i]) ; 
    totalMegnField[i] = (externalMegnField[i] + internalMegnField[i]) ;
  }
}

void Particle::GetMatrix()         // to obtain the transfer matrix for each particle
{
//   particleTranferMap.      resize(6);
//   particleTranferMap.cols() = 6;
   
   
}
