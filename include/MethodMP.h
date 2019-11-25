#ifndef METHODMP_H
#define METHODMP_H

#include <string>
#include "Method.h"
#include "Beam.h"
#include "Lattice.h"
#include "Mesh.h"


class MethodMP: public Method
{
  public: 
   ~MethodMP();
    MethodMP();

    void Push(Lattice &,Mesh &,Beam &);
    
    double timeStep=1,timeNow=0;
    int    stepMax =1;
    int    spaceChargeFlag=0;   
    double win;
    int    printStep; 
    int     traj;
    string pusherType;

  private:
};

#endif
