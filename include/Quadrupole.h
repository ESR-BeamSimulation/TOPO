#ifndef QUADRUPOLE_H
#define QUADRUPOLE_H
#include "Element.h"
#include "Particle.h"
#include <iostream>
#include <vector>
#include <string>

class Quadrupole : public Element 
{
  public:
    void GetField(double time,Particle &, double synPhase);
    void GetMatrix(Particle &);
    void SetParameter(double posit,const vector<string> & );
  private:
    double   gradient;
};

#endif //QUADRUPOLE_H
