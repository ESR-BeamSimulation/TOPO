#ifndef RFQ_H
#define RFQ_H
#include "Element.h"
#include "Particle.h"
#include <iostream>
#include <vector>

typedef vector<double> v1d;
typedef vector<vector<double> > v2d;
class RFQ : public Element 
{
  public:
    void GetField(double time,Particle &,double synPhase);
    void GetMatrix(Particle &);
    void SetParameter(double posit,const vector<string> & );
    double getradius(double z);
    double Getfrequency();
  protected:
    v2d factor; 
    v1d cellLength;
    v1d cellPosition;
    v1d cellPositionEnd;
    v1d radius;
    v1d UL;
    v1d mod;
    v1d aperture;
  private:
    double frequency;
    double phiatEntrance;
    vector<double> zero3d={0,0,0};
    double ERR=1.0e-16;
    int getEleNum(double z);

    
    //RFQ::getField private variable
    double synPhase;


};

#endif
