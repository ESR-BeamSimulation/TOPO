#ifndef FIELD_H 
#define FIELD_H

#include <string>
#include <vector>

#include "Element.h"
#include "Beam.h"
using std::vector;
using v1d = std::vector<double>;
using v2d = vector<vector<double> > ;

class Field:public Element
{
public:
  ~Field();
  Field();
  void SetParameter(double posit, const vector<string> & param);
  void GetField(double time,Particle &,double synPhase);
  void GetMatrix(Particle &);
  void Synchronize(double time,Beam &beam);
  void SetOffset(double time, Beam &beam);
  void ScanOffset(double time, Beam &beam);
  
  //Synchronize
  int    synchronized=0;
protected:
  v2d cor;
private:
  string filename;
  int type=1;
  double ke,kb;
  double locationEnd;
  double frq;
  v1d zero3d=v1d{0,0,0};
  v1d zero6d=v1d{0,0,0};
  
  int count=0;
  
  //type0 get field
  double d1,d2;
  v2d gridLocation;
  v2d eField;
  v2d mField;
  v1d fieldofParticle=v1d{0,0,0,0,0,0};

  //type1
  int    nx,ny,nz;
  double locationMax;
  double xmin,xmax,ymin,ymax,zmax;
  v2d    field;
  double gridNum;
  double synPhase,offsetPhase;

  double spacex,spacey,spacez;

  //type2
  int nr;	//z had been defined before
  double rmax;
  double spacer;

};

#endif
