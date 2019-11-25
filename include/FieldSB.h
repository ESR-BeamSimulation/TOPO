#ifndef FIELDSB_H 
#define FIELDSB_H

#include <string>
#include <vector>

#include "Element.h"
#include "Beam.h"
using std::vector;
using v1d = std::vector<double>;
using v2d = vector<vector<double> > ;

class FieldSB:public Element
{

public:
  ~FieldSB();
  FieldSB();
  void SetParameter(double posit, const vector<string> & param);
  void GetField(double time,Particle &,double synPhase);


private:
  string filename;
  int type=1;
  double kb;
  double locationEnd;
  double locationMax;
  v1d zero3d=v1d{0,0,0};
  v1d zero6d=v1d{0,0,0};
  v2d    field;
  double gridNum;
  

  int nr,nz;//z had been defined before
  double rmax, zmax;
  double spacer, spacez;

  

};

#endif
