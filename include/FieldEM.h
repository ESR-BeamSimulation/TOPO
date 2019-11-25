#ifndef  FIELDEM_H
#define  FIELDEM_H

#include <string>
#include <vector>
#include <Element.h>
#include "Beam.h"
using std::vector;
using v1d = std::vector<double>;
using v2d = vector<vector<double> > ;

class FieldEM : public Element
{
public:
  ~FieldEM();
  FieldEM();
  void SetParameter(double posit, const vector<string> & param);
  void GetField(double time,Particle &,double synPhase);
private:
  string filename;
  int type=1;
  double ke,kb;
  double locationEnd;
  double frq;
  v1d zero3d=v1d{0,0,0};
  v1d zero6d=v1d{0,0,0};
  
  int count=0;
  


  //type1
  int    nx,ny,nz;
  double locationMax;
  double xmin,xmax,ymin,ymax,zmax;
  v2d    field;
  double gridNum;
  double synPhase,offsetPhase;
  double spacex,spacey,spacez;



};
#endif
