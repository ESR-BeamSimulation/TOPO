#ifndef ELEMENT_H
#define ELEMENT_H

#include <string>
#include <vector>
#include "Particle.h"
#include "Eigen/Dense"
using namespace Eigen;
using namespace std;
class Element
{
  public:
    virtual void GetField(double time,Particle &,double synPhase)       {}
    virtual void GetMatrix(Particle &)                                  {}
    virtual void SetParameter(double posit, const vector<string> & )    {}
    void         SetBasicParameter(double posit, const vector<string> & );
    virtual double getradius(double z)                                  {}
    virtual double Getfrequency ()                                      {}
    
    vector<string>      parameterString;        //store the information of element in one line in "lattice.txt"
    
    double              apertureX;              //to store the aperture of element
    double              apertureY;
    int                 shape;                  //0 means a round element, 1: means a rect element
    string              name="drift";
    double              position;               //to store the start position of each element
    
    MatrixXd             transferMatrixM = MatrixXd::Zero(6,6); 
 
 
  protected:
    vector<double>      parameter;


    double              length;                 //to store the length of element
};
#endif
