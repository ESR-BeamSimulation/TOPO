#include "Drift.h"
#include <iostream>
#include "Eigen/Dense"
using namespace Eigen;
using namespace std;


void Drift::GetField(double time,Particle & particle, double synPhase)
{
  particle.externalElecField = {0,0,0};
  particle.externalMegnField = {0,0,0};
}

void Drift::GetMatrix(Particle & particle)
{
    double gamma2; 
    
    gamma2   =   1+ pow(1+particle.xyz[5],2); 
    
    transferMatrixM.resize(6,6);

    transferMatrixM.row(0)  <<   1, length/gamma2, 0, 0,               0,  0;
    transferMatrixM.row(1)  <<   0, 1,             0, 0,               0,  0;
    transferMatrixM.row(2)  <<   0, 0,             1, length/gamma2,   0,  0;
    transferMatrixM.row(3)  <<   0, 0,             0, 1,               0,  0;
    transferMatrixM.row(4)  <<   0, 0,             0, 0,               1,  length/gamma2;
    transferMatrixM.row(5)  <<   0, 0,             0, 0,               0,  1;

    return;
}

void Drift::SetParameter(double posit, const vector<string> & param)
{  
  if(param.size()<4)
  {
    cout<<"Error! parameter not enough in element : \n"<<endl;
    for(int i = 0;i<param.size();++i)
    {
      cout<<param.at(i)<<"\t";
    }
    cout<<endl;
    return;
  }
  SetBasicParameter(posit,param);
}

//double Drift::getradius(double z)
//{
//    return aperturex;

//}
