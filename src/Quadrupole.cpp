#include "Quadrupole.h"
#include "Global.h"
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;


void Quadrupole::GetField(double time,Particle &particle,double synPhase)
{
  double r, theta;
  double br, bt;
  double bx,by;
  double G2;
  
  r     =   sqrt  ( pow(particle.xyz[0],2) + pow(particle.xyz[2],2) );
  theta =   atan2 ( particle.xyz[2], particle.xyz[0]);

  G2    =   gradient;
  
  br    =   G2 *r * sin (2 * theta);
  bt    =   G2 *r * cos (2 * theta);    
  
  
  bx    =  br * cos(theta) - bt * sin(theta);
  by    =  br * sin(theta) + bt * cos(theta);
  
    
  
  particle.externalElecField = { 0, 0, 0};
  particle.externalMegnField = {bx,by, 0};
  
  
  
//        for(int j=0;j<6;++j)
//        {
//            cout<<particle.xyz[j]<<"\t";
//        }
//        cout<<endl;
//        for(int j=0; j<3; ++j)
//        {
//           cout<<particle.externalMegnField[j]<<"\t";
//        }
//        cout<<endl;
//        
//        getchar();
  
  
}

void Quadrupole::GetMatrix(Particle &particle)
{

    double gamma2; 
    
    gamma2   =   1+ pow(1+particle.xyz[5],2); 
    
    transferMatrixM.resize(6,6);
    
    double G2, k, theta;
    G2  =  gradient;


    
    if(G2>0)
    {
        k       =   G2 / (BaseMassInMeV*10e6 *particle.xyz[5] / particle.qOverMass / BaseCharge / C_light);
        theta   =   k  * length ;
        transferMatrixM.row(0)  <<    cos(theta), sin(theta)/k,       0,                0,               0,  0;
        transferMatrixM.row(1)  << -k*sin(theta), cos(theta)  ,       0,                0,               0,  0;
        transferMatrixM.row(2)  <<             0,            0,      cosh(theta),       sinh(theta)/k,   0,  0;
        transferMatrixM.row(3)  <<             0,            0,    k*sinh(theta),       cosh(theta)  ,   0,  0;
    }
    else if(G2==0)
    {
        transferMatrixM.row(0)  <<    1, length,       0,                0,               0,  0;
        transferMatrixM.row(1)  <<    0,      1,       0,                0,               0,  0;
        transferMatrixM.row(2)  <<    0,      0,       1,           length,               0,  0;
        transferMatrixM.row(3)  <<    0,      0,       0,               1 ,               0,  0;
    
    }
    else if(G2<0)
    {
        k       =  -G2 / (BaseMassInMeV*10e6 *particle.xyz[5] / particle.qOverMass / BaseCharge / C_light);
        theta   =   k  * length ;
        transferMatrixM.row(0)  <<   cosh(theta), sinh(theta)/k,       0,                0,               0,  0;
        transferMatrixM.row(1)  << k*sinh(theta), cosh(theta)  ,       0,                0,               0,  0;
        transferMatrixM.row(2)  <<             0,            0,      cos(theta),       sin(theta)/k,   0,  0;
        transferMatrixM.row(3)  <<             0,            0,   -k*sin(theta),       cos(theta)  ,   0,  0;
    
    }
    
   transferMatrixM.row(4)  <<             0,            0,      0,       0,   1,  length/gamma2;
   transferMatrixM.row(5)  <<             0,            0,      0,       0,   0,  1;
    
}

void Quadrupole::SetParameter(double posit, const vector<string> & param)
{
  if(param.size()<5)
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
  gradient =    stod( param.at(4) );
}

