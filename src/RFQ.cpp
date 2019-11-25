#include <cmath>
#include <fstream>
#include "RFQ.h"
#include "MyFunc.h"
#include<vector>

double RFQ::Getfrequency () 
{
    return frequency;
}

void RFQ::GetField(double time,Particle &particle,double synPhase)
{
  double Er, Ez, Et, Ex, Ey;
  
  int  eleNum = getEleNum(particle.xyz[4]);
  if(eleNum==-1)
  {
      particle.externalElecField=zero3d;
      particle.externalMegnField=zero3d;
  }
  double theta = atan2( particle.xyz[2],    particle.xyz[0] );
  double radiu = sqrt ( particle.xyz[0]*particle.xyz[0] + particle.xyz[2]*particle.xyz[2]);
  double z     = particle.xyz[4]-cellPosition[eleNum];
  double k     = M_PI/cellLength[eleNum];
  
  double A10, A30, A03, A01, A12, A32, A21, A23;
  
         A01 = factor[eleNum][0];
         A10 = factor[eleNum][1];
         A03 = factor[eleNum][2];
         A12 = factor[eleNum][3];
         A21 = factor[eleNum][4];
         A23 = factor[eleNum][5];
         A30 = factor[eleNum][6];
         A32 = factor[eleNum][7];
  
  double Ul, r0;
         Ul  = UL[eleNum];
         r0 = radius[eleNum];       

  if(radiu<ERR)
  {
    Er =0.0;
    Et =0.0;
    Ez = -abs(Ul) /2.0 * (A10 *sin(  k*z) *Bessin(0,  k*radiu)     * k 
                        + A30 *sin(3*k*z) *Bessin(0,3*k*radiu) * 3 * k ) * pow(-1.0,eleNum);
  }
  else
  {
    double b41,b43,b22,b62;
    
//    b11=Bessin(1,  k*radiu);
//    b13=Bessin(1,3*k*radiu);
//    b51=Bessin(5,  k*radiu);
//    b53=Bessin(5,3*k*radiu);
//    b32=Bessin(3,2*k*radiu);
//    b72=Bessin(7,2*k*radiu);
//    b41=Bessin(4,  k*radiu);
//    b43=Bessin(4,3*k*radiu);
//    b22=Bessin(2,2*k*radiu);
//    b62=Bessin(6,2*k*radiu);
//    b01=Bessin(0,  k*radiu);
//    b03=Bessin(0,3*k*radiu);
    
    
    
    Er =   abs(Ul)/2.0* (   A01*2*radiu           / pow(r0,2) * cos(2*theta)
                         +  A03*6*pow(radiu,5)    / pow(r0,6) * cos(6*theta)
                        
                         +  A10 * Bessin(1,  k*radiu)  *k *cos(  k*z)               *pow(-1.0,eleNum)
                         +  A30 * Bessin(1,3*k*radiu)*3*k *cos(3*k*z)               *pow(-1.0,eleNum)
                         +  A12 * Bessin(5,  k*radiu)  *k *cos(  k*z) *cos(4*theta) *pow(-1.0,eleNum)
                         +  A32 * Bessin(5,3*k*radiu)*3*k *cos(3*k*z) *cos(4*theta) *pow(-1.0,eleNum)
                         +  A21 * Bessin(3,2*k*radiu)*2*k *cos(2*k*z) *cos(2*theta) *pow(-1.0,eleNum)
                         +  A23 * Bessin(7,2*k*radiu)*2*k *cos(2*k*z) *cos(6*theta) *pow(-1.0,eleNum));
  

    Et =  -abs(Ul)/2.0* (   A01 * pow(radiu/r0,2)                    *sin(2*theta)*2 
                         +  A03 * pow(radiu/r0,6)                    *sin(6*theta)*6 
                         
                         +  A12 * Bessin(4,  k*radiu) *cos(  k*z)    *sin(4*theta)  *4  *pow(-1.0,eleNum)
                         +  A32 * Bessin(4,3*k*radiu) *cos(3*k*z)    *sin(4*theta)  *4  *pow(-1.0,eleNum)
                         +  A21 * Bessin(2,2*k*radiu) *cos(2*k*z)    *sin(2*theta)  *2  *pow(-1.0,eleNum)
                         +  A23 * Bessin(6,2*k*radiu) *cos(2*k*z)    *sin(6*theta)  *6  *pow(-1.0,eleNum)
                         ) /radiu;
    
    Ez=   -abs(Ul)/2.0* (   A10 * Bessin(0,  k*radiu)*sin(  k*z)  *k
                         +  A30 * Bessin(0,3*k*radiu)*sin(3*k*z)*3*k 
                         +  A12 * Bessin(4,  k*radiu)*sin(  k*z)  *k *cos(4* theta)
                         +  A32 * Bessin(4,3*k*radiu)*sin(3*k*z)*3*k *cos(4* theta)
                         +  A21 * Bessin(2,2*k*radiu)*sin(2*k*z)*2*k *cos(2* theta)
                         +  A23 * Bessin(6,2*k*radiu)*sin(2*k*z)*2*k *cos(6* theta))*pow(-1.0,eleNum);

 
// test for two terms expression
//      Er =   abs(Ul)/2.0* (   A01 * 2*radiu /pow(r0,2) * cos(2*theta)
//                            
//                            + A10 * Bessin(1  ,k*radiu)*cos(  k*z)*k  *pow(-1.0,eleNum)     );
//                                 
//      Et =  -abs(Ul)/2.0*    A01 * pow(radiu/r0,2)  * sin(2*theta) *2 /radiu;
//      Ez =  -abs(Ul)/2.0* (  A10 *  Bessin(0,  k*radiu)*sin(k*z)*k )  *pow(-1.0,eleNum);   
  }
  
  
  
  
  
  
  double phaseCos=cos(2*M_PI*frequency*time + synPhase);
  
  Ex=Er*cos(theta)-Et*sin(theta);
  Ey=Er*sin(theta)+Et*cos(theta);
  
  particle.externalElecField[0] = Ex*phaseCos;
  particle.externalElecField[1] = Ey*phaseCos;
  particle.externalElecField[2] = Ez*phaseCos;
  particle.externalMegnField    = zero3d;
}

void RFQ::GetMatrix(Particle &)
{

}
void RFQ::SetParameter(double posit,const vector<string> & param)
{

  if(param.size()!=7)
  {
    cout<<"Error! parameters are wrong in RFQ set : \n"<<endl;
    
    for(int i = 0;i<param.size();++i)
    {
      cout<<param.at(i)<<"\t";
    }
    cout<<endl;
    return;
  } 
  

  
  SetBasicParameter(posit,param);
  frequency     = stod(param.at(4));
  synPhase     = stod(param.at(5))/180.0*M_PI;
 

  
  ifstream rfStream;
  rfStream.open(param.at(6));
  if(!rfStream.is_open())
  {
    cout<<"RFQ Structure path ERROR!"<<endl;
  }
  string str;
  vector<string> strVec;
  v1d vtemp(8);
  for(int i = 0;rfStream.peek()!=EOF;++i)
  {
    getline(rfStream,str);
    StringSplit(str,strVec);

    if(i==0)
    {
      cellPosition.push_back(position);    
    }
    else
    {
      cellPosition.push_back(cellPosition.back()+cellLength.back());
    }

    
    cellLength      .push_back(  stod(strVec[1])/100.0     );
    cellPositionEnd .push_back(  cellLength.back() + cellPosition.back() );
    UL              .push_back(  stod(  strVec[2]) *1e3    );    //kV to V
    radius          .push_back(  stod(  strVec[3]) /100    );    //cm to m
    aperture        .push_back(  stod(  strVec[4]) /100    );
    mod             .push_back(  stod(  strVec[5])         );

    for(int j=0;j<8;j++)
    {
        vtemp[j]=stod(strVec.at(j+6));
    }
                 
    factor.push_back(vtemp);
  }

  cout<<"RFQ cell number : "<<UL.size()<<endl;
  cout<<"RFQ length      : "<<cellPositionEnd.back()-cellPosition.front()<<endl;

    rfStream.close();
}

double RFQ::getradius(double z)
{
  int eleNum=	getEleNum(z);
  //cout<<aperture[eleNum]<<endl;
  if(eleNum==-1)  
  {
    return 1e99;
  }
  else
  {
    return aperture[eleNum];
  }  
}

int RFQ::getEleNum(double z)
{
  int count = 0;
    
  while(count < cellPosition.size() && z > cellPositionEnd[count])
  {
    ++count;
  }
  if(count==cellPosition.size())
  {
    count=-1;
    cout<<"Particle in RFQ error: exceed the RFQ length but still in RFQ"<<endl;
    // this actually is not supposed to happen in the simulation
    // if particle is not in RFQ, it will call other Element *ele->Getfield and loss critera.
  }
  return count;
}




