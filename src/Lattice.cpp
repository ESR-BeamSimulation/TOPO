#include "Lattice.h"
#include "Field.h"
#include "Beam.h"
#include<iostream>
#include"Method.h"
#include"MethodMP.h"
#include "Initial.h"
#include <vector>
#include "Pusher.h"
#include "LeapFrog.h"
using namespace std;

Lattice::Lattice()
{
}


Lattice::~Lattice()
{
  //need to release memory which are allcated in InitLattice function
  for(auto a = element.begin();a!=element.end();++a)
  {
//    delete *a;
  }
}

int Lattice::SetElement(vector<string> strVec, Element *ele)
{

  if(position.size()>0)
  {
    position.push_back( position.back() + length.back() );
  }
  else
  {
    position.push_back(0.0);
  }

  length.push_back(stod( strVec[1] ));  //bug fixed 
  double totalLength = position.back() + length.back();
  positionBack.push_back(totalLength);

  ele -> SetParameter(position.back(),strVec);
 
  element.push_back(ele);
}


int Lattice::GetExternalField(double time, Beam & beam)
{
  int particleNumber  = beam.particle.size();

  
  for(unsigned int i=0; i < particleNumber; ++i)
  {
    //1. if the particle is not in this lattice or lost in trensversely
    if(   beam.particle[i].xyz[4] < 0
        ||beam.particle[i].xyz[4] > positionBack.back()
        ||beam.particle[i].transLossFlag==1)
    {
      beam.particle[i].externalElecField = {0,0,0};
      beam.particle[i].externalMegnField = {0,0,0};
      continue;
    }
        

    //2. get which element one particle is in 
    if(beam.particle[i].xyz[4]>position[beam.particle[i].indexElem])
    { 
      while(beam.particle[i].xyz[4] > positionBack[beam.particle[i].indexElem])
      {
        ++beam.particle[i].indexElem;
      }
    }
    else if (beam.particle[i].xyz[4]< position[beam.particle[i].indexElem])
    {
      while(beam.particle[i].xyz[4] < position[beam.particle[i].indexElem] )
      {
        --beam.particle[i].indexElem;
        
        if(beam.particle[i].indexElem<0)
        {
          beam.particle[i].indexElem=0;
          cout<<"Warning! beam.particle[i].indexElem<0"<<endl;
        }
      }
      
    }
//      cout<<beam.particle[i].indexElem<<endl;
//      cout<<"this"<<endl;
//  getchar();
       
    //3. set the field offset, if necessary
//    if(element[beam.particle[i].indexElem]->parameterString[0]=="field")
//    {
//      Field *field = (Field*)element[beam.particle[i].indexElem];

//      if(field->synchronized==0)
//      {
//        field->Synchronize(time,beam);
//      }
//     
//      field=NULL;
//    }
//    
    //4. get the field by calling the GetField function of the corresponding element
    
    element[beam.particle[i].indexElem] -> GetField(time, beam.particle[i], beam.cossyn);


    
//    if(i==0)
//    {
//        cout<<beam.particle[0].xyz[4]<<"\t"<<element[beam.particle[0].indexElem]-> getradius(beam.particle[0].xyz[4])<<endl;
//    }
  }
  
//  5. cout something, if necessary
//    for(int i=0;i<beam.particle.size();++i)
//    {   
//        if(1)
//        {
//            cout<<i<<"  "<<beam.particle[i].indexElem<<"  "<<beam.particle[i].xyz[4]<<"  "<<position[beam.particle[i].indexElem]<<"  "<<positionBack[beam.particle[i].indexElem]<<"  "<<beam.particle[i].externalMegnField[0]<<element[beam.particle[i].indexElem]->parameterString[0]<<endl;
//            getchar();
//        }
//    }
  
  
  return 1;
}

int Lattice::GetTransferMatrix(double stepinElement, Beam &beam)
{

    element[beam.particle[0].indexElem] -> GetMatrix(beam.particle[0]);


}

int Lattice::LostCount(Beam &beam)
{
    //only if the particle hit the transverse pipe can it be lost.
  
  
  double rfperiod   =   C_light/beam.frq * beam.average[5];

    
  for(int i = 0;i<beam.particleNumber;++i)
  {
 
   if(beam.particle[i].transLossFlag )
   {
        continue;
   }

   else if ( abs(beam.particle[i].xyz[4]-beam.average[4])>0.5*rfperiod)
   {
         beam.particle[i].longiLossFlag  =1;
         continue;
   }
    

    if(element[beam.particle[i].indexElem]->shape==0)   //shape=0 means round pipe with aperture apertureX;
    {
      double radius  = sqrt    ( beam.particle[i].xyz[0] * beam.particle[i].xyz[0]
                                +beam.particle[i].xyz[2] * beam.particle[i].xyz[2] );
      double aperture= element[beam.particle[i].indexElem]->apertureX;
      if( radius >= aperture)
      {
//        cout<<radius<<"  "<<aperture<<endl;
        beam.particle[i].transLossFlag  =1;
        beam.particle[i].lossFlag       =1;
      }
    }
    else if(element[beam.particle[i].indexElem]->shape==1)
    {
      if( beam.particle[i].xyz[0] >= element[beam.particle[i].indexElem]->apertureX)
      {
        beam.particle[i].transLossFlag  =1;
        beam.particle[i].lossFlag       =1;
      }
      if( beam.particle[i].xyz[2] >= element[beam.particle[i].indexElem]->apertureY)
      {
        beam.particle[i].transLossFlag  =1;
        beam.particle[i].lossFlag       =1;
      }
    }
  }
}

        
        


    









