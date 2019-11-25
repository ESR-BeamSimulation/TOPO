#include <iostream>
#include <string>
#include <vector>

#include "Beam.h"
#include "Lattice.h"
#include "Initial.h"
#include "MethodMP.h"
#include "MethodMT.h"

using namespace std;

int main()
{
  Initial               initialize;
  Beam                  beam;
  Lattice               lattice;
  initialize.InitLattice(lattice);
  initialize.SetSynPhase(lattice);  


  initialize.InitParticle(beam);
 
  
  
  cout<<beam.particle.size()<<endl; 
  beam.         CaculateEmittance     (               ); // to obtain the initial parameter before revolution
  cout<<beam.particle.size()<<endl; 




  if( initialize.GetMethod()=="MP")
  {
    Mesh mesh;
    initialize.InitMesh(mesh);

    
    MethodMP    MPrun;
    initialize.InitMethodMP(MPrun);
                       

    if (initialize.GetInternalFieldMethodethod()=="PIC")
    {
      MPrun.Push(lattice,mesh,beam);
    }
    else if (initialize.GetInternalFieldMethodethod()=="PTP")
    {
      MPrun.Push(lattice,mesh,beam);
    }
  }
  
  else if( initialize.GetMethod()=="MT")
  {
    MethodMT    MTrun;
    MTrun.Push(beam,lattice);
  }
  cout<<"Complete!"<<endl;
  return 0;
  
//  for(int i=0;i<beam.particle.size()&&false;++i)
//  {
//    for(int j=0;j<6;++j)
//     cout<<beam.particle.at(i).xyz[j]<<"  ";
//    cout<<endl;
//  }
}

