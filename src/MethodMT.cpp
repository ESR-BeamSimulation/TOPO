#include <iostream>
#include "MethodMT.h"
#include"LinearMT.h"
#include "Eigen/Dense"
#include "Lattice.h"
using namespace std;







void MethodMT::Push(Beam &beam,Lattice &lattice)
{
    LinearMT    linearMT;
    
    
    
    for(int i=0;i<lattice.element.size();++i)
    {
        beam.SetMatrix();
        lattice.element[i]->GetMatrix(beam.particle[0]);
        linearMT.   Update(beam,lattice.element[i]);
    
    }

}
