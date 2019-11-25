#ifndef METHODMT_H
#define METHODMT_H

#include <iostream>
#include <vector>
#include "Method.h"
#include "Beam.h"
#include "Lattice.h"


class MethodMT: public Method
{
public:
    void Push(Beam &, Lattice &);

private:
    double stepinElement;
      
      //functions for Matrix Munipulation  
     

};



#endif
