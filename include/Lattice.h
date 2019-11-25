#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include "Global.h"
#include "Element.h"
#include "Beam.h"

 
using std::vector;

class Lattice
{
  public: 
        Lattice();
        ~Lattice();
        int SetElement(vector<string> strVec,Element *ele);
        int GetExternalField(double time, Beam &);
        int GetTransferMatrix(double stepinElement, Beam &);  
        int LostCount(Beam &);
        
        vector<Element * > element; // to store the stucture of the accelerator, every element in vector is one object of the Element class
        vector<double> positionBack;    //the position of the ending for a element

        vector<double> position;        //the position of the beginning for each element
  private:

        vector<double> length;          // list of length for each element
        unsigned int beamSize;
//        Partilce synParticle;

};
#endif
