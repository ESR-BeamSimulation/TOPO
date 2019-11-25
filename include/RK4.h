#ifndef RK4_H
#define RK4_H

#include "Beam.h"
#include "Pusher.h"
#include <cmath>

class RK4 :public Pusher
{
  public:
    void Update(double timeStep, Beam &beam); 
};

#endif
