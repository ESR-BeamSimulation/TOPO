#ifndef LEAPFROG_H
#define LEAPFROG_H

#include "Beam.h"
#include "Pusher.h"


class LeapFrog :public Pusher
{
  public:
    void Update(double timeStep, Beam &beam, Lattice &lattice); 
};

#endif
