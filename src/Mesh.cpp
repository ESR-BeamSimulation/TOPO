#include "Mesh.h"
#include <iostream>
using namespace std;

Mesh::Mesh()
{
}

Mesh::~Mesh()
{

}

int Mesh::Initial(vector<int> &numberGrid,double particleNum)
{
  numberOfGrid  = numberGrid;
  particleNumber= particleNum;

  countx        = v1i(particleNumber);
  county        = v1i(particleNumber); 
  countz        = v1i(particleNumber);
  meshx         = v1d(numberOfGrid[0]);
  meshy         = v1d(numberOfGrid[1]);
  meshz         = v1d(numberOfGrid[2]);

  rho           = v3d(numberOfGrid[0],   v2d(numberOfGrid[1],v1d(numberOfGrid[2])));
  phi           = v3d(numberOfGrid[0],   v2d(numberOfGrid[1],v1d(numberOfGrid[2])));
  ex            = v3d(numberOfGrid[0],   v2d(numberOfGrid[1],v1d(numberOfGrid[2])));
  ey            = v3d(numberOfGrid[0],   v2d(numberOfGrid[1],v1d(numberOfGrid[2])));
  ez            = v3d(numberOfGrid[0],   v2d(numberOfGrid[1],v1d(numberOfGrid[2])));
  bx            = v3d(numberOfGrid[0],   v2d(numberOfGrid[1],v1d(numberOfGrid[2])));
  by            = v3d(numberOfGrid[0],   v2d(numberOfGrid[1],v1d(numberOfGrid[2])));
  bz            = v3d(numberOfGrid[0],   v2d(numberOfGrid[1],v1d(numberOfGrid[2])));

}
int Mesh::MeshGenerator(Beam & beam)
{

}
