#ifndef MESH_H
#define MESH_H

#include <vector>
#include "Beam.h"

using std::vector;
using v1i= vector<int>;
using v1d= vector<double> ;
using v2d= vector<vector<double> > ;
using v3d= vector<vector<vector<double> > > ;

class Mesh 
{
  public:
    Mesh();
    ~Mesh();
    int Initial(vector<int> &numberOfGrid,double particleNum);
    
    vector<int> numberOfGrid;
    int particleNumber;
    
    //which mesh the particle locate
    vector<int> countx;         //=vector<int>(particleNumber);
    vector<int> county;         //=vector<int>(particleNumber); 
    vector<int> countz;         //=vector<int>(particleNumber);
    
    //the coordinate of mesh point
    vector<double> meshx;       //=vector<double>(numberOfGridX);
    vector<double> meshy;       //=vector<double>(numberOfGridY);
    vector<double> meshz;       //=vector<double>(numberOfGridZ);

    //information on mesh point, each is a 3d vector
    v3d rho;            //=v3d(numberOfGridX,   v2d(numberOfGridY,v1d(numberOfGridZ)));
    v3d phi;            //=v3d(numberOfGridX,   v2d(numberOfGridY,v1d(numberOfGridZ)));
    v3d ex;             //=v3d(numberOfGridX,   v2d(numberOfGridY,v1d(numberOfGridZ)));
    v3d ey;             //=v3d(numberOfGridX,   v2d(numberOfGridY,v1d(numberOfGridZ)));
    v3d ez;             //=v3d(numberOfGridX,   v2d(numberOfGridY,v1d(numberOfGridZ)));
    v3d bx;             //=v3d(numberOfGridX,   v2d(numberOfGridY,v1d(numberOfGridZ)));
    v3d by;             //=v3d(numberOfGridX,   v2d(numberOfGridY,v1d(numberOfGridZ)));
    v3d bz;             //=v3d(numberOfGridX,   v2d(numberOfGridY,v1d(numberOfGridZ)));
    int MeshGenerator(Beam &);
};

#endif  //MESH_H
