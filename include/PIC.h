#ifndef PIC_H
#define PIC_H

#include <vector>
#include <cmath>
#include <fftw3.h>
#include "Beam.h"
#include "Mesh.h" 
#include "InternalField.h"
#include "Global.h"
#include "Lattice.h"
using std::vector;
using v1i= vector<int> ;
using v1d= vector<double> ;
using v2d= vector<vector<double> > ;
using v3d= vector<vector<vector<double> > > ;

class PIC : public InternalField
{
  public:
    ~PIC();
    PIC(Mesh &mesh,Beam &beam);
    PIC();
    int Initial(Mesh &mesh,Beam &beam);                //This mesh is used to store the mesh initialzied by Initial::InitialMesh()

    void GetInternalFiled(Mesh &mesh, Beam &beam, Lattice &lattice)
    {
      Weight(mesh,beam,lattice);
      FFT(mesh);
      Allocate(mesh,beam);
    }

  protected:
    //some default numbers

    //Property of grid
    int numofgridx;
    int numofgridy;
    int numofgridz;
    double stepx,stepy,stepz,vt;       //vt=stepx*stepy*stepz;//unit:m

    double particleInMeshNumber;
    v2d weighOfGrid;         //store assigned volume for each particle in certain mesh 
    v1i particleOutOfMesh;

  private:
    int Weight(Mesh &mesh, Beam &beam, Lattice &lattice);
    int FFT(Mesh &mesh);
    int Allocate(const Mesh &mesh, Beam &beam);
    
    //fftw
    double              *out;
    double              *in,            *outf;  
    double              *in_xy,         *out_xy;
 
    fftw_complex        *in_z,          *out_z,         *outf_z;
    fftw_plan           rhoGridToFourier_XY,    phiFourierToGrid_XY;
    fftw_plan           rhoGridToFourier_Z,     phiFourierToGrid_Z;
    fftw_r2r_kind       kind_xy_forward[2] = {FFTW_RODFT00,FFTW_RODFT00};
    fftw_r2r_kind       kind_xy_backward[2]= {FFTW_RODFT00,FFTW_RODFT00};

};
#endif
