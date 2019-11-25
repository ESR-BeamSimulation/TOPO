#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>
#include <omp.h>

#include "MyFunc.h"
#include "PIC.h"
using v1d= vector<double> ;
using v2d= vector<vector<double> > ;
using v3d= vector<vector<vector<double> > > ;

using namespace std;

PIC::PIC()
{
}
PIC::PIC(Mesh &mesh,Beam &beam)
{
  Initial(mesh,beam);
}
PIC::~PIC()
{
  fftw_free(in);
  fftw_free(out);

  fftw_free(in_xy);
  fftw_free(out_xy);

  fftw_free(in_z);
  fftw_free(out_z);

  fftw_free(outf); 
  fftw_free(outf_z);
 

  fftw_destroy_plan(rhoGridToFourier_XY);
  fftw_destroy_plan(phiFourierToGrid_XY);

  fftw_destroy_plan(rhoGridToFourier_Z);
  fftw_destroy_plan(phiFourierToGrid_Z);
 
}

int PIC::Initial(Mesh &mesh,Beam &beam)
{
  cout<<"PIC Initialing..."<<endl;

  numofgridx = mesh.numberOfGrid[0];
  numofgridy = mesh.numberOfGrid[1];
  numofgridz = mesh.numberOfGrid[2];


  weighOfGrid       =  v2d(beam.particle.size(), v1d(8));   //store assigned volume for each particle in certain mesh
  particleOutOfMesh =  v1i(beam.particle.size(), 0);

  in  	= (double*)             fftw_malloc(sizeof(double)      * numofgridx * numofgridy * numofgridz);
  out 	= (double*)             fftw_malloc(sizeof(double)      * numofgridx * numofgridy * numofgridz);
  outf 	= (double*)             fftw_malloc(sizeof(double)      * numofgridx * numofgridy * numofgridz);
  in_xy = (double*)             fftw_malloc(sizeof(double)      * numofgridx * numofgridy * numofgridz);
  out_xy= (double*)             fftw_malloc(sizeof(double)      * numofgridx * numofgridy * numofgridz);
  in_z  = (fftw_complex*)       fftw_malloc(sizeof(fftw_complex)* numofgridx * numofgridy * numofgridz);
  out_z = (fftw_complex*)       fftw_malloc(sizeof(fftw_complex)* numofgridx * numofgridy * numofgridz);
  outf_z= (fftw_complex*)       fftw_malloc(sizeof(fftw_complex)* numofgridx * numofgridy * numofgridz);


 



  //1
  fftw_plan_with_nthreads(omp_get_max_threads());
  int n_xy[]={numofgridx,numofgridy};

  rhoGridToFourier_XY=fftw_plan_many_r2r(2, n_xy, numofgridz,

      in_xy  ,n_xy,
      numofgridz, 1,
      out_xy ,n_xy,
      numofgridz, 1,
      kind_xy_forward, FFTW_MEASURE);
  //2
  fftw_plan_with_nthreads(omp_get_max_threads());
  int n_z[]={numofgridz};
  rhoGridToFourier_Z=fftw_plan_many_dft(1, n_z, numofgridx*numofgridy,
      in_z , n_z ,
      1,numofgridz,
      out_z, n_z ,
      1,numofgridz,
      FFTW_FORWARD, FFTW_MEASURE);
  //3
  fftw_plan_with_nthreads(omp_get_max_threads());
  phiFourierToGrid_Z=fftw_plan_many_dft(1, n_z, numofgridx*numofgridy,
      out_z  , n_z ,
      1,numofgridz,
      outf_z, n_z ,
      1,numofgridz,
      FFTW_BACKWARD, FFTW_MEASURE);
  //4
  fftw_plan_with_nthreads(omp_get_max_threads());
  phiFourierToGrid_XY=fftw_plan_many_r2r(2, n_xy, numofgridz,
      out_xy,n_xy,
      numofgridz, 1,
      in_xy,n_xy ,
      numofgridz, 1,
      kind_xy_backward, FFTW_MEASURE);
      
  return 1;
}

int PIC::Weight(Mesh &mesh,Beam &beam, Lattice &lattice)
{
    if(beam.particle.size()==0)
    {
      cout<<"No particle exist!"<<endl;
      return -1;
    }

  for(int i = 0; i < mesh.numberOfGrid[0]; ++i)
  {
    for(int j = 0; j < mesh.numberOfGrid[1]; ++j)
    {
      for(int k = 0; k < mesh.numberOfGrid[2]; ++k)
      {
        mesh.rho[i][j][k]=0;
      }
    }
  }
  
  for(int i=0;i<numofgridx;i++)
  {
    mesh.meshx[i]   = 0;
  }
  for(int i=0;i<numofgridy;++i)
  {
    mesh.meshy[i]   =0;
  }
  for(int i=0;i<numofgridz;++i)
  {
    mesh.meshz[i]   =0;
  }
  
 for (int i=0; i<beam.particle.size();++i)
 {
      particleOutOfMesh[i]=0;              
      beam.particle[i].outSCMeshFlag=0;
      beam.particle[i].longiLossFlag=0;
 } 
  

  double rfperiod   =  1.1 * C_light/beam.frq * beam.average[5];

  stepx =   beam.sigma[0]*2*5/(numofgridx-1);
  stepy =   beam.sigma[2]*2*5/(numofgridy-1);
  stepz =   rfperiod/(numofgridz-1);
  


  if(0)//&&.name=="rfq")
  {
    stepx=  2*lattice.element[beam.particle[0].indexElem]->apertureX / (numofgridx-1);
    stepy=  2*lattice.element[beam.particle[0].indexElem]->apertureY / (numofgridx-1);
    stepz = 2* beam.sigma[4]*2*2*2/(numofgridz-1);
  }
  
   
  for(int i=0;i<numofgridx;i++)
  {
    mesh.meshx[i]=beam.average[0]+stepx*(i - numofgridx/2.0+0.5);
  }
  for(int i=0;i<numofgridy;++i)
  {
    mesh.meshy[i]=beam.average[2]+stepy*(i - numofgridy/2.0+0.5);
  }
  for(int i=0;i<numofgridz;++i)
  {
    mesh.meshz[i]=beam.average[4]+stepz*(i  -numofgridz/2.0+0.5);
  }
  vt=stepx*stepy*stepz;
        
 
  int idx,idy,idz;
  
  for(int i=0;i<beam.particle.size();i++)
  {
    if(beam.particle[i].lossFlag)   // if particle is lost transversely, do not deal with PIC
    {
      continue;
    }
    
    mesh.countx[i]=floor((beam.particle[i].xyz[0]-mesh.meshx[0])/stepx);
    mesh.county[i]=floor((beam.particle[i].xyz[2]-mesh.meshy[0])/stepy);
    mesh.countz[i]=floor((beam.particle[i].xyz[4]-mesh.meshz[0])/stepz);

    if(   mesh.countx[i]>=numofgridx-1         //Take a notice about FFT and FFT_sin in the mannul of FFTW lib.       
        ||mesh.county[i]>=numofgridy-1           
        ||mesh.countz[i]>=numofgridz          
        ||mesh.countx[i]<0                              
        ||mesh.county[i]<0
        ||mesh.countz[i]<0)
    {
      particleOutOfMesh[i]=1;              
      beam.particle[i].outSCMeshFlag=1;
      for(int j=0;j<8;++j)
      {
        weighOfGrid[i][j]=  0;
      }
    }
    else
    {
      particleOutOfMesh[i]=0;
      beam.particle[i].outSCMeshFlag=0;
      for(int j=0;j<8;++j)
      {
        if(mesh.countz[i]==numofgridz-1)
        {
           idz=(j <=3)     ? -mesh.countz[i] : 0 ;
        }
        else
        {
            idz=(j <=3) 	? 1 : 0;
        }

        weighOfGrid[i][j] =  
        abs( beam.particle[i].xyz[0] - mesh.meshx[mesh.countx[i]+ ((j%2==0)	? 1 : 0) ]) * 
        abs( beam.particle[i].xyz[2] - mesh.meshy[mesh.county[i]+ ((j%4< 2) ? 1 : 0) ]) *
        abs( beam.particle[i].xyz[4] - mesh.meshz[mesh.countz[i]+ idz ])
        /vt;
      }
    }
  }

  particleInMeshNumber=0;
  for(int i=0;i<beam.particle.size();++i)
  {
    if(beam.particle[i].lossFlag or particleOutOfMesh[i]) continue;
    
    ++particleInMeshNumber;
    mesh.rho[mesh.countx[i]  ][mesh.county[i]  ][mesh.countz[i]]      += weighOfGrid[i][0] * beam.particle[i].charge;
    mesh.rho[mesh.countx[i]+1][mesh.county[i]  ][mesh.countz[i]]      += weighOfGrid[i][1] * beam.particle[i].charge;
    mesh.rho[mesh.countx[i]  ][mesh.county[i]+1][mesh.countz[i]]      += weighOfGrid[i][2] * beam.particle[i].charge;
    mesh.rho[mesh.countx[i]+1][mesh.county[i]+1][mesh.countz[i]]      += weighOfGrid[i][3] * beam.particle[i].charge;
    if(mesh.countz[i]==numofgridz-1)
    {
      mesh.rho[mesh.countx[i]  ][mesh.county[i]  ][0]   += weighOfGrid[i][4]*beam.particle[i].charge;
      mesh.rho[mesh.countx[i]+1][mesh.county[i]  ][0]   += weighOfGrid[i][5]*beam.particle[i].charge;
      mesh.rho[mesh.countx[i]  ][mesh.county[i]+1][0]   += weighOfGrid[i][6]*beam.particle[i].charge;
      mesh.rho[mesh.countx[i]+1][mesh.county[i]+1][0]   += weighOfGrid[i][7]*beam.particle[i].charge;
    }
    else
    {
      mesh.rho[mesh.countx[i]  ][mesh.county[i]  ][mesh.countz[i]+1]    += weighOfGrid[i][4]*beam.particle[i].charge;   //unit is coumb
      mesh.rho[mesh.countx[i]+1][mesh.county[i]  ][mesh.countz[i]+1]    += weighOfGrid[i][5]*beam.particle[i].charge;
      mesh.rho[mesh.countx[i]  ][mesh.county[i]+1][mesh.countz[i]+1]    += weighOfGrid[i][6]*beam.particle[i].charge;
      mesh.rho[mesh.countx[i]+1][mesh.county[i]+1][mesh.countz[i]+1]    += weighOfGrid[i][7]*beam.particle[i].charge;
    }
  }
      if (particleInMeshNumber<int (0.1*beam.particle.size()))
      {
        cout<<"the SC solver is not right, more than 90/% particles are out of SC mesh"<<endl;
        exit(0);
      }
      
//      cout<<"how many particle are still in mesh  "<<particleInMeshNumber<<endl;
//    cout<<"the electric field particle feel     ";
//    for (int i=0; i<3;++i)
//    cout<<beam.particle[10].internalElecField[i]<<"\t";
//    cout<<beam.particle[10].transLossFlag<<"\t"<<particleOutOfMesh[10]<<endl;


  return 1;
}



int PIC::FFT(Mesh &mesh)
{


  //#pragma omp for
  for(int i=0;i<numofgridx;++i)
  {
    for(int j=0;j<numofgridy;++j)
    {
      for(int k=0;k<numofgridz;++k)
      {
//       mesh.rho[i][j][k]=0;
        in_xy[i*numofgridy*numofgridz+j*numofgridz+k]=mesh.rho[i][j][k]/stepx/stepy/stepz;   // c/m^3
      }
    }
  }
  
//  
//  mesh.rho[32][32][32]=1.0;
//  in_xy   [32*numofgridy*numofgridz+32*numofgridz+32]=1.0;

 
  //#pragma omp barrier
  //#pragma omp single 
  fftw_execute(rhoGridToFourier_XY);                //  in_xy->out_xy
  //#pragma omp barrier
    

  int coun365;
  //#pragma omp for
  for(int i=0;i<numofgridx;++i)
  {
    for(int j=0;j<numofgridy;++j)
    {
      for(int k=0;k<numofgridz;++k)
      {
        coun365=i*numofgridy*numofgridz+j*numofgridz+k;
        in_z[coun365][0]=out_xy[coun365];          //   out_xy-> in_z
        in_z[coun365][1]=0;
      }
    }
  }

  //#pragma omp barrier
  //#pragma omp single 
  fftw_execute(rhoGridToFourier_Z);                 //  in_z->out_z
  //#pragma omp barrier
 

  double K_rho2phi=0,knx=0,kny=0,knz=0;
  int coun445;
  //#pragma omp for
  for(int i=0;i<numofgridx;++i)
  {
    for(int j=0;j<numofgridy;++j)
    {
      for(int k=0;k<numofgridz;++k)
      {

        knx=M_PI/2*(i+1)/(numofgridx+1);
        kny=M_PI/2*(j+1)/(numofgridy+1);
        knz=M_PI  *(k)  /(numofgridz);          // chech the formular used in FFTW mannul in function FFT,  
                                                // The SK function is uniquely decided/////
        K_rho2phi=   pow(2*sin(knx)/stepx,2)
                    +pow(2*sin(kny)/stepy,2)
                    +pow(2*sin(knz)/stepz,2);   // unit is 1/m^2
        
        coun445=i*numofgridy*numofgridz+j*numofgridz+k;
        out_z[coun445][0] = out_z[coun445][0] /K_rho2phi /Epsilon;   // unit is  (C/m^3) / (1/m^2) /(C/m*V)==V                
        out_z[coun445][1] = out_z[coun445][1] /K_rho2phi /Epsilon;
      }
    }
  }
  
  //#pragma omp barrier
  //#pragma omp single 
  fftw_execute(phiFourierToGrid_Z);             // out_z  -> outf_z
  //#pragma omp barrier
 
  //#pragma omp for
  for(int i=0;i<numofgridx;++i)
  {
    for(int j=0;j<numofgridy;++j)
    {
      for(int k=0;k<numofgridz;++k)
      {
            coun365=i*numofgridy*numofgridz+j*numofgridz+k;
            out_xy[coun365]=outf_z[coun365][0];                     // outf_z -> out_xy
      }
    }
  }
  
  //#pragma omp barrier
  //#pragma omp single 
  fftw_execute(phiFourierToGrid_XY);                                //  out_xy -> in_xy
  //#pragma omp barrier

  //#pragma omp for
  for(int i=0;i<numofgridx;i++)
  {
    for(int j=0;j<numofgridy;j++)
    {
      for(int k=0;k<numofgridz;++k)
      {
        mesh.phi[i][j][k] = in_xy[i*numofgridy*numofgridz+j*numofgridz+k]
                          / (4*(numofgridx+1)*(numofgridy+1)*(numofgridz));
      }
    }
  }



  //Diferentiate to get the efield
  //#pragma omp for
  for(int i=0;i<numofgridx;++i)
  {
    for(int j=0;j<numofgridy;++j)
    {
      for(int k=0;k<numofgridz;++k)
      {
        if(i==0)
        {
           mesh.ex[i][j][k]        =    -(mesh.phi[i+1][j][k]  - 0)     /(2*stepx);                  //unit:V/m;
        }
        else if(i==numofgridx-1)
        {
           mesh.ex[i][j][k]        =    -(0  -mesh.phi[i-1][j][k])      /(2*stepx);                  //unit:V/m;
        }
        else
        {
           mesh.ex[i][j][k]        =    -(mesh.phi[i+1][j][k]  -mesh.phi[i-1][j][k])    /(2*stepx);
        }
        
        if(j==0)
        {
           mesh.ey[i][j][k]        =    -(mesh.phi[i][j+1][k]  -0)                      /(2*stepy);
        }
        else if(j==numofgridy-1)
        {
           mesh.ey[i][j][k]        =    -(0  -mesh.phi[i][j-1][k])                      /(2*stepy);
        }
        else
        {
           mesh.ey[i][j][k]        =    -(mesh.phi[i][j+1][k]  -mesh.phi[i][j-1][k])    /(2*stepy);
        }
        
        if(k==0)
        {
          mesh.ez[i][j][k]         =   -(mesh.phi[i][j][k+1]  -mesh.phi[i][j][numofgridz-1])/(2*stepz);
        }
        else if(k==numofgridz-1)
        {
          mesh.ez[i][j][k]         =   -(mesh.phi[i][j][0]    -mesh.phi[i][j][k-1])      /(2*stepz);
        }
        else
        {
          mesh.ez[i][j][k]         =   -(mesh.phi[i][j][k+1]  -mesh.phi[i][j][k-1])      /(2*stepz);
        }
      }
    }
  }

  return 0;
}

int PIC::Allocate(const Mesh &mesh, Beam &beam)
{
  double e1=0,e2=0,e3=0,b1=0,b2=0,b3=0;
  double gamma,beta;
  int nParticleInMesh=0;

  for(int i=0;i<beam.particle.size();++i)
  {
    if(beam.particle[i].transLossFlag||beam.particle[i].outSCMeshFlag||beam.particle[i].longiLossFlag)           //particle are lost transversely or out of space chare mesh  
    {                                                                                           //
      beam.particle[i].internalElecField[0]=0;
      beam.particle[i].internalElecField[1]=0;
      beam.particle[i].internalElecField[2]=0;
      beam.particle[i].internalMegnField[0]=0;
      beam.particle[i].internalMegnField[1]=0;
      beam.particle[i].internalMegnField[2]=0;
      continue;
    }
    else                                                                                        // only particle are inside of the space charge mesh are calculated
    {
        ++ nParticleInMesh;
      e1  =   mesh.ex[mesh.countx[i]  ][mesh.county[i]  ][mesh.countz[i]  ]     * weighOfGrid[i][0]
	+     mesh.ex[mesh.countx[i]+1][mesh.county[i]  ][mesh.countz[i]  ]     * weighOfGrid[i][1]
	+     mesh.ex[mesh.countx[i]  ][mesh.county[i]+1][mesh.countz[i]  ]     * weighOfGrid[i][2]
	+     mesh.ex[mesh.countx[i]+1][mesh.county[i]+1][mesh.countz[i]  ]     * weighOfGrid[i][3]
	+     mesh.ex[mesh.countx[i]  ][mesh.county[i]  ][mesh.countz[i]+1]     * weighOfGrid[i][4]
	+     mesh.ex[mesh.countx[i]+1][mesh.county[i]  ][mesh.countz[i]+1]    	* weighOfGrid[i][5]
	+     mesh.ex[mesh.countx[i]  ][mesh.county[i]+1][mesh.countz[i]+1]    	* weighOfGrid[i][6]
	+     mesh.ex[mesh.countx[i]+1][mesh.county[i]+1][mesh.countz[i]+1]     * weighOfGrid[i][7];
      e2  =   mesh.ey[mesh.countx[i]  ][mesh.county[i]  ][mesh.countz[i]  ]     * weighOfGrid[i][0]
	+     mesh.ey[mesh.countx[i]+1][mesh.county[i]  ][mesh.countz[i]  ]    	* weighOfGrid[i][1]
	+     mesh.ey[mesh.countx[i]  ][mesh.county[i]+1][mesh.countz[i]  ]    	* weighOfGrid[i][2]
	+     mesh.ey[mesh.countx[i]+1][mesh.county[i]+1][mesh.countz[i]  ]     * weighOfGrid[i][3]
	+     mesh.ey[mesh.countx[i]  ][mesh.county[i]  ][mesh.countz[i]+1]     * weighOfGrid[i][4]
	+     mesh.ey[mesh.countx[i]+1][mesh.county[i]  ][mesh.countz[i]+1]     * weighOfGrid[i][5]
	+     mesh.ey[mesh.countx[i]  ][mesh.county[i]+1][mesh.countz[i]+1]     * weighOfGrid[i][6]
	+     mesh.ey[mesh.countx[i]+1][mesh.county[i]+1][mesh.countz[i]+1]     * weighOfGrid[i][7];
      e3  =   mesh.ez[mesh.countx[i]  ][mesh.county[i]  ][mesh.countz[i]  ]     * weighOfGrid[i][0]
	+     mesh.ez[mesh.countx[i]+1][mesh.county[i]  ][mesh.countz[i]  ]    	* weighOfGrid[i][1]
	+     mesh.ez[mesh.countx[i]]  [mesh.county[i]+1][mesh.countz[i]  ]    	* weighOfGrid[i][2]
	+     mesh.ez[mesh.countx[i]+1][mesh.county[i]+1][mesh.countz[i]  ]     * weighOfGrid[i][3]
	+     mesh.ez[mesh.countx[i]  ][mesh.county[i]  ][mesh.countz[i]+1]     * weighOfGrid[i][4]
	+     mesh.ez[mesh.countx[i]+1][mesh.county[i]  ][mesh.countz[i]+1]     * weighOfGrid[i][5]
	+     mesh.ez[mesh.countx[i]  ][mesh.county[i]+1][mesh.countz[i]+1]     * weighOfGrid[i][6]
	+     mesh.ez[mesh.countx[i]+1][mesh.county[i]+1][mesh.countz[i]+1]     * weighOfGrid[i][7];


      //Lorentz Transform verir
      gamma  =  sqrt(1 + pow(beam.particle[i].xyz[1],2) +pow(beam.particle[i].xyz[3],2) +pow(beam.particle[i].xyz[5],2));
      beta      = sqrt(gamma*gamma-1)/gamma;
      
//      cout<< beta<<" "<<gamma<<endl;
//      e1        =e1/C_light;
//      e2        =e2/C_light;
//      e3        =e3/C_light;
//      beam.particle[i].internalElecField[0]= C_light    *gamma*(e1+beta*b2);
//      beam.particle[i].internalElecField[1]= C_light    *gamma*(e2-beta*b1);
//      beam.particle[i].internalElecField[2]= C_light    *e3;
//      beam.particle[i].internalMegnField[0]=            gamma*(b1-beta*e2);
//      beam.particle[i].internalMegnField[1]=            gamma*(b2+beta*e1);
//      beam.particle[i].internalMegnField[2]=            b3;
      

      beam.particle[i].internalElecField[0]= e1 / gamma;
      beam.particle[i].internalElecField[1]= e2 / gamma;
      beam.particle[i].internalElecField[2]= e3 / gamma;
      beam.particle[i].internalMegnField[0]= 0;
      beam.particle[i].internalMegnField[1]= 0;
      beam.particle[i].internalMegnField[2]= 0;
      
    }

  }

}   

