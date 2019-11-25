#ifndef BEAM_H
#define BEAM_H

#include <vector>
#include <ostream>
#include "Particle.h"
#include "Eigen/Dense"
using namespace Eigen;
using std::vector;
using v1d= vector<double> ;
using v2d= vector<vector<double> > ;



class Beam
{
  public:
    void SetField();
    void SetMatrix();
    void CaculateEmittance();
    void CaculateSigma();
    void Print(ostream &par);
    
    v1d      total           = v1d(6);             
    v1d      average         = v1d(6);              //6 average size of <x> <y>...

    v2d      total2          = v2d(6,v1d(6));
    v2d      average2        = v2d(6,v1d(6));


    v1d      emittance       = v1d(3,0);
    v1d      sigma           = v1d(6,0);            //6 rms sizes average <x^2>, <p_x^2>....
    v1d      alpha           = v1d(3,0);
    v1d      beta            = v1d(3,0);  

    vector<Particle>    particle;
    
//v2d(6,vector(6))
//vector<vector<double> >


    MatrixXd sigmaMatrix = MatrixXd::Zero(6,6);     
    
    //caculate_emittance
    double frq;
    int particleNumber;
    int particleLiveNumber=0;
    int particleLiveInMesh=0;
    double cossyn;
//    vector<double> 	eGamma;
//    vector<double> 	eBeta;
//    vector<double> 	eBetaGamma;
//    vector<double> 	dw;
//    vector<double> 	phi;
//    double eBetaTotal     =0	,eBetaAverage     =0	,eBetaSigma     =0; 
//    double eGammaTotal    =0	,eGammaAverage    =0	,eGammaSigma    =0; 
//    double eBetaGammaTotal=0	,eBetaGammaAverage=0	,eBetaGammaSigma=0;
//    double phiTotal   =0,phiAverage   =0,phiSigma   =0;
//    double PhiGamma=0;
//    double PhiBetaGamma=0;
//    double sigmaz,sigmadz,dwwmax,dppmax,dbbmax;
//    double TBetax,TBetay,TBetaz,TAlphax,TAlphay,TAlphaz,TGammax,TGammay,TGammaz;
//    double emittancex,emittancey,emittancez,correlationx,correlationy,correlationz;
//    double x_sigma,y_sigma,z_sigma,dx_sigma,dy_sigma,dz_sigma,xdx_sigma,ydy_sigma,zdz_sigma;
    
    vector<vector<double> > ParaTransfer;
  private:

};

#endif

