#ifndef PARTICLE_H
#define PARTICLE_H
#include "Eigen/Dense"
#include <vector>
using namespace std;
using namespace Eigen;

class Particle
{
  public:
    Particle();
    vector<double>      xyz;

    vector<double>      externalElecField;
    vector<double>      externalMegnField;

    vector<double>      internalElecField;
    vector<double>      internalMegnField;

    vector<double>      totalElecField;
    vector<double>      totalMegnField;

    // variables for transfer matrix of each particle
    MatrixXd particleTranferMap;

    int     indexElem;                       //defuldt vlue is zero,  inilize is not nessceary
    int    outSCMeshFlag;
    int    transLossFlag;                   //1 means loss
    int    longiLossFlag;                   //1 means loss
    int         lossFlag;                   //1 means loss
    int    buncherNmber;                   //to know which bunch the particle are located.
    double charge;
    double mass;
    double qOverMass;

    void GetTotalField();               // to obtain the total field
    void GetMatrix();                   // to obtain the transfer matrix for each particle
};


#endif



