#include"LinearMT.h"
#include<vector>
#include<iostream>
#include"Eigen/Dense"
#include"Lattice.h"

using namespace Eigen;
using namespace std;



void LinearMT::Update(Beam &beam, Element *elem)
{
    
    MatrixXd transpos  = MatrixXd::Zero(6,6);
    MatrixXd transfer  = MatrixXd::Zero(6,6);  
    MatrixXd sigmaMat0 = MatrixXd::Zero(6,6);
    MatrixXd sigmaMat1 = MatrixXd::Zero(6,6); 
    
    transfer    =   elem->transferMatrixM;
    transpos    =   transfer.transpose();
    
    sigmaMat0   =   beam.sigmaMatrix;
    
    sigmaMat1   =   transfer * sigmaMat0 * transpos;
    
    beam.sigmaMatrix    =   sigmaMat1;
    
    return;
}

