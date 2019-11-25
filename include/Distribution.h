#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include <vector>
using std::vector;

void KV2DTwiss(const vector<double> &param,vector<vector<double> > &xy);
void KV2DSigma(const vector<double> &param,vector<vector<double> > &xy);
void KV2DAxial(const vector<double> &param,vector<vector<double> > &xy);

void PB2DTwiss(const vector<double> &param,vector<vector<double> > &xy);
void PB2DSigma(const vector<double> &param,vector<vector<double> > &xy);
void PB2DAxial(const vector<double> &param,vector<vector<double> > &xy);

void WB2DTwiss(const vector<double> &param,vector<vector<double> > &xy);
void WB2DSigma(const vector<double> &param,vector<vector<double> > &xy);
void WB2DAxial(const vector<double> &param,vector<vector<double> > &xy);

void GS2DTwiss(const vector<double> &param,vector<vector<double> > &xy);
void GS2DSigma(const vector<double> &param,vector<vector<double> > &xy);
void GS2DAxial(const vector<double> &param,vector<vector<double> > &xy);

void KV3DTwiss(const vector<double> &param,vector<vector<double> > &xy);

vector<double>  TwissToSigma2D(const vector<double> &param);
vector<double>  SigmaToAxial2D(const vector<double> &param);

vector<double>  distribution_base1(const vector<double> &param);
void            distribution_base2(const vector<double> &param,vector<double> &xy);

double Lambert(double x);
#endif
