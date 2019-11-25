#include "MyFunc.h"
#include <fstream>
#include <iostream>
#include <iomanip>

#include <iterator>
#include <cstdio> 
#include <cmath>

using namespace std;

void StringSplit(const string &s,vector<string> &vec)
{
  if(vec.size()>0)
  vec.clear();
  int length = s.length();
  int start=0;
  char splitchar=' ';
  char splitchar2 = '	';
  for(unsigned int i=0;i<length;i++)
  {
    if (s[i]=='!')
    {
      i=length;
    }
    else if((s[i] == splitchar || s[i] == splitchar2) && i==start)
    {
      start++;
    }
    else if(s[i] == splitchar || s[i] == splitchar2)
    {
      vec.push_back(s.substr(start,i - start));
      start = i+1;
    }
    else if(i == length-1)
    {
      vec.push_back(s.substr(start,i+1 - start));
    }
  }
}

int GetLineNumber(const string &p)
{
  ifstream file(p);
  string str;
  vector<string> strVec;
  int count = 0;
  while (file.peek()!=EOF)
  {
    getline(file, str);
    StringSplit(str,strVec);
    if (strVec.size() > 0) 
    {
      count ++;
    }
  }
  file.close();
  return count;
}

double Bessi0(double x)
{
  double result=0;
  double nx=1;
  for(int i=0;i<20;++i)
  {
    if(i>0)
    {
      nx/=i;
    }
    result+=pow(x/2,2*i)/(tgamma(i+1))*nx;
  }
  return result;
}

double Bessi1(double x)
{
  double result=0;
  double nx=1;
  for(int i=0;i<20;++i)
  {
    if(i>0)
    {
      nx/=i;
    }
    result+=pow(x/2,2*i+1)/(tgamma(i+1+1))*nx;
  }
  return result;
}

double Bessin(int n,double x)
{
  double result=0;
  double nx=1;
  for(int i=0;i<20;++i)
  {
    if(i>0)
    {
      nx/=i;
    }
    result+=pow(x/2,2*i+n)/(tgamma(i+1+n))*nx;
  }
  return result;
}

void PrintV1(ostream &par,const vector<double>  &pv)
{
  ostream_iterator<double> outpar(par,"\t");
  copy(pv.begin(),pv.end(),outpar);
  par<<endl;
}

void PrintV2(ostream &par,const vector<vector<double>  > &pv)
{
  ostream_iterator<double> outpar(par,"\t");
  //for_each(pv.begin(),pv.end(), [&par,&outpar](const vector<double> &a){  copy(a.begin(),a.end(),outpar);par<<endl;});
  for(const auto &a:pv) 
  {
    copy(a.begin(),a.end(),outpar);
    par<<endl;
  }
}

