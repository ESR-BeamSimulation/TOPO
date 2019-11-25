#include "FieldSB.h"
#include "MyFunc.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <Global.h>

using namespace std;

FieldSB::~FieldSB()
{
}

FieldSB::FieldSB()
{
}

void FieldSB::SetParameter(double posit, const vector<string> & param)
{
    if(param.size()<7)
    {
        cout<<"Error! parameter not enough in element : \n"<<endl;
        for(int i = 0;i<param.size();++i)
        {
            cout<<param.at(i)<<"\t";
        }
        cout<<endl;
        exit(0);
    }
  
    SetBasicParameter(posit,param);
  
    locationEnd   =       0;

    locationEnd   =       position+length;
    type          =       stoi(param[4]);
    kb            =       stod(param[5]);
    filename      =       param[6];
      
    ifstream fin;
    string                str;
    vector<string>        strVec;
      
    vector<string> filedex={".bsz",".bsr"};
      
    for(int dim=0;dim<filedex.size();dim++)
    {
        string filepath=filename+filedex[dim];
          
        fin.open(filepath);
          
        if(!fin.is_open())
        {
            cout<<"field type2 file error : "<<filepath<<endl;
            exit(0);
        }
          
        if(dim==0)
        {
            fin>>nz>>zmax;
            fin>>nr>>rmax;
            field.resize((nz+1)*(nr+1),vector<double>(2));
            locationMax=position+zmax;
            spacer=(rmax     )/nr;
            spacez=(zmax     )/nz;
            getline(fin,str);
            getline(fin,str);
        }
        else
        {
            for(int i = 0;i<3;++i)      
             {
              //to skip several lines
              getline(fin,str);
             }
        }
        if(stoi(str)!=1)
        {
            cout<<"mark error : "<<filepath<<"  "<<str<<endl;
        }
          
        double temp=0;
        int index=0;
        for(int j=0;j<=nz;++j)
        {
            for(int i=0;i<=nr;++i)
            {
              fin>>temp;
              index=j*(nr+1)+i;
              field[index][dim]=temp*kb;
            }
        }
        
        fin.close();
        gridNum=field.size();
      }
      
      string fout=filename+"_fieldout"+".dat";
      ofstream fieldout(fout);
      
      
    for(int j=0;j<=nz;++j)
    {
      for(int i=0;i<=nr;++i)
      {
        fieldout<<i*spacer<<"  ";
        fieldout<<j*spacez<<"  ";
        double index=j*(nr+1)+i;
        for(int dim=0;dim<2;dim++)
        {
          fieldout<<field[index][dim]<<"  ";
        }
        fieldout<<endl;
      }
    }

}

void FieldSB::GetField(double time,Particle &particle,double synPhaseoffse)
{
    double r    =sqrt(particle.xyz[0]*particle.xyz[0]+particle.xyz[2]*particle.xyz[2]);
    double theta=atan2(particle.xyz[2],particle.xyz[0]);
    if( r>rmax
        ||particle.xyz[4]<position
        ||particle.xyz[4]>locationEnd)
    {
      particle.externalElecField=zero3d;
      particle.externalMegnField=zero3d;
      return;
    }
    else if(particle.xyz[4]>locationMax)
    {
      particle.externalElecField=zero3d;
      particle.externalMegnField=zero3d;
      return;
    }

    double rn,zn;
    int rGrid,zGrid,indexGrid;
    v1d v(4);
    v2d fieldOnGrid(2,v1d(4));
    v1d fieldzr(2,0.0); 
    rn=r                /spacer;
    zn=(particle.xyz[4]-position)/spacez;

    if(r==0)  // when particle are loacted on axis
    {
          
        rGrid =   floor(rn);
        zGrid =   floor(zn);
          
        indexGrid=zGrid*(nr+1)+rGrid;
        
        fieldzr[1]  =   0;
        fieldzr[0]  =   abs(zn-floor(zn)-1) * field[indexGrid][0] +  abs(zn-ceil(zn)) *  field[(zGrid+1)*(nr+1)+rGrid][0];
        
    }
    else
    {   
        for(int i=0;i<4;++i)
        {
            rGrid =   ((i%2< 1)?ceil(rn):ceil(rn)-1);
            zGrid =   ((i  < 2)?ceil(zn):ceil(zn)-1);

          indexGrid=zGrid*(nr+1)+rGrid;
          
          
          if(ceil(rn)==rn)
          {
             rGrid +=1;
          }
          if(ceil(zn)==zn)
          {
             zGrid +=1;
          }      

          v[i]  =   abs(rn-rGrid) * abs(zn-zGrid);


          for(int dim=0;dim<2;++dim)
          {
            fieldOnGrid[dim][i]=v[i]*field[indexGrid][dim];
            fieldzr[dim]+=fieldOnGrid[dim][i];
          }
        }
    }
    particle.externalMegnField.at(0)=fieldzr[1]*cos(theta);
    particle.externalMegnField.at(1)=fieldzr[1]*sin(theta);
    particle.externalMegnField.at(2)=fieldzr[0];
    
    //set the external field of synparticle to zero. Once bending, what should be done???? 
    if (particle.xyz[0]==0&&particle.xyz[2]==0)
    {
        particle.externalMegnField.at(0)    =0;
        particle.externalMegnField.at(1)    =0;
    }

    return;
}























