#include "FieldEM.h"
#include "MyFunc.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <Global.h>

using namespace std;

FieldEM::~FieldEM()
{
}

FieldEM::FieldEM()
{
}

void FieldEM::SetParameter(double posit, const vector<string> & param)
{
    if(param.size()<9)
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
    type          =       stoi(param[3]);
    frq           =       stod(param[4]);
    synPhase      =       stod(param[5])/180*M_PI;
    ke            =       stod(param[6]);
    kb            =       stod(param[7]);
    filename      =       param[8];
    
    ifstream fin;
    string                str;
    vector<string>        strVec;
    
    vector<string> filedex={".edx",".edy",".edz",".bdx",".bdy",".bdz"};
    
    for(int dim=0;dim<6;dim++)
    {
        string filepath=filename+filedex[dim];
        fin.open(filepath);
        
        if(!fin.is_open())
        {
            cout<<"field type1 file error : "<<filepath<<endl;
            exit(0);
        }
        
       if(dim==0)
       {
            fin>>nz>>zmax;
            fin>>nx>>xmin>>xmax;
            fin>>ny>>ymin>>ymax;
            field.resize((nx+1)*(ny+1)*(nz+1),vector<double>(6,0));
            locationMax=position+zmax;
            spacex=(xmax-xmin)/nx;
            spacey=(ymax-ymin)/ny;
            spacez=(zmax     )/nz;
            getline(fin,str);
            getline(fin,str);
      }
      else
      {
          for(int i = 0;i<4;++i)
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
      int    index=0;
      
      for(int k=0;k<=nz;++k)
      {
        for(int j=0;j<=ny;++j)
        {
          for(int i=0;i<=nx;++i)
          {
            fin>>temp;
            index=k*(ny+1)*(nx+1) + j*(nx+1) + i;
            field[index][dim]=temp;
          }
        }
      }
      
      fin.close();
      gridNum=field.size();
    
    }
    
    
    // obtain the data for latter use
    for(int i=0;i<gridNum;++i)
    {
        for(int j=0;j<3;++j)
        {
            field[i][j]=field[i][j]*1.0e6*ke;       //MV -> V
        }
          
        for(int j=3;j<6;++j)
        {
            field[i][j]=field[i][j]*kb;             // T ->T
        }
    }
    
    
    
    // print the date to data file
    string fout=filename+"_fieldout"+".dat";
    ofstream fieldout(fout); 
    for(int i=0;i<=nx;++i)
    {
        for(int j=0;j<=ny;++j)
        {
            for(int k=0;k<=nz;++k)
            {
                fieldout<<i*spacex+xmin<<"  ";
                fieldout<<j*spacey+ymin<<"  ";
                fieldout<<k*spacez<<"  ";
                
                double index=k*(ny+1)*(nx+1)+j*(nx+1)+i;
                for(int dim=0;dim<6;dim++)
                {
                    fieldout<<field[index][dim]<<" \t";
                }
                
                fieldout<<endl;
            }
        }
    }

     return;
}

void FieldEM::GetField(double time,Particle &particle,double synPhaseoffset)
{
    synPhaseoffset = synPhase + synPhaseoffset/180*M_PI;
    
    if(   particle.xyz[0]<xmin     || particle.xyz[0]>xmax
        ||particle.xyz[2]<ymin     || particle.xyz[2]>ymax
        ||particle.xyz[4]<position || particle.xyz[4]>locationEnd)
    {
        particle.externalElecField=zero3d;
        particle.externalMegnField=zero3d;
        return;
    }
    else if(particle.xyz[4]>locationMax)    // depend on date given by CST!!!! 
    {                                
      particle.externalElecField=zero3d;
      particle.externalMegnField=zero3d;
      return;
    }
    
    double xn,  yn,   zn;
    int   xGrid,yGrid,zGrid,indexGrid;
    v1d   v(8,0), fieldtemp(6,0);
    v2d   fieldOnGrid(6,v1d(8));
    double vx,vy,vz;
    
    xn=(particle.xyz[0]-xmin)    /spacex;
    yn=(particle.xyz[2]-ymin)    /spacey;
    zn=(particle.xyz[4]-position)/spacez;
    
    
    if( particle.xyz[0]==0 &&  particle.xyz[2]==0 ) 
    {
        // only the Ez filed is calculated in fieldEM
        xGrid   =   floor(xn);
        yGrid   =   floor(yn);
        zGrid   =   floor(zn);
        indexGrid=zGrid*(ny+1)*(nx+1)+yGrid*(nx+1)+xGrid;
        
        fieldtemp[2]  = abs(zn-floor(zn)-1) * field[indexGrid][2] + abs(zn-floor(zn)) *  field[(zGrid+1)*(ny+1)*(nx+1)+yGrid*(nx+1)+xGrid][2];
        fieldtemp[5]  = abs(zn-floor(zn)-1) * field[indexGrid][5] + abs(zn-floor(zn)) *  field[(zGrid+1)*(ny+1)*(nx+1)+yGrid*(nx+1)+xGrid][5];
        fieldtemp[0]  = 0;
        fieldtemp[1]  = 0;
        fieldtemp[3]  = 0;
        fieldtemp[4]  = 0;
    }
    else
    {
        for(int i=0;i<8;++i)
        {
            xGrid   =   ((i%2< 1)? ceil(xn)  : ceil(xn)-1);
            yGrid   =   ((i%4< 2)? ceil(yn)  : ceil(yn)-1);
            zGrid   =   ((i  < 4)? ceil(zn)  : ceil(zn)-1);  
            
            indexGrid=zGrid*(ny+1)*(nx+1)+yGrid*(nx+1)+xGrid;
            
              //do no consider the cases that  once particles is loacted on mesh 
//            if(ceil(xn)==xn)
//            {
//                vx  =   1;
//            }
//            else
//            {
//                vx  =   abs(xn-xGrid);
//            }
//            if(ceil(yn)==yn)
//            {
//                vy   =   1;
//            }
//            else
//            {
//                vy  =   abs(yn-yGrid);
//            }
//            if(ceil(zn)==zn)
//            {
//               vz   =   1;
//            }
//            else
//            {
//                vz  =   abs(zn-zGrid);
//            }
            
            v[i]=abs(xn-xGrid) * abs(yn-yGrid) * abs(zn-zGrid);
        
            for(int dim=0;dim<6;++dim)
            {
                fieldOnGrid[dim][i] =   v[i]*field[indexGrid][dim];
                fieldtemp[dim]     +=   fieldOnGrid[dim][i];
            } 
        }
    }
    
    double phaseCos=cos(2*M_PI*frq*time + synPhaseoffset);  //    -offsetPhase);
   
    for(int dim=0;dim<3;++dim)
    {
      particle.externalElecField[dim]=fieldtemp[dim  ] * phaseCos;
      particle.externalMegnField[dim]=fieldtemp[dim+3] * phaseCos;
    }
    
    if(particle.xyz[0]==0&&particle.xyz[2]==0)
    {
        particle.externalElecField[0]   = 0;
        particle.externalElecField[1]   = 0;
        particle.externalMegnField[0]   = 0;
        particle.externalElecField[1]   = 0;
    }
    
    
    return;  
}


    
    
 
