#include "Field.h"
#include "MyFunc.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <Global.h>

using namespace std;

Field::~Field()
{
}

Field::Field()
{
}

void Field::SetParameter(double posit, const vector<string> & param)
  //(const string &path,double position,double length,double phase,double bfactor,double efactor,int type)
{
  if(param.size()<10)
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
  frq           =       stod(param[5]);
  synPhase      =       stod(param[6])/180*M_PI;
  ke            =       stod(param[7]);
  kb            =       stod(param[8]);
  filename      =       param[9];

  ifstream fin;
  string                str;
  vector<string>        strVec;
  if(type==0)                           //(x,y,z, ex,ey,ez,bx,by,bz) in one datafile 
  {
    fin.open(filename);
    if(!fin.is_open())
    {
      cout<<"field type0 file error : "<<filename<<endl;
      exit(0);
    }
    vector<double>      xyzTemp(3);
    vector<double>      eFieldTemp(3);
    vector<double>      mFieldTemp(3);
    for(int i = 0;i<2;++i)
    {
      //to skip several lines
      getline(fin,str);
    }
    for(int i = 0;fin.peek()!=EOF;++i)
    {
      getline(fin,str);
      //cout<<str<<endl;
      StringSplit(str,strVec);
      xyzTemp[0]        =       stod(   strVec[0])      /1000;          //mm to m
      xyzTemp[1]        =       stod(   strVec[1])      /1000;
      xyzTemp[2]        =       stod(   strVec[2])      /1000   +position;
      gridLocation.push_back(xyzTemp);

      eFieldTemp[0]     =       stod(   strVec[3])*ke   /1000   ;
      eFieldTemp[1]     =       stod(   strVec[4])*ke   /1000;
      eFieldTemp[2]     =       stod(   strVec[5])*ke   /1000;
      eField.push_back(eFieldTemp);

      mFieldTemp[0]=    stod(   strVec[6])*kb   /1000;
      mFieldTemp[1]=    stod(   strVec[7])*kb   /1000;
      mFieldTemp[2]=    stod(   strVec[8])*kb   /1000;
      mField.push_back(mFieldTemp);
    }
    fin.close();
  }
  else if(type == 1)            ////(ex,ey,ez,bx,by,bz) in six datafiles seperately vector(z,(y,(x)))  
  {
    vector<string> filedex={".edx",".edy",".edz",".bdx",".bdy",".bdz"};
    for(int dire=0;dire<6;dire++)
    {
      string filepath=filename+filedex[dire];
      fin.open(filepath);
      if(!fin.is_open())
      {
        cout<<"field type1 file error : "<<filepath<<endl;
        exit(0);
      }
      if(dire==0)
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
      int index=0;

      for(int k=0;k<=nz;++k)
      {
        for(int j=0;j<=ny;++j)
        {
          for(int i=0;i<=nx;++i)
          {
            fin>>temp;
            index=k*(ny+1)*(nx+1) + j*(nx+1) + i;
            field[index][dire]=temp;
          }
        }
      }
      fin.close();
      gridNum=field.size();
    }


    for(int i=0;i<gridNum;++i)
    {
      for(int j=0;j<3;++j)
      {
        field[i][j]=field[i][j]*1.0e6*ke;
      }
      for(int j=3;j<6;++j)
      {
        field[i][j]=field[i][j]*kb;
      }
    }

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
          for(int dire=0;dire<6;dire++)
          {
            fieldout<<field[index][dire]<<" \t";
          }
          fieldout<<endl;
        }
      }
    }
  }
  
  else if(type==2)                                  //solenoid filed  
  {
    
    vector<string> filedex={".bsz",".bsr"};
    for(int dire=0;dire<filedex.size();dire++)
    {
      string filepath=filename+filedex[dire];
      
      fin.open(filepath);
      
      if(!fin.is_open())
      {
        cout<<"field type2 file error : "<<filepath<<endl;
        exit(0);
      }
      if(dire==0)
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
          field[index][dire]=temp*kb;
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
        for(int dire=0;dire<2;dire++)
        {
          fieldout<<field[index][dire]<<"  ";
        }
        fieldout<<endl;
      }
    }
//    int printFlag=0;
//    if(printFlag==1)
//    {
//      ofstream fieldout("fieldout.dat");
//      for(int k=0;k<=nz;++k)
//      {
//        for(int j=0;j<=ny;++j)
//        {
//          for(int i=0;i<=nx;++i)
//          {
//            fieldout<<i*spacex+xmin<<"  ";
//            fieldout<<j*spacey+ymin<<"  ";
//            fieldout<<k*spacez<<"  ";
//            double index=k*(ny+1)*(nx+1)+j*(nx+1)+i;
//            for(int dire=0;dire<6;dire++)
//              fieldout<<field[index][dire]<<"  ";
//            fieldout<<endl;
//          }
//        }
//      }
//    }
  }
  else
  {
    cout<<"Field type error!"<<endl;
  }
}

void Field::GetField(double time,Particle &particle,double synPhaseoffset)
{

    synPhaseoffset = synPhase + synPhaseoffset/180*M_PI;
    
   
//   cout<<synPhaseoffset<<"the synPhaseoffset is "<<endl;
//   getchar();
  
     
  if(type==1)
  {
    if(   particle.xyz[0]<=xmin     || particle.xyz[0]>=xmax
        ||particle.xyz[2]<=ymin     || particle.xyz[2]>=ymax
        ||particle.xyz[4]<=position || particle.xyz[4]>=locationEnd)
    {
      particle.externalElecField=zero3d;
      particle.externalMegnField=zero3d;
      return;
    }
    else if(particle.xyz[4]>locationMax)                    // lies in the date file give by CST, whihc is not exactly as long as field map
    {
      particle.externalElecField=zero3d;
      particle.externalMegnField=zero3d;
      return;
    }
   

        
    double xn,yn,zn;
    int xGrid,yGrid,zGrid,indexGrid;
    v1d v(8),fieldtemp(6);
    v2d fieldOnGrid(6,v1d(8));

//    particle.xyz[0] =   0;
//    particle.xyz[2] =   0;
    
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
          xGrid=((i%2< 1)? ceil(xn)  : ceil(xn)-1);
          yGrid=((i%4< 2)? ceil(yn)  : ceil(yn)-1);
          zGrid=((i  < 4)? ceil(zn)  : ceil(zn)-1);
        
          indexGrid=zGrid*(ny+1)*(nx+1)+yGrid*(nx+1)+xGrid; 
        

//          // once particle on mesh, need to consider that part.
//          if(ceil(xn)==xn)
//          {
//             xGrid +=1;
//          }
//          if(ceil(yn)==yn)
//          {
//             yGrid +=1;
//          }
//          if(ceil(zn)==zn)
//          {
//             zGrid +=1;
//          }


      
          
          v[i]=abs(xn-xGrid) * abs(yn-yGrid) * abs(zn-zGrid);
        
            
    //      cout<<xn<<"\t"<<yn<<"\t"<<zn<<"\t"<<xGrid<<"\t"<<yGrid<<"\t"<<zGrid<<"\t"<<v[i]<<"\t\t"<< field[indexGrid][0]
    //                 <<"\t\t"<<field[indexGrid][1]<<"\t\t"<<field[indexGrid][2]<<endl;
    //getchar();
         

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
    
    //set the external field of synparticle to zero. Once bending, what should be done???? 
    if(particle.xyz[0]==0&&particle.xyz[2]==0)
    {
        for (int i=0;i<2;++i)
        {
            particle.externalMegnField.at(i)    =0;
            particle.externalElecField.at(i)    =0;
        }
    }

    
//    cout<<particle.externalElecField[2]<<endl;
  }
  else if(type==2)
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
        fieldzr[0]  =   abs(zn-floor(zn)-1) * field[indexGrid][0] +  abs(zn-floor(zn)) *  field[(zGrid+1)*(nr+1)+rGrid][0];
    }
    else
    {
        for(int i=0;i<4;++i)
        {
          rGrid=((i%2< 1)?ceil(rn):floor(rn));
          zGrid=((i  < 2)?ceil(zn):floor(zn));
          v[i]=abs(rn-rGrid) * abs(zn-zGrid);

          rGrid=((i%2==1)?ceil(rn):floor(rn));
          zGrid=((i  >=2)?ceil(zn):floor(zn));
          indexGrid=zGrid*(nr+1)+rGrid;

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
    if(particle.xyz[0]==0&&particle.xyz[2]==0)
    {
        for (int i=0;i<2;++i)
        {
            particle.externalMegnField.at(i)=0;
        }
    }
    
  }
}

void Field::Synchronize(double time,Beam &beam)
{
  beam.CaculateSigma();
  double gamma, beta;
  gamma =   sqrt(1+pow(beam.average[5],2)+ pow(beam.average[3],2)+pow(beam.average[1],2));
  beta  =   sqrt(pow(gamma,2)-1)/gamma;

  double distance       = (position+length/2)-beam.average[4];
  double deltaTime      = distance/beta/C_light;
  offsetPhase           =-2*M_PI*frq*deltaTime;
  synchronized=1;

}

void Field::ScanOffset(double time, Beam &beam)
{

}


void Field::GetMatrix(Particle &)
{
}
