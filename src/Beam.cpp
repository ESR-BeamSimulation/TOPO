#include <cmath>
#include <iostream>
#include <fstream>
#include "iomanip"

#include "Beam.h"
#include "MyFunc.h"
#include "Global.h"



void Beam::SetField()
{
  for(unsigned int i=0; i < particleNumber;++i)
  {
    particle[i].GetTotalField();
  }
}
void Beam::SetMatrix()
{   
    for(int i=0; i<6;++i)
    {
        for(int j=0;j<6;++j)
        {
            sigmaMatrix(i,j)  =   sigmaMatrix(i,j);
        }
    }
}

void Beam::CaculateEmittance()
{
  double  particleNumber= particle.size();
 
//  cout<<particleNumber<<endl;
    
    
    v2d  xyz   = v2d(particleNumber,v1d(6)); //to store the particle[i].xyz[j] 
    

    for (int i=0;i<particleNumber;++i)
    {
     
        xyz[i][1] =  particle[i].xyz[1] ;                              
        xyz[i][3] =  particle[i].xyz[3] ;                              
        xyz[i][5] =  particle[i].xyz[5] ;                           
        xyz[i][0] =  particle[i].xyz[0] ;
        xyz[i][2] =  particle[i].xyz[2] ;
        xyz[i][4] =  particle[i].xyz[4] ;     
    }
 
    particleLiveNumber  =   0;
    particleLiveInMesh  =   0;
    
    
    for (int i=0;i<6;++i)
    {
        for (int j=0;j<6;j++)
        {
            total2[i][j]=0;
//            sigmaMatrix(6*i +j)=0;
            sigmaMatrix(i,j)=0;
        }
        total[i]=0;
        sigma[i]=0;
    }
    
    for(int i=0;i<3;++i)
    {
        alpha[i]    =   0;
        beta[i]     =   0;
        emittance[i]=   0;    
    }


//  only particle in the main bunch are taken into account for emittance caluclation
    for(int i=0;i<particleNumber;i++)
    {
        if(particle[i].transLossFlag)
        {
            continue;
        }
        
        ++particleLiveNumber;
        
        if (particle[i].longiLossFlag)
        {
            continue;
        }        
        ++particleLiveInMesh;
        
        for(int j=0;j<6;++j)
        {
            total[j] = total[j]   + xyz[i][j];         
        }
    }   
  
  
    for(int j=0;j<6;++j)
    {
        average[j]  =   total[j]   /particleLiveInMesh;
    }  

    
    
    
    for (int i=0;i<particleNumber;++i)
    {
        xyz[i][1] =  (xyz[i][1] - average[1]) / xyz[i][5];                    // px/pz ->rad        
        xyz[i][3] =  (xyz[i][3] - average[3]) / xyz[i][5];                    // py/pz ->rad  
        xyz[i][5] =  (xyz[i][5] - average[5]) / xyz[i][5];                    // pz/pz ->rad  
        xyz[i][0] =  (xyz[i][0] - average[0]) ;
        xyz[i][2] =  (xyz[i][2] - average[2]) ;
        xyz[i][4] =  (xyz[i][4] - average[4]) ;
    }





    for(int i=0;i<particleNumber;i++)
    {
        if(particle[i].transLossFlag||particle[i].longiLossFlag)
        {
            continue;
        }
        
        for(int j=0;j<6;++j)
        {
            for(int k=0;k<6;++k)
            {
                total2[j][k] = total2[j][k] + xyz[i][j] *  xyz[i][k];
            }
        }
    }


//    for(int j=0;j<6;++j)
//    {
//        for(int k=0; k<6;++k)
//        {
//             if(k==j)
//             {cout<<total2[k][k]<<"\t"; }           
//        }
//    }
//    cout<<endl;



    for(int j=0;j<6;++j)
    {
        for(int k=0; k<6;++k)
        {
            total2[j][k] = total2[j][k] /particleLiveInMesh;
        }
    }

 
 
    for(int i=0;i<3;++i)
    {
        emittance[i] = sqrt(total2[2*i][2*i] * total2[2*i+1][2*i+1] - total2[2*i][2*i+1]*total2[2*i+1][2*i]);
    } 
    
    for (int i=0; i<3;i++)
    {
        beta [i]    =        total2[2*i  ][2*i  ] /emittance[i];
        alpha[i]    =   -1 * total2[2*i  ][2*i+1] /emittance[i];
    }
    
    for(int i=0;i<3;++i)
    {
        emittance[i]    =   emittance[i]    * particle[0].xyz[5];
    }

   

    for (int i=0;i<6;++i)
    {
        for (int j=0;j<6;++j)
        {
//            sigmaMatrix(6*i +j) = total2[i][j];
            sigmaMatrix(i,j)=total2[i][j];
        }    
    }

    for (int i=0;i<6;i++)
    {
        sigma[i]    =   sqrt(total2[i][i]);
//        cout<< sigma[i]<<"\t"<<total2[i][i]<<"\t";
    }
//        cout<<endl;

  
//  eGammaTotal		=0;
//  eBetaGammaTotal	=0;
//  phiTotal		=0;
//  eGammaSigma		=0;
//  eBetaGammaSigma	=0;
//  phiSigma		=0;
//  PhiGamma		=0;
//  PhiBetaGamma		=0;

//  for(int i=0;i<particleNumber;i++)
//  {
//    if(particle[i].lossFlag)
//    {
//      continue;
//    }
//    eBeta [i]    = particle[i].xyz[5]/C_light;
//    eGamma[i]	 = 1.0/sqrt(  1-pow(eBeta[i],2)  );
//    eBetaGamma[i]= eBeta[i]*eGamma[i];
//  }

//  for(int i=0;i<particleNumber;i++)
//  {
//    if(particle[i].lossFlag)
//    {
//      continue;
//    }
//    for(int j=0;j<6;++j)
//    {
//      total[j]+=particle[i].xyz[j];
//    }
//    eGammaTotal		+=eGamma[i];
//    eBetaGammaTotal	+=eBetaGamma[i];
//    ++particleLiveNumber;
//  }

//  eBetaTotal	= total[5]/C_light;

//  if(particleLiveNumber==0) 
//  {
//    cout<<"No particle LIVE, send in caculate_emittance"<<endl;
//    return;
//  }
//  for(int i=0;i<6;++i)
//  {
//    average[i]          =       total[i]/particleLiveNumber;
//  }
//  eBetaAverage          =       eBetaTotal	/particleLiveNumber;
//  eGammaAverage         =       eGammaTotal	/particleLiveNumber;
//  eBetaGammaAverage     =       eBetaGammaTotal	/particleLiveNumber;

//  dwwmax=0;
//  dppmax=0;
//  dbbmax=0;

//  for(int i=0;i<particleNumber;i++)
//  {
//    if(particle[i].lossFlag)
//    {
//      continue;
//    }
//    for(int j=0;j<6;j=j+2)
//    {
//      sigma[j]		+= pow( particle[i].xyz[j]-average[j]    	                ,2);
//      sigma[j+1]	+= pow((particle[i].xyz[j+1]-average[j+1])/particle[i].xyz[5]    ,2);
//    }
//    for(int j=0;j<3;++j)
//    {
//      sigmaxdx[j]	+= (particle[i].xyz[2*j]-average[2*j])	*	(particle[i].xyz[2*j+1]-average[2*j+1])	/	particle[i].xyz[5];
//    }
//    eGammaSigma		+= pow( eGamma[i]-eGammaAverage	,2);
//    PhiGamma		+= (eGamma[i]-eGammaAverage)	*	(particle[i].xyz[4]  -average[4]);
//    eBetaGammaSigma	+= pow( eBetaGamma[i]-eBetaGammaAverage	,2);
//    PhiBetaGamma	+= (eBetaGamma[i]-eBetaGammaAverage)*	(particle[i].xyz[4]  -average[4]);	
//    if(abs(eGamma[i]-eGammaAverage)/(eGammaAverage-1)>dwwmax) 
//    {
//      dwwmax=abs(eGamma[i]-eGammaAverage)/(eGammaAverage-1);
//    }
//    if(abs(eBetaGamma[i]-eBetaGammaAverage)/(eBetaGammaAverage)>dppmax) 
//    {
//      dppmax=abs(eBetaGamma[i]-eBetaGammaAverage)/(eBetaGammaAverage);
//    }
//    if(abs(eBeta[i]-eBetaAverage)/(eBetaAverage)>dbbmax) 
//    {
//      dbbmax=abs(eBeta[i]-eBetaAverage)/(eBetaAverage);
//    }
//  }
//  sigmaz	=	eGammaSigma	/ pow(eGammaAverage-1	,2)/particleLiveNumber;
//  eGammaSigma	=	eGammaSigma	/ pow(eGammaAverage	,2);
//  eBetaGammaSigma=	eBetaGammaSigma	/ pow(eBetaGammaAverage	,2);
//  PhiGamma	=	PhiGamma   	/    (eGammaAverage	)	*   (1/average[5]*frq*360  );
//  PhiBetaGamma	=	PhiBetaGamma	/    (eBetaGammaAverage)	*   (1/average[5]*frq*360  );
//  phiSigma	=	sigma[4]   					*pow(1/average[5]*frq*360,2);

//  for(int j=0;j<6;++j)
//  {
//    sigma[j]	/=	particleLiveNumber;
//  }
//  for(int j=0;j<3;++j)
//  {
//    sigmaxdx[j]	/=	particleLiveNumber;
//  }
//  phiSigma	 /=particleLiveNumber;
//  eGammaSigma	 /=particleLiveNumber;
//  PhiGamma	 /=particleLiveNumber;
//  eBetaGammaSigma/=particleLiveNumber;
//  PhiBetaGamma	 /=particleLiveNumber;

//  for(int i=0;i<3;++i)
//  {
//    emit[i]=sqrt(sigma[2*i]*sigma[2*i+1]-sigmaxdx[i]*sigmaxdx[i]);
//  }
//  //x
//  x_sigma	=sigma[0];
//  dx_sigma	=sigma[1];
//  xdx_sigma	=sigmaxdx[0];
//  emittancex	=emit[0];

//  //y
//  y_sigma	=sigma[2];
//  dy_sigma	=sigma[3];
//  ydy_sigma	=sigmaxdx[1];
//  emittancey	=emit[1];

//  //Z
//  z_sigma	=phiSigma;
//  dz_sigma	=eBetaGammaSigma;
//  //dz_sigma	=eGammaSigma;
//  //sigmaz	=eGammaSigma;
//  sigmadz	=sigma[5];
//  zdz_sigma	=PhiBetaGamma;
//  //zdz_sigma	=PhiGamma;
//  emittancez	=sqrt(z_sigma*dz_sigma-zdz_sigma*zdz_sigma);

//  //dz_sigma	=eGammaSigma;
//  //zdz_sigma	=PhiGamma;
//  //emittancez	=sqrt(z_sigma*dz_sigma-zdz_sigma*zdz_sigma);

//  TBetax	=x_sigma/emittancex;
//  TBetay	=y_sigma/emittancey;
//  TBetaz	=z_sigma/emittancez;
//  TGammax	=dx_sigma/emittancex;
//  TGammay	=dy_sigma/emittancey;
//  TGammaz	=dz_sigma/emittancez;
//  if(TBetax*TGammax>1)
//  {
//    TAlphax=(xdx_sigma<0?1:-1)*sqrt(TBetax*TGammax-1);
//  }
//  else
//  {
//    TAlphax=0;
//  }
//  if(TBetay*TGammay>1)
//  {
//    TAlphay=(ydy_sigma<0?1:-1)*sqrt(TBetay*TGammay-1);
//  }
//  else
//  {
//    TAlphay=0;
//  }
//  if(TBetaz*TGammaz>1)
//  {
//    TAlphaz=(zdz_sigma<0?1:-1)*sqrt(TBetaz*TGammaz-1);
//  }
//  else
//  {
//    TAlphaz=0;
//  }










//  ParaTransfer.push_back({
//      average[4] *100,                               //beam position Z       (cm)
//      rgamma * BaseMassInMeV,                   //Kinetic               (MeV)
//      sigma[0] *100,                                 //rms size in X         (cm)
//      sigma[2] *100,                                 //rms size in Y         (cm)
//      sigma[4] *100,                                 //rms size in Z         (cm)
//      sigma[1] *1000,                                 //rms size in Px        (mrad)
//      sigma[3] *1000,                                 //rms size in Py        (mrad)
//      sigma[5] *1000,                                 //rms size in Pz        (mrad)
//      emittance[0]*10e5,                             //emit. nor. x          (cm*mrad)
//      emittance[1]*10e5,                             //emit. nor. y          (cm*mrad)
//      emittance[2]*10e5,                             //emit. nor. z          (cm*mrad)
//      beta[0] *100,                                   //beta                     (cm/rad) 
//      beta[1] *100,
//      beta[2] *100,
//      alpha[0],
//      alpha[1],
//      alpha[2],
//      particleLiveNumber
//      dppmax,
//      dwwmax,
//      dbbmax
//      });
}

void Beam::CaculateSigma()
{
        
    
//  vector<double> total(6);
//  int particleLiveNumber=0;
//  particleNumber=particle.size();
//  for(int i=0;i < particleNumber;i++)
//  {
//    if(particle[i].lossFlag)
//    {
//      continue;
//    }
//    for(int j=0;j<6; j=j+2)
//    {
//      total[j]+=particle[i].xyz[j];
//    }
//    ++particleLiveNumber;
//  }

//  for(int i = 0 ; i < 6 ; i=i+2)
//  {
//    average1[i] = total[i]/particleLiveNumber;
//  }

//  sigma=vector<double>(6,0);
//  for(int i=0;i<particleNumber;i++)
//  {
//    if(particle[i].lossFlag)
//    {
//      continue;
//    }
//    for(int j = 0 ; j<6 ;j=j+2)
//    {
//      sigma[j]   += pow( particle[i].xyz[j],2 );  - average[j],2);
//    }
//  }

//  for(int i = 0 ; i < 6 ; i=i+2)
//  {
//    sigma[i] /= particleLiveNumber;
//    sigma[i]  =  sqrt(sigma[i]);
//  }
}

void Beam::Print(ostream &par)
{  
  
  double rgamma = sqrt( pow( average[5], 2)+1);
  double rbeta  = sqrt( pow(rgamma,2)-1  ) / rgamma;
  double temp, temp1, temp2;
  
  
  temp  = 360 / (C_light/frq) * BaseMassInMeV * 10e3 * average[5]* emittance[2] ; // z emittance keV *ns
  temp1 = 360 / rbeta / (C_light/frq) *  sigma[4] ;                              // z rms size in deg
  temp2 =  pow(rbeta,2) * rgamma * BaseMassInMeV *10e3 *  sigma[5];              // pz rms size in keV

  





  par<<setw(15)<<average[4] *100                               //beam position Z       (cm)
     <<setw(15)<<(rgamma-1) * BaseMassInMeV                        //Kinetic               (MeV)
     <<setw(15)<<sigma[0] *100                                 //rms size in X         (cm)
     <<setw(15)<<sigma[2] *100                                 //rms size in Y         (cm)
     <<setw(15)<<sigma[4] *100                                 //rms size in Z         (cm)
     <<setw(15)<<temp1                                         //rms size in z         (deg)
     <<setw(15)<<sigma[1] *1000                                //rms size in Px        (mrad)
     <<setw(15)<<sigma[3] *1000                                //rms size in Py        (mrad)
     <<setw(15)<<sigma[5] *1000                                //rms size in Pz        (mrad)
     <<setw(15)<<emittance[0]*10e5                             //emit. nor. x          (cm*mrad)
     <<setw(15)<<emittance[1]*10e5                             //emit. nor. y          (cm*mrad)
     <<setw(15)<<emittance[2]*10e5                             //emit. nor. z          (cm*mrad)
     <<setw(15)<<temp                                          //emit. nor.            (Kev*ns) 
     <<setw(15)<<beta[0] *100                                  //beta                     (cm/rad) 
     <<setw(15)<<beta[1] *100                                  //beta                     (cm/rad)
     <<setw(15)<<beta[2] *100                                  //beta                     (cm/rad)
     <<setw(15)<<alpha[0]
     <<setw(15)<<alpha[1]
     <<setw(15)<<alpha[2]
     <<setw(15)<<particleLiveNumber
     <<setw(15)<<particleLiveInMesh<<endl;
  
  

  
  
//  PrintV1(par,ParaTransfer.back());
}


