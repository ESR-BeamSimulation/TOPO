#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>

#include "Beam.h"
#include "Lattice.h"
#include "Initial.h"
#include "Distribution.h"
#include "MyFunc.h"
#include "AllElement.h"
#include "Global.h"
#include "LeapFrog.h"
#include <sys/time.h>   
#include <strstream>

using namespace std;

int Initial::Read()
{
  ifstream in("input.txt");	//read from file and
  //generate the partilce for give type.
  if(!in.is_open())
  {
    cout<<"File \"input.txt\" not opened!"<<endl;
    exit(0);
  }


  vector<string> strVec;
  string         str;

  vector<double> param;
  
  cout<<"       Chao Li      /IMP/IHEP~~~~lichao@ihep.ac.cn"<<endl;
  cout<<"       Zhicong Liu  /IHEP/LBNL"<<endl;
  cout<<"       Yaliang Zhao /IHEP"<<endl;
  cout<<"       Qing Qin     /IHEP"<<endl;
  cout<<"*************************************************************************"<<endl;

  cout<<"        check the initial parameter of input                             "<<endl;
  
  while(in.peek()!=EOF)         //firstly, to know 2D or 3D,type
  {
    getline(in,str);            //read one line fron the file
    StringSplit(str,strVec);    //Split the line to several string and store in strVec
    if(strVec.empty())
    {
      continue;
    }
    for(int i=0;i<strVec.size();++i)
    {
      transform   (strVec[0].begin(), strVec[0].end(), strVec[0].begin(), ::tolower);     //change to lower bucket
    }

    input.push_back(str);

    //the following section is for input parameter according to different key word.
    if(strVec[0]=="method")
    {
      method = strVec[1];
      cout<<"the method                                         :"<< method<<endl;
    }
    else if(strVec[0]=="scfieldmethod")
    {
      internalFieldMethod = strVec[1];
      cout<<"the space charge field is obtained with            :"<<internalFieldMethod<<endl;
    }
    else if(strVec[0]=="dimension")
    {
      if(strVec[1]=="2")
      {
        dimension=2;
        cout<<"the space charge is limit to                       :"<<dimension<<"D continous beam"<<endl;
      }
      else if(strVec[1]=="3")
      {
        dimension=3;
        cout<<"the space charge is limit to                       :"<<dimension<<"D bunched beam"<<endl;
      }
      else
      {
        cout<<"Dimension Error!"<<endl;
        cout<<"The Format is"<<endl;
        cout<<"  Dimension  2D/3D"<<endl;
        exit(0);
      }
    }
    else if(strVec[0]=="distribution")
    {
      if(strVec.size()!=3)
      {
        cout<<"we have two indexs for distribution generation";
        cout<<"separately in transvers and longitudinal"<<endl;
        exit(0);
      }
      for(int i=0;i<2;++i)
      {
        distributionType.at(i)=strVec[i+1];
        string distype=strVec[i+1];
        if(distype!="KV" && distype!="WB"&& distype!="PB"&& distype!="GS")
        {
          cout<<"Distribution Type Error!"<<endl;
          cout<<"The Format is"<<endl;
          cout<<"  Distribution  KV/WB/PB/GS  KV/WB/PB/GS"<<endl;
          exit(0);
        }
      }
        cout<<"the transverse distrbution                         :"<<distributionType[0]<<" type"<<endl;
        cout<<"the longitudinal distrbution                       :"<<distributionType[1]<<" type"<<endl;
    }
    else if(strVec[0]=="particle")
    {
      numofCharge = stod(strVec[1]);
      numofMass   = stod(strVec[2]);
      cout<<"the charge of each particle                        :"<<numofCharge<<endl;
      cout<<"the mass number of each paritlce                   :"<<numofMass<<endl;
    }
    else if(strVec[0]=="twissx")
    {
      if(strVec.size()!=4)
      {
        cout<<"twissx have 3 parameters"<<endl;
        cout<<"emittance, beta, alpha"<<endl;
        exit(0);
      }
      emitx     = stod(strVec[1]);
      betax     = stod(strVec[2]);
      alphax    = stod(strVec[3]);
      for(int i=0;i<3;++i)
      {
        twiss[i]=stod(strVec[i+1]);
      }
      cout<<"the normalized emittancex                          :"<<emitx<<" m *mrad"<<endl;
      cout<<"the twiss alphax                                   :"<<alphax<<endl;
      cout<<"the twiss betax                                    :"<<betax<<" m/rad"<<endl;
    }
    else if(strVec[0]=="twissy")
    {
      if(strVec.size()!=4)
      {
        cout<<"twissy have 3 parameters"<<endl;
        cout<<"emittance, beta, alpha"<<endl;
        exit(0);
      }
      emity     = stod(strVec[1]);
      betay     = stod(strVec[2]);
      alphay    = stod(strVec[3]);
      for(int i=0;i<3;++i)
      {
        twiss[i+3]=stod(strVec[i+1]);
      }
      cout<<"the normalized emittancey                          :"<<emity<<" m *mrad"<<endl;
      cout<<"the twiss alphay                                   :"<<alphay<<endl;
      cout<<"the twiss betay                                    :"<<betay<<" m/rad"<<endl;
    }
    else if(strVec[0]=="twissz")
    {
      if(strVec.size()!=4)
      {
        cout<<"twissz have 3 parameters"<<endl;
        cout<<"emittance, betax, alpha"<<endl;        
        exit(0);
      }
      emitz     = stod(strVec[1]);
      betaz     = stod(strVec[2]);
      alphaz    = stod(strVec[3]);
      for(int i=0;i<3;++i)
      {
        twiss[i+6]=stod(strVec[i+1]);
      }
      cout<<"the normalized emittancez                          :"<<emity<<" m *mrad"<<endl;
      cout<<"the twiss alphaz                                   :"<<alphay<<endl;
      cout<<"the twiss betaz                                    :"<<betay<<" m/rad"<<endl;
    }
    else if(strVec[0]=="particlenumber")
    {
      particleNumber    = stoi(strVec[1]);
      cout<<"the number of macroparticles are                   :"<<particleNumber<<endl;
    }
    else if(strVec[0]=="displacepos")
    {
      for(int i=0;i<strVec.size()-1;++i)
      {
        displacePos[i]=stod(strVec[i+1]);
      }
      cout<<"the initial colletive displacemetn of x y z        :"<<displacePos[0]<<" "<<displacePos[1]<<" "<<displacePos[2]<<" m"<<endl;
    }
    else if(strVec[0]=="displacedpos")
    {
      for(int i=0;i<strVec.size()-1;++i)
      {
        displaceDPos[i]=stod(strVec[i+1]);
      }
      cout<<"the initial collective tile displace of x y z      :"<<displaceDPos[0]<<" "<<displaceDPos[1]<<" "<<displaceDPos[2]<<" rad"<<endl;
    
    }
    
    else if(strVec[0]=="frequency")
    {
      frequency = stod(strVec[1]);          //Hz
      cout<<"the baisc frequency                                :"<<frequency<<" Hz"<<endl;
    }
    
    else if(strVec[0]=="kneticenergy")
    {
      KneticE   =   stod(strVec[1]);        //MeV
      dw        =   stod(strVec[2]);        //energy spread
      cout<<"the initial beam average energy and energy spread  :"<<KneticE<<" MeV "<<dw*100<<"% "<<endl;
    }
    else if(strVec[0]=="beamlength")
    {
      beamLength[0] = stoi(strVec[1]);
      beamLength[1] = stod(strVec[2]);
      cout<<"The generated beamlength is                         :"<<beamLength[0] <<"RFperiods"<<endl;
    }
    else if(strVec[0]=="numofgrid")
    {
      for(int i=0;i<strVec.size()-1;++i)
      {
        numofGrid[i]=stoi(strVec[i+1]);
      } 
    }
    else if(strVec[0]=="current")
    {
      current   = stod(strVec[1]);
      cout<<"the beam current                                   :"<<current<<" Amp"<<endl; 
      if(strVec[1]=="0")
      {
        spaceChargeFlag	=false;
      }
      else if(stod(strVec[1])>0)
      {
        spaceChargeFlag	=true;
        cout<<"the mesh numerb for space charge nx, ny, nz        :"<<numofGrid[0]<<" "<<numofGrid[1]<<" "<<numofGrid[2]<<endl;
      }
      else
      {
        cout<<"current cannot be set less than 0"<<endl;
      }
    }
    else if(strVec[0]=="pusher")
    {
      pusher        = strVec[1];                //LeapFrog RK4 etc...
      transform(pusher.begin(), pusher.end(), pusher.begin(), ::tolower);
      cout<<"the particle push subrutine                        :"<<pusher<<endl;
    }
    else if(strVec[0]=="steppercycle")
    {
      stepPerCycle  = stoi(strVec[1]);
       cout<<"The number of time step in each RF period          :"<<stepPerCycle<<endl;
    }
//    else if(strVec[0]=="maxstepnumber")
//    {
//      maxStepNumber = stoi(strVec[1]);
//      cout<<"The maximum time step in each run                   :"<<maxStepNumber<<endl;
//    }
    else if(strVec[0]=="bunchnumber")
    {
      bunchNumber = stoi(strVec[1]);
      cout<<"The total bunch number                               :"<<bunchNumber<<endl;
    }
    else if(strVec[0]=="printstep")
    {
      printStep = stoi(strVec[1]);
      cout<<"the time steps between each print                    :"<<printStep<<endl;
    }
    else if(strVec[0]=="runperiod")
    {
      runPeriodNumber = stoi(strVec[1]);
      cout<<"the running period number is limited to               :"<<runPeriodNumber<<endl;
    }
    else if(strVec[0]=="traj")
    {
      traj = stoi(strVec[1]);
      cout<<"the number of particles for trajectories recording    :"<<runPeriodNumber<<endl;
    }
    else if(strVec[0]=="setsyn")
    {
      setsynflag = stoi(strVec[1]);
      cout<<"only syn particle is used to set the syn phase of cav  :"<<setsynflag<<endl;
    }
    else
    {
      cout<<"error, unknown item in input file"<<endl;
    }
  }
  cout<<"*************************************************************************"<<endl;
    
    getchar();
    in.close();
    
}

int Initial::InitParticle(Beam &beam)
{
  
  Read();

  beam.frq              =frequency;
  

//  beam.eGamma           .resize(particleNumber);
//  beam.eBeta            .resize(particleNumber);
//  beam.eBetaGamma       .resize(particleNumber);
//  beam.dw               .resize(particleNumber);
//  beam.phi              .resize(particleNumber);
  
  double gamma          =       1.0+KneticE/(BaseMassInMeV);
  double beta           =       sqrt(1-1.0/gamma/gamma);
  double speed          =       beta*C_light;
  double pz0            =       beta*gamma;
  double lambda         =       C_light/beam.frq;
  double rfPeriodLength =       beta *lambda;
  
  double macroNumber;
  

  
  if(setsynflag==1)
  {
    this -> particleNumber=1;
    beam.particleNumber =1;
    beam.particle.resize(beam.particleNumber);
    beam.particle[0].xyz[4] = 0;
    beam.particle[0].xyz[5] = pz0; 
    beam.particle[0].xyz[0] = 0;           
    beam.particle[0].xyz[1] = 0;    
    beam.particle[0].xyz[2] = 0; 
    beam.particle[0].xyz[3] = 0;
    beam.particle[0].indexElem        = 0;
    beam.particle[0].transLossFlag    = 0;
    beam.particle[0].longiLossFlag    = 0;
    beam.particle[0].lossFlag         = 0;
    beam.particle[0].charge           = numofCharge * BaseCharge;             //unit is charge C
    beam.particle[0].mass             = numofMass   * BaseMass;               //unit is kg
    beam.particle[0].qOverMass        = numofCharge / numofMass  * BaseCharge /BaseMass;    //unit is C/KG
    beam.particle[0].outSCMeshFlag    = 0;  
    beam.cossyn                       = 0;
           
    cout<<beam.particle.size()<<"ssss"<<endl;
    getchar();
    return 0;
    
  }

     beam.particleNumber   =    this->particleNumber;
    
 
// cout<<beam.particleNumber<<"what is wrong"<<endl; 
// getchar();
 
  beam.particle.resize(particleNumber);

  if(spaceChargeFlag)
  {
    if(dimension==3)
        {
            macroNumber    =       current  / frequency  / particleNumber / BaseCharge;  //
        }
    else if(dimension==2)
        {
            macroNumber    =       current / frequency  /  particleNumber / rfPeriodLength;
        }
    else
    {
        cout<<"space charge dimension is wrong, it is 2 or 3"<<endl;
    }
  }
  else
  {
    macroNumber    =1;
  }
  
  twiss[0]=twiss[0]/(beta*gamma)/1000;    //emittance is normalized and units change form m*mrad->m*rad
  twiss[3]=twiss[3]/(beta*gamma)/1000;
  twiss[6]=twiss[6]/(beta*gamma)/1000;

 

  vector<vector<double> > xyz(particleNumber,vector<double>(6));
  

  if(     distributionType[0]=="KV")
  {
    KV2DTwiss(twiss,xyz); 
  }
  else if(distributionType[0]=="WB")
  {
    WB2DTwiss(twiss,xyz);
  }
  else if(distributionType[0]=="GS")
  {
    GS2DTwiss(twiss,xyz);
  }
  else if(distributionType[0]=="PB")
  {
    PB2DTwiss(twiss,xyz);
  }
  else
  {
    //do nothing
  }
  for(int i=0;i<particleNumber;++i)
  {
    for(int j=0;j<4;++j)
    {
      beam.particle[i].xyz[j] = xyz[i][j];
    }
    beam.particle[i].indexElem        = 0;
    beam.particle[i].transLossFlag    = 0;
    beam.particle[i].longiLossFlag    = 0;
    beam.particle[i].lossFlag         = 0;
    beam.particle[i].charge           = macroNumber * numofCharge * BaseCharge;             //unit is charge C
    beam.particle[i].mass             = macroNumber * numofMass   * BaseMass;               //unit is kg
    beam.particle[i].qOverMass        = numofCharge / numofMass  * BaseCharge /BaseMass;    //unit is C/KG
    beam.particle[i].outSCMeshFlag    = 0;  
  }
 
  

  if(dimension==3)
  {
    vector<double> twissLongi(9);
    
    cout<<"3D bunched beam is generated with input twiss parameters"<<endl;
   
    for(int i=0;i<3;++i)
    {
      twissLongi[i]=twiss[i+6];
      twissLongi[i+3]=twiss[i+6];
    }
    if( distributionType[1]=="KV")
    {
      KV2DTwiss(twissLongi,xyz);
    }
    else if(distributionType[1]=="WB")
    {
      
      WB2DTwiss(twissLongi,xyz);
    }
    else if(distributionType[1]=="GS")
    {
     
      GS2DTwiss(twissLongi,xyz);
    }
    else if(distributionType[1]=="PB")
    {
      
      PB2DTwiss(twissLongi,xyz);
    }
    else
    {
      //do nothing
    }
    
    for(int i=0;i<particleNumber;++i)
    {
      for(int j=0;j<2;++j)
      {
        beam.particle[i].xyz[j+4] = xyz[i][j];
      }
    }
  }
 


  //  double gammatemp,betatemp;
  //  default_random_engine generator;
  //  uniform_real_distribution<double> distribution2(1-dw,1+dw);
 
  

  double temp, temp1;

  if(beamLength[1]!=0)
  {
    cout<<"the z distribtion is replaced by EK, dw, beamlength is  "<< beamLength[1]<<" RFperiod"<<endl;
   
    for(int i=0;i<particleNumber;i++)
    {
       temp= (double)rand() / RAND_MAX;
       temp1=(double)rand() / RAND_MAX;
       

       
       beam.particle[i].xyz[4] = 0;
       beam.particle[i].xyz[5] = 0;       
       
       beam.particle[i].xyz[4] = beam.particle[i].xyz[4] + (temp-0.5) *beamLength[1] *rfPeriodLength + displacePos [2];;   //  \delta z-> cm    
       beam.particle[i].xyz[5] = gamma/(gamma+1) * dw * 2*(temp1-0.5) ;                                 //  pz/pz0  -> rad 
       
       beam.particle[i].xyz[0] = beam.particle[i].xyz[0] + displacePos [0] + displaceDPos[0]* (beam.particle[i].xyz[4]-beam.particle[0].xyz[4]);  // \deltax  ->cm
       beam.particle[i].xyz[1] = beam.particle[i].xyz[1] ;  //  Px/Pz   ->rad 
       beam.particle[i].xyz[2] = beam.particle[i].xyz[2] + displacePos [1] + displaceDPos[1]* (beam.particle[i].xyz[4]-beam.particle[0].xyz[4]) ;  // \deltay  ->cm
       beam.particle[i].xyz[3] = beam.particle[i].xyz[3] ;  //  Py/Pz   ->rad 
        

    
    }    
  }
  else
  {
    for(int i=1;i<particleNumber;i++)
    {   
        beam.particle[i].xyz[4] = beam.particle[i].xyz[4] + displacePos [2];   //  \delta z-> cm
        beam.particle[i].xyz[5] = beam.particle[i].xyz[5] ;   //  pz/pz0  -> rad
        
        beam.particle[i].xyz[0] = beam.particle[i].xyz[0] + displacePos [0] + displaceDPos[0]* (beam.particle[i].xyz[4]-beam.particle[0].xyz[4]);  // \deltax  ->cm
        beam.particle[i].xyz[1] = beam.particle[i].xyz[1] ;  //  Px/Pz   ->rad 
        beam.particle[i].xyz[2] = beam.particle[i].xyz[2] + displacePos [1] + displaceDPos[1]* (beam.particle[i].xyz[4]-beam.particle[0].xyz[4]) ;  // \deltay  ->cm
        beam.particle[i].xyz[3] = beam.particle[i].xyz[3] ;  //  Py/Pz   ->rad 
    }
  }
   

  beam.CaculateSigma();



//print out the initial particle distribution in *dis.dat
  int particleOutFlag=1;
  if(particleOutFlag)
  {
    ofstream pout("initial particle dis.dat");
    pout<<setiosflags(ios::showpos);
    pout.setf(ios::showpoint);
    pout<<"   posx(m)   "<<"     px/pz(rad)   "<<"      posy(m)   "<<"    py/pz(rad)   "
        <<"      posz(m)   "<<"        pz/pz(rad)   "<<endl;
    
    for(int i = 0;i<beam.particle.size();++i)
    {
      for(int j =0;j<6;++j)
      {
        pout<<setprecision(8)<<setw(12)<<beam.particle[i].xyz[j]<<"       "; 
      }
        pout<<endl;
    }
    pout.close();
  }


// change the coordinate to (x,px,y,py,z,pz) for the coordinate updating
//  ofstream pout1("randamtest.dat");

  for(int i=0;i<particleNumber;++i)
  {
        beam.particle[i].xyz[4] = beam.particle[i].xyz[4];
        beam.particle[i].xyz[5] = beam.particle[i].xyz[5] * pz0 +pz0; 
        beam.particle[i].xyz[0] = beam.particle[i].xyz[0];                              //m
        beam.particle[i].xyz[1] = beam.particle[i].xyz[1] * beam.particle[i].xyz[5];    //beta*gama
        beam.particle[i].xyz[2] = beam.particle[i].xyz[2]; 
        beam.particle[i].xyz[3] = beam.particle[i].xyz[3] * beam.particle[i].xyz[5] ;       
//        for (int j=0;j<6;j++)      
//        {
//            pout1<<beam.particle[i].xyz[j]<<"\t"; 
//        }
//        pout1<<endl;
  }
  
  // set the first particle as the syn-paticle for tracking  
  beam.particle[0].xyz[0]   =   0;
  beam.particle[0].xyz[1]   =   0;
  beam.particle[0].xyz[2]   =   0;
  beam.particle[0].xyz[3]   =   0;
  beam.particle[0].xyz[4]   =   0 + displacePos [2];
  beam.particle[0].xyz[5]   =   pz0;
  
  
  
  
  
  
  
  return 0;
}

int Initial::InitLattice(Lattice &lattice)
{
  ifstream infile("lattice.txt");
  if(!infile.is_open())
  {
    cout<<"Error! Cannot open \"lattice.txt\".\nexit."<<endl;
    exit(0);
  }
  vector<string>        strVec;
  string                str;
  cout<<"initial lattice"<<endl;
  
  int lattcIndex=0;
  
  while(infile.peek()!=EOF)
  {
    getline(infile,str);
    StringSplit(str,strVec);   //element information is store in serVec
    if(strVec.size()==0)
    {
      continue;
    }
    
    if(strVec[0]=="drift")
    {
      ++lattcIndex;
      Element *ele      = new Drift;
      lattice.SetElement(strVec,ele);
      ele   = NULL;
    }
    else if(strVec[0]=="quadm")
    {
      Element *ele      = new Quadrupole;
      lattice.SetElement(strVec,ele);
      ele   = NULL;
      
    }
    else if(strVec[0]=="rfq")
    {    
      Element *ele      = new RFQ;
      lattice.SetElement(strVec,ele);
      ele   = NULL;               
    }
    else if(strVec[0]=="field")
    {
      Element *ele      = new Field;
      lattice.SetElement(strVec,ele);
      ele   = NULL;
    }
    else if(strVec[0]=="fieldem")
    {
      Element *ele      = new FieldEM;
      lattice.SetElement(strVec,ele);
      ele   = NULL;
    }
    else if(strVec[0]=="fieldsb")
    {
      Element *ele      = new FieldSB;
      lattice.SetElement(strVec,ele);
      ele   = NULL;
    }
//    else if(strVec[0]=="fieldse")
//    {
//      Element *ele      = new FieldSE;
//      lattice.SetElement(strVec,ele);
//      ele   = NULL;
//    } 
    
    else if(strVec[0]=="end")
    {
        cout<<"the totoal length is "<< lattice.positionBack.back() <<" m"<<endl;
        getchar();
        return 0;
    }
    

  }

  if( lattice.element.empty() )
  {
    cout<<"Error! \"lattice.txt\" has no effective element.\nexit."<<endl;
    exit(0);
  }
  

  
  for(int i=0;i<lattice.element.size();++i)
  {
     cout<<i<<"\t"<<lattice.position[i]<<"\t"<< lattice.element.at(i)->name<<endl;
  }
  
  
    cout<<"the totoal length is "<< lattice.positionBack.back() <<" m"<<endl;
    getchar();
  
}


int Initial::InitMesh(Mesh &mesh)
{
  mesh.Initial(numofGrid,particleNumber);
}

int Initial::InitMethodMP(MethodMP &MP)
{

  MP.timeNow        =   0;
  MP.timeStep      =    1.0/frequency/stepPerCycle;
  MP.stepMax       =    stepPerCycle * runPeriodNumber;
  MP.printStep     =    printStep;
  MP.spaceChargeFlag=   spaceChargeFlag;
  MP.pusherType    =    pusher;
  MP.traj          =    traj; 
}


int Initial::SetSynPhase(Lattice &lattice)
{

// now, one synparticle is used to set the initial phase for each cavity. The initial phase is saved in the synphaseset.dat.
// the result is not that good, should try to get the initial phase with beam instead of synparticle. 


    Beam          synBeam;
    InitParticle (synBeam);
    MethodMP      synParticleRun;
    InitMethodMP (synParticleRun);
    LeapFrog    leapfrog;

    if(setsynflag==0)
    {
        return 0;
    }
            
           
    double freq;
    double time2,time3;

    int    synPhaseset  = 0;
    int    numoffiled   = 0;
    double desigw[16]={3.2, 3.2, 3.494, 3.825, 4.2, 4.624,5.091,5.469,5.848, 6.106, 6.498, 7.158, 7.851, 8.565, 9.294, 10.036};

    
    
    ofstream test("single_particle_opti.dat");
    ofstream synphaseset("synphaseset.dat");
     
   
    
    for (int i=0;i<synParticleRun.stepMax;++i)
    {
        

        lattice.        GetExternalField    (synParticleRun.timeNow, synBeam);
        synBeam.        SetField            (                               );
        leapfrog.       Update              (synParticleRun.timeStep, synBeam, lattice);
        synParticleRun.timeNow  += synParticleRun.timeStep;
        
        cout<<i<<"\t"<<synBeam.particle[0].xyz[4]<<"\t"<<lattice.element[synBeam.particle[0].indexElem]->name<<"\t"
            <<synBeam.particle[0].xyz[3]<<"\t"<<synBeam.particle[0].xyz[2]<<endl;
            
        
        if(lattice.element[synBeam.particle[0].indexElem]->name=="field")
        {
            freq    =  stod(lattice.element[synBeam.particle[0].indexElem]->parameterString.at(5));
            if(freq!=0)
            {

                double zpos       =   synBeam.particle[0].xyz[4];
                double pz0        =   synBeam.particle[0].xyz[5];
                double win        =   (sqrt(1+pow(pz0,2)) -1) * BaseMassInMeV;
                synBeam.cossyn    = stod(lattice.element[synBeam.particle[0].indexElem]->parameterString.at(6));
                double wou        = 0 ;
                double pz1        = 0 ;
                int    deltastep  = 0 ;
                double deltatime  = 0 ;
                int     count     = 0 ;  
            
                ostringstream s1;
                string s2;
                s1 << "phase_energy_filed" << numoffiled << ".dat";
                s2 = s1.str();
                s1.str("");
                s1.clear();
                ofstream output (s2);
                

//               
                do
                {
                    synBeam.particle[0].xyz[5]  = pz0;
                    synBeam.particle[0].xyz[4]  = zpos;
                    synBeam.particle[0].xyz[3]  = 0;
                    synBeam.particle[0].xyz[2]  = 0;
                    synBeam.particle[0].xyz[1]  = 0;
                    synBeam.particle[0].xyz[0]  = 0;
                    double time1                = synParticleRun.timeNow;
                    deltatime                   = 0;
                    deltastep                   = 0;
                    double cosvalue             = 0;
                    double sinvalue             = 0;

 
                    do
                    {
                        lattice       .GetExternalField           (time1, synBeam  );
                        synBeam       .SetField                   (                                 );
                        leapfrog      .Update                     (synParticleRun.timeStep, synBeam, lattice);
                        
                        time1            +=    synParticleRun.timeStep;
                        pz1               =    synBeam.particle[0].xyz[5];
                        wou               =   (sqrt(1+pow(pz1,2)) -1) * BaseMassInMeV;
                        deltatime        +=    synParticleRun.timeStep;
                        deltastep        +=    1;
                        
                        double gammatemp  =   sqrt(1+pow(pz1,2));
                        double betatemp   =   sqrt(pow(gammatemp,2)-1)/ gammatemp;


                        cosvalue         +=   synBeam.particle[0].externalElecField[2] * cos(synBeam.cossyn/180.0*M_PI + 2*M_PI*freq*time1 ) * betatemp * C_light *synParticleRun.timeStep ;
                        sinvalue         +=   synBeam.particle[0].externalElecField[2] * sin(synBeam.cossyn/180.0*M_PI + 2*M_PI*freq*time1 ) * betatemp * C_light *synParticleRun.timeStep ;  
                        

//                       test <<setw(15)<<time1
//                            <<setw(15)<<time1   /   synParticleRun.timeStep
//                            <<setw(15)<<synBeam.cossyn
//                            <<setw(15)<<synBeam.particle[0].xyz[0]
//                            <<setw(15)<<synBeam.particle[0].xyz[1]
//                            <<setw(15)<<synBeam.particle[0].xyz[2]
//                            <<setw(15)<<synBeam.particle[0].xyz[3]
//                            <<setw(15)<<synBeam.particle[0].xyz[4]
//                            <<setw(15)<<synBeam.particle[0].xyz[5]
//                            <<setw(15)<<synBeam.particle[0].externalElecField[0]
//                            <<setw(15)<<synBeam.particle[0].externalElecField[1]
//                            <<setw(15)<<synBeam.particle[0].externalElecField[2]
//                            <<endl; 
                    
                    }while(lattice.element[synBeam.particle[0].indexElem]->name=="field");
                    
//                    cout<<"stop here!!!!to check"<<endl;
//                    getchar();
                    
                    output<< setw(15)<<synBeam.cossyn<<setw(15)<<wou<<"\t"<<cosvalue<<"\t"<<sinvalue<<"\t"<<atan2(cosvalue,sinvalue)/M_PI*180  <<endl;
                    
                    pz1        =   synBeam.particle[0].xyz[5];
                    wou        =   (sqrt(1+pow(pz1,2)) -1) * BaseMassInMeV;
                    synBeam.cossyn    =   synBeam.cossyn + 0.01;
                    
                    cout<<wou<<"\t"<<synBeam.cossyn <<"\t"<<desigw[numoffiled]<<endl;
                    
                    if(abs(wou-desigw[numoffiled])<1.0e-4)
                    {
                        synphaseset<<setw(15)<<numoffiled
                                   <<setw(15)<<synBeam.cossyn
                                   <<setw(15)<<wou<<endl;
                         cout<<"this synphase for this field "<<synBeam.cossyn<<"\t"<<count<<endl;  
//                         getchar();        
//                         ++count;
                       
//                       if(count>1) 
//                       {
//                            break;
//                       }
                        
                    }
                    

                }while(abs(synBeam.cossyn-360)>0.1 );               //while(abs(wou-desigw[numoffiled])>10e-5);
                
                
                synParticleRun.timeNow  =   synParticleRun.timeNow + deltatime;
                i  =   i + deltastep;
                cout<<i<<"the "<< numoffiled<<"filed synphase: "<<synBeam.cossyn<<"\t"<<wou<<endl;
                synBeam.particle[0].xyz[5] =sqrt( pow(1 + desigw[numoffiled] / BaseMassInMeV, 2) -1 );
                numoffiled             += 1 ;
                output.close();
 
                
                
//                getchar();
            }
        }
        
        if(synBeam.particle[0].xyz[4]>lattice.positionBack.back())   //beam are out of lattice
        {
            cout<<"lattice is ran out, and the total length is  "<<lattice.positionBack.back()<<"m"<<endl;
            break;
        }
        
    }
    

 
    cout<<"the  syn-phase have been set for cavities"<<endl;
    exit(0);

}











