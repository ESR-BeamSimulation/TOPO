#include<fstream>
#include<vector>
#include<string>
#include <iostream>
#include <iomanip>

#include <iterator>
#include <cstdio> 

#include <cmath>
#include <strstream>
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




int main()
{
 ostringstream s1;
// istringstream stream1;
//  string string1 = "25";
//  stream1.str(string1);
//  int i;
//  stream1 >> i;
//  cout << i << endl;  // displays 25



//    ifstream fin("buncher.dat");

//  string                str;
//  vector<string>        strVec;

//  double ke=1;
//  double kb=1;
//  double position=0;


//    vector<double>      xyzTemp(3);
//    vector<double>      eFieldTemp(3);
//    vector<double>      mFieldTemp(3);

//    for(int i = 0;fin.peek()!=EOF;++i)
//    {
//      getline(fin,str);
//      cout<<str<<endl;
//      StringSplit(str,strVec);
//      xyzTemp[0]        =       strtod(   strVec[0])      /1000;          //mm to m
////      xyzTemp[1]        =       stod(   strVec[1])      /1000;
////      xyzTemp[2]        =       stod(   strVec[2])      /1000   +position;
////      gridLocation.push_back(xyzTemp);

////      eFieldTemp[0]=    stod(   strVec[3])*ke   /1000   ;
////      eFieldTemp[1]=    stod(   strVec[4])*ke   /1000;
////      eFieldTemp[2]=    stod(   strVec[5])*ke   /1000;
////      eField.push_back(eFieldTemp);

////      mFieldTemp[0]=    stod(   strVec[6])*kb   /1000;
////      mFieldTemp[1]=    stod(   strVec[7])*kb   /1000;
////      mFieldTemp[2]=    stod(   strVec[8])*kb   /1000;
////      mField.push_back(mFieldTemp);
//getchar();
//    }
//    fin.close();






return 0;
}





