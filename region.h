#include <iostream>
#include <map>
#include <vector>
#include <string>
#include "segment.h"
using namespace std;

#ifndef SK_REGION_H
#define SK_REGION_H

class Region
{
public:  
 map<string,Segment*> rgn;
 virtual void read()=0;
 virtual void add(const char* chr,int start, const char* strand, double readcnt)=0;
 
 
   
};


class dnase:public Region
{
 public:

   dnase(string fn)
   {  //rgn.name="DNase";
      fname=fn;
   }
   void read();
   void add(const char* chr,int start, const char* strand, double readcnt){};
   virtual ~dnase(){};
 private:
  string fname;
};


class lad:public Region
{
 public:
  //map<string, Segment> rgn;
  lad(string fn)
  {   //rgn.name="LAD";
      fname=fn;
   }

   	void read();
    void add(const char* chr,int start, const char* strand, double readcnt){};

  private:
	string fname;
};

class cpgi:public Region
{
 public:
  map<string, vector<int> > cpgNum;
  map<string, vector<int> > gcNum;
  map<string, vector<float> > obsExp;
  cpgi(string fn)
  {
    //rgn.name="CpGIsland";
    fname=fn;
  }
  void read();
  void add(const char* chr,int start, const char* strand, double readcnt){};

private:
  string fname;
};


#endif
