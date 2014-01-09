
#include <stdlib.h>
#include <sstream>
#include <fstream>
#include "region.h"
using namespace std;
/*
All the file type here are special bed file, with/without bin-number
so I will change all start to 1-based by adding 1 to them. 2013.10.30

 */
//here is binary,need to check
void dnase::read()
{
  if (fname.empty())
  {
    return;
    
  }
  ifstream fs;
  fs.open(fname.c_str());
  long lid=0;
  string line;
  
  
  while(getline(fs,line))
  {
    istringstream is(line);
    vector<string> cols;
    lid++;
    
    string tmp;
    while(getline(is,tmp,'\t'))
    {
      cols.push_back(tmp);
    }
    
    stringstream ss;
    ss<<"dnase"<<lid;
    string sid=ss.str();
    
    string chrom;
    chrom.append(cols[0]);
    map<string,Segment*>::iterator it=rgn.find(chrom);
    
    if (it != rgn.end())
    {
      //int pos=atol(cols[(std:size_t)2]);
      
      it->second->chromStart.push_back(atoi(cols[(size_t)1].c_str())+1);
      it->second->chromEnd.push_back(atoi(cols[(size_t)2].c_str()));
      it->second->id.push_back(sid);
      it->second->value.push_back(atoi(cols[(size_t)6].c_str()));
      //it->second->chromStart.push_back(atol(cols[(std:size_t)2]))
    }
    else
    {
      Segment* sgm=new Segment();
      sgm->chrom=chrom;
      sgm->chromStart.push_back(atol(cols[(size_t)1].c_str())+1);
      sgm->chromEnd.push_back(atol(cols[(size_t)2].c_str()));
      sgm->value.push_back(atoi(cols[(size_t)6].c_str()));
      sgm->name="DNase";
      sgm->id.push_back(sid);
      rgn.insert(pair<string,Segment*>(chrom,sgm));
    }
    
  }
  
  /*
  fs.open(fname.c_str(),iostream::in|iostream::binary);
  string chrom;
  uint chromStart;
  uint chromEnd;
  string name;
  uint score;
  char strand[1];
  float sigVal;
  float pVal;
  float qVal;
  int peak;
  while(fs.good())
  {
    
    fs.read((char*)&chrom,sizeof(string));
    fs.read((char*)&chromStart,sizeof(uint));
    fs.read((char*)&chromEnd,sizeof(uint));
    fs.read((char*)&name,sizeof(string));
    fs.read((char*)&score,sizeof(uint));
    fs.read(strand,sizeof(char)*1);
    fs.read((char*)&sigVal,sizeof(float));
    fs.read((char*)&pVal,sizeof(float));
    fs.read((char*)&qVal,sizeof(float));
    fs.read((char*)&peak,sizeof(int));

    map<string,Segment*>::iterator it;
    it=rgn.find(chrom);
    if(it!=rgn.end())
    {
      it->second->chromStart.push_back(chromStart);
      it->second->chromEnd.push_back(chromEnd);
      //it->second->name.push_back("DNase");
    }else
    {
      Segment* sgm=new Segment();
      sgm->chrom=chrom;
      sgm->chromStart.push_back(chromStart);
      sgm->chromEnd.push_back(chromEnd);
      sgm->name="DNase";
      rgn.insert(pair<string,Segment*>(chrom,sgm));
      
    }
    

  }
  */
  fs.close();
  
}


void lad::read()
{
  fstream fs;
  fs.open(fname.c_str(),iostream::in);
  string line;
  //getline(fs,line);//remove header line
  long lid=0;
  
  while(getline(fs,line))
  {
    istringstream is(line);
    string tmp;
    lid++;
    
    int i=0;
    string chrom;
    int start,end;
    int bin;
    
    while(getline(is,tmp,'\t'))
    {
      switch (i){
      case 0://bin
        bin=atoi(tmp.c_str());
        
        break;
      case 1://chr
        chrom=tmp;
        
        break;
      case 2://st
        start=atol(tmp.c_str());
        
        break;
      case 3://end
        end=atol(tmp.c_str());
        
        break;
      default:
        
        break;
      }
      i++;
    }

    stringstream ss;
    ss<<"lad"<<lid;
    string sid=ss.str();
    
    map<string,Segment*>::iterator it;
    it=rgn.find(chrom);
    if (it!=rgn.end())
    {
      it->second->chromStart.push_back(start+1);
      it->second->chromEnd.push_back(end);
      it->second->bin.push_back(bin);
      it->second->id.push_back(sid);
      
      
    }else//new segment
    { Segment* sgm=new Segment();
      sgm->chrom=chrom;
      sgm->chromStart.push_back(start+1);
      sgm->chromEnd.push_back(end);
      sgm->bin.push_back(bin);
      sgm->id.push_back(sid);
      
      sgm->name="LADs";
      
      rgn.insert(pair<string,Segment*>(chrom,sgm));
      
    }
    
    
    
  }
  
  fs.close();
  
}

/*

 */
void cpgi::read()
{

  if(fname.empty())
  {
    return;
  }

  fstream fs;
  fs.open(fname.c_str(),istream::in);
  string line;
  //if have header,remove
  //getline(fs,line);
  long lid=0;
  
  while(getline(fs,line))
  {
    istringstream is(line);
    vector<string> cols;
    lid++;
    
    string tmp;
    while(getline(is,tmp,'\t'))
    {
      cols.push_back(tmp);
    }
    //read to rgn
    stringstream ss;
    ss<<"cpgi"<<lid;
    string sid=ss.str();
    
    string chrom;
    chrom.append(cols[1]);
    map<string,Segment*>::iterator it=rgn.find(chrom);
    
    if (it != rgn.end())
    {
      //int pos=atol(cols[(std:size_t)2]);
      
      it->second->chromStart.push_back(atol(cols[(size_t)2].c_str())+1);
      it->second->chromEnd.push_back(atol(cols[(size_t)3].c_str()));
      it->second->id.push_back(sid);
      
      //it->second->chromStart.push_back(atol(cols[(std:size_t)2]))
    }
    else
    {
      Segment* sgm=new Segment();
      sgm->chrom=chrom;
      sgm->chromStart.push_back(atol(cols[(size_t)2].c_str())+1);
      sgm->chromEnd.push_back(atol(cols[(size_t)3].c_str()));
      sgm->name="CpGIsland";
      sgm->id.push_back(sid);
      rgn.insert(pair<string,Segment*>(chrom,sgm));
    }
    }
  
  

  fs.close();
  

}
