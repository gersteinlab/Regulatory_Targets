#include <fstream>
#include <sstream>
#include <map>
#include <vector>
/*#if __GNUC__ ==4 && __GNUC_MINOR__ > 5
//#  include <features.h>
//#  if __GNUC_PREREQ(4,6)
#  include <regex> //c++0x
#else*/
#  include <boost/regex.hpp>
//using namespace boost;
//#endif
//#endif
#include "gff.h"
#include "segment.h"

gff::gff(string fn)
{
  fname=fn;
}

void gff::read()
{
  if(fname.empty())
    return;
  ifstream fs;
  fs.open(fname.c_str());
  long lid=0;
  string line;
  while(getline(fs,line))
  {
    /* 
#if __GNUC__ ==4 && __GNUC_MINOR__ >5 
    regex skip("^#.+$");
    smatch sm;
    if(regex_search(line,sm,skip))
    {
      cout<<line<<"\n";
      continue;
    }
#else
    */
    boost::regex skip("^#.+$");
    boost::smatch sm;
    if(boost::regex_search(line,sm,skip))
    {
      cout<<line<<"\n";
      continue;
    }
    //#endif
    istringstream iss(line);
    string tmp;
    lid++;
    
    vector<string> cols;
    while(getline(iss,tmp,'\t'))
    {
      cols.push_back(tmp);
    }
    //
    stringstream ss;
    ss<<"gff"<<lid;
    string sid=ss.str();
    string chrom=cols[0];
    string txid("");
    //comment=cols[8];
    //#if __GNUC_PREREQ(4,6)
    /*#if __GNUC__ == 4 && __GNUC_MINOR__ > 5
    regex ee("transcript_id\\s+\"([^\\.]+)\\.\\d+\";.+gene_name\\s\"([^\"]+)\";");
    //boost::regex ee("transcript_id\\s+\"([^\"]+)\";.+gene_name\\s\"([^\"]+)\";");
    
    smatch matches;
    
    if(regex_search(cols[8],matches,ee))
      {
        txid=string(matches[1]).append("\t").append(matches[2]).append("\t").append(cols[2]).append(":").append(cols[6]);
      }
      
      #else*/
    boost::regex ee("transcript_id\\s+\"([^\\.]+)\\.\\d+\";.+gene_name\\s\"([^\"]+)\";");
    //boost::regex ee("transcript_id\\s+\"([^\"]+)\";.+gene_name\\s\"([^\"]+)\";");
    boost::smatch matches;
    if(boost::regex_search(cols[8],matches,ee))
      {
        txid=string(matches[1]).append("\t").append(matches[2]).append("\t").append(cols[2]).append(":").append(cols[6]);
      }
    //#endif
    map<string,Segment*>::iterator it=rgn.find(chrom);
    if(it!=rgn.end())
    {
      it->second->chromStart.push_back(atoi(cols[3].c_str()));
      it->second->chromEnd.push_back(atoi(cols[4].c_str()));
      it->second->id.push_back(txid);
      it->second->strand.push_back(cols[6]);
    }else
    {
      Segment* sgm=new Segment();
      sgm->chrom=chrom;
      sgm->chromStart.push_back(atoi(cols[3].c_str()));
      sgm->chromEnd.push_back(atoi(cols[4].c_str()));
      sgm->id.push_back(txid);
      sgm->strand.push_back(cols[6]);
      sgm->name="Gencode";
      rgn.insert(pair<string,Segment*>(chrom,sgm));
      
    }
      
    

    
  }
  
    fs.close();
    
}


