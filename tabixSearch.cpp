#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
#include "tabixSearch.h"
#include "tabix.h"

#define M 162176005  //total number of pairs 
#define A 19.4814084220388  //sum(1/(1:M))

double dabs(double val)
{
  if (isnan(val))
    return val;
  return (val > 0)? val: (-1.0*val);
}

tabSearch::tabSearch(const char* gzfile, const char* gzidx, int readcnt, int loadnow)
{
  ifstream f1(gzfile);
  ifstream f2(gzidx);  
  if (!f1.good())
  {
    std::cout << "Error, "<< gzfile << "cannot open or not exists\n";
    return;
  }
  else if (!f2.good())
  {
    std::cout << "Error, " << gzidx << "cannot open or not exists\n";
    return;
  }
  f1.close();
  f2.close();
  //line=new char[256];
  totalreads=readcnt;
  gz=gzfile;
  index=gzidx;
  tab=ti_open(gz,index);
  //idxload();
  if(loadnow)
  {  
    if(idxload())
      std::cout<< "error, can not load the index of "<< gzfile << " and "<< gzidx << "\n";
  }
}

int tabSearch::idxload()
{
  return ti_lazy_index_load(tab); 
}

ti_iter_t tabSearch::search(const char* chrome,int start, int end)
{

  ti_iter_t iter;
  qStart=start;
  qEnd=end;
  //qRawp=pval;
  
  strncpy(qChrome,chrome,strlen(chrome));
  qChrome[strlen(chrome)]='\0';
  iter=ti_query(tab,chrome, start, end);
  return iter;
  
}

/* 
Bonferroni.adj[i]= M*rawp[i];
Holm.adj[i]=max(adj[i-1],min((M-i+1)*rawp[i],1)) or min(M*rawp[1],1)
BH.adjp[i]=min(adj[i+1],min((M/i)*rawp[i],1,na.rm=TRUE),na.rm=TRUE) or 
BY.adjp[i]=min(adj[i+1],min((M*A/i)*rawp[i],1,na.rm=TRUE),na.rm=TRUE)
chr, st, ed, rank, DRM","GeneId","TSS","corr","rawp","rawp", Bonf, Holm, BY, BH, 
datatype:
begin column(bc) is floor(-log(p))-1
end column (ec) is ceiling(-log(p))+1
ec+1 column is the rawp
ec+6,7 rawp, 8,9,10,11,bon,holm, bh,by

 */
void tabSearch::coverage_adjustp(double pval, sk_conf_t datatype,vector<double> &adjp,ti_iter_t iter)
{
  int count=0;
  if (iter ==NULL)// some errors appear, usually search with an unknow chrom
    return;
  //vector<double> val1(10,0.0);
    const char* line=new char[256];
  //char line[256];
  int len=256;
  double  tmprawp=0.0, tmpadjHolm=0.0,adjBH1=1.0, adjBY1=1.0; //rawp and adjp for i+1 rank
  int rank=0;
  
  while(line=ti_read(tab, iter, &len))
  {
    string tmp;
    vector<string> cols;
    //istringstream ss(string(line,strlen(line)));
    istringstream ss(line);
    while(getline(ss, tmp, '\t'))
      cols.push_back(tmp);

    string seqid=cols[datatype.sc];
    int sst=atoi(cols[datatype.bc].c_str()); //current 0-based and convert later
    int sed=atoi(cols[datatype.ec].c_str());
    double rawp=atof(cols[datatype.ec+2].c_str());
    adjBH1=atof(cols[datatype.ec+3].c_str());//i+1 BH1
    adjBY1=atof(cols[datatype.ec+4].c_str());//i+1 BY1
    tmpadjHolm=atof(cols[datatype.ec+5].c_str());
    double adjBonf=atof(cols[datatype.ec+6].c_str());
    double adjHolm=atof(cols[datatype.ec+7].c_str());
    double adjBH=atof(cols[datatype.ec+8].c_str());
    double adjBY=atof(cols[datatype.ec+9].c_str());
    
    rank=atoi(cols[datatype.ec+1].c_str());
    if (count==0)
      tmprawp=rawp;

    if (rank==1 && pval < rawp)
    {
      adjp[0]=min(1.0,M*pval);
      adjp[1]=min(1.0,M*pval);
      adjp[2]=min(adjBH1, min(M*pval,1.0));
      adjp[3]=min(adjBY1, min((M*A)*pval,1.0)); 
      break;
    }
    
    if (rawp==pval)
    { //assign the largest same adjust-p value to adjp
      adjp[0]=adjBonf;
      adjp[1]=adjHolm;
      adjp[2]=adjBH;
      adjp[3]=adjBY;
      break;
    }else if (pval < tmprawp && pval >rawp) //
    {
      adjp[0]=min(1.0,M*pval);
      adjp[1]=max(adjHolm,min((M-(rank+1)+1)*pval,1.0));
      adjp[2]=min(adjBH1, min((M/(rank+1))*pval,1.0));
      adjp[3]=min(adjBY1, min((M*A/(rank+1))*pval,1.0)); 
      break;
    }else if (pval > tmprawp)//seems if only add 1 
    {
      adjp[0]=min(1.0, M*pval);
      adjp[1]=max(adjHolm,min((M-(rank+1)+1)*pval,1.0));
      adjp[2]=min(adjBH1, min((M/(rank+1))*pval,1.0));
      adjp[3]=min(adjBY1, min((M*A/(rank+1))*pval,1.0)); 
      break;
    }
    
    tmprawp=rawp;//retein the previous one
    count++;
  }

  if (pval < tmprawp && rank >1)//recalculate adjp[1]
  {
    adjp[0]=min(1.0, M*pval);
    adjp[1]=max(tmpadjHolm,min((M-(rank)+1)*pval,1.0));
    adjp[2]=min(adjBH1, min((M/(rank))*pval,1.0));
    adjp[3]=min(adjBY1, min((M*A/(rank))*pval,1.0));
  }
  
}

/*
sk_conf_t.present: 0 for methylation, 1 for h3k4me1 and 2 for h3k27ac
need math.h: NAN
 */
double tabSearch::coverage(sk_conf_t datatype,ti_iter_t iter)
{
  double val1=0.0, val2=0.0;
  double tmprawp=1.0;
  const char* line=new char[256];
  int len = 256;
  while(line=ti_read(tab, iter, &len))
  {
    string tmp;
    vector<string> cols;
    istringstream ss(string(line,strlen(line)));
    while(getline(ss, tmp, '\t'))
      cols.push_back(tmp);

    string seqid=cols[datatype.sc];
    int sst=atoi(cols[datatype.bc].c_str()); //current 0-based and convert later
    int sed=atoi(cols[datatype.ec].c_str());
    double v1=0.0,v2=0.0;
    if (seqid != string(qChrome))
    {
      std::cout<< "Error: " << seqid << " != " << qChrome << "\n";
      return NAN;
    }
    
    if (datatype.present==0)//0 meth, 1, h3k4me1, 2,h3k27ac 
    {
      v1=atof(cols[9].c_str());//readcount
      v2=atof(cols[10].c_str())*v1/100;//methyl-ed readcnt
      val1+=v1;
      val2+=v2;
    }
    else if(datatype.present==1 || datatype.present==2)
    {
      sst+=1;//convert to 1-based
      int readlen=sed-sst+1;
      if (qStart <= sst && qEnd >=sst)
      {
        val1+= (qEnd >= sed)?1.0 : ((qEnd - sst +1) * 1.0 /readlen);//readcnt set as 1.0
        continue;
      }else if (qStart > sst && qStart <= sed)
      {
        val1 += (qEnd <= sed)? ((qEnd - qStart +1)*1.0/readlen): ((sed - qStart +1)*1.0/readlen);
        
      }  
    }   
  }
  
  if (datatype.present==0)
    return (val1==0)?NAN:(val2/val1); //math.h
  else if (datatype.present==1 || datatype.present==2)
    return (val1*pow(10,6)/totalreads);
}

void tabSearch::unloadidx()
{
  ti_close(tab);
}
