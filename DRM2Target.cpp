/*
Author: Shaoke LOU
Date:   2013-11-28



*/


#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <algorithm>
#include "tabixSearch.h"
#include "drmPair.h"
#include "vmstat.h"
#include "tabix.h"
#include "gff.h"
#include "stat.h"
#include "adjustp.h"
#include <cmath>

//#define VERSION "1.02"
#define DEBUG

using namespace std;
class Fragment
{
 public:
  string chrom;
  vector<int> chromStart;
  vector<int> chromEnd;
  vector<string> strand;
  vector<string> id;
  vector<vector<double> > levels;//meth,20, h3k4me1,19, h3k27ac,19
};

class GZInfo
{
public:
  vector<int> totalreads;
  vector<string> gzfile;
  string sample;  
};
typedef map<string, GZInfo*> GZInfoMap;
typedef map< string, Fragment*> Range;
typedef struct 
{
  char chrom[50];
  int start;
  int end;
  double pval;
} RAW2ADJ_conf_t;

void updateGZInfo(GZInfoMap& gzinfo, string tissue, int totalreads, string gzfile);
void genRawp2AdjpConf(vector<double> pval,int dtype, int corrtype, int tailtype, RAW2ADJ_conf_t &conf);
void updateDRMcoverage(char* data_dir,sk_conf_t conf,GZInfoMap gzinfo, int st_id, Range &drm);
void usage(char *title);

void corrSigTest(vector<double> &a, vector<double> &b, int corrtype, vector<double> &pvalue);

int main(int argc, char** argv)
{

  //variable
  char *pgz=NULL;
  char *pidx=NULL;//TODO need to declear, pvalue index file
  int maxDist=1000000;//-M
  int minDist=10000;//-m
  char *drmbed=NULL, *annotgtf=NULL,*rpmfile=NULL; //-b, -g, -r
  int MAX_DRM_NUM=3000;//-N
  char* data_dir=NULL; //-d directory for bgzip and index file
  char* output=NULL;// -o 
  char* gzlist=NULL;//-l
  char *cflag=NULL,*tflag=NULL;
  double pv_cutoff=0.05;
  int padj_bool=0;
  int padj_new=0;
  int meth_corrtype=1,h3k4me1_corrtype=1,h3k27ac_corrtype=1;//input a list, -c, 0 for pearson, 1 for spearman
  int meth_tailtype=0,h3k4me1_tailtype=1,h3k27ac_tailtype=1;//input a list, -t, 0 1 2, 0 left, 1 right, 2 both
  int opt;
  
  int opt_b=0, opt_g=0, opt_r=0, opt_d=0,opt_o=0, opt_l=0, opt_c=0, opt_t=0;
  int opt_s=0;
  int tmplen=0;

  while((opt = getopt(argc, argv,"p:jJs:m:M:b:g:r:N:d:l:c:t:o:h"))!=-1)
  {
    switch(opt)
    {
    case 'p':pv_cutoff=atof(optarg);break;
    case 'j':padj_bool=1;break;
    case 'J':padj_new=1;break;
    case 's':pgz=optarg;opt_s=1;break;
    case 'm':minDist=atoi(optarg);break;
    case 'M':maxDist=atoi(optarg);break;
    case 'b':drmbed=optarg;opt_b=1;break;
    case 'g':annotgtf=optarg;opt_g=1;break;
    case 'r':rpmfile=optarg;opt_r=1;break;
    case 'N':MAX_DRM_NUM=atoi(optarg);break;
    case 'd':data_dir=optarg;opt_d=1;break;
    case 'l':gzlist=optarg;opt_l=1;break;
    case 'c':
      cflag=optarg;
      if (strlen(cflag)>=5)
      {
        meth_corrtype=cflag[0]-'0';
        h3k4me1_corrtype=cflag[2]-'0';
        h3k27ac_corrtype=cflag[4]-'0';
      }      
      break;
    case 't':
      tflag=optarg;
      if (strlen(tflag)>=5)
      {
        meth_tailtype=tflag[0]-'0';
        h3k4me1_tailtype=tflag[2]-'0';
        h3k27ac_tailtype=tflag[4]-'0';
      }
      break;
    case 'o':output=optarg;opt_o=1;break;
    case 'h':usage(argv[0]); exit(1);
    default: usage(argv[0]); break;
    }
  }
  //printf("data_dir:%s",data_dir);
  if (!(opt_s && opt_b && opt_g && opt_r && opt_d && opt_o && opt_l ))
  {
    usage(argv[0]);exit(1);
  }
  if (data_dir[strlen(data_dir)-1]-'/'!=0)
  {
    strcat(data_dir,"/");
  }
  
  pv_cutoff=(pv_cutoff > 0.05 && padj_bool==1)? 1.0:pv_cutoff;//if adjusting p-value using pre-calculated data, pv_cutoff <=0.05 or output all results
  padj_bool=(padj_new)?0:padj_bool;//if adjust pvalue using user's data, it will not adjust based on pre-calculated data

  pidx=new char[strlen(pgz)+5];
  strcpy(pidx,pgz);
  strcat(pidx,".tbi");
  
  Range drm, annot;
  GZInfoMap meth_gz, h3k4me1_gz, h3k27ac_gz;//in order of rpm sample
  GZInfoMap::iterator gziter;
  map<string, vector<double> > gencode_rpm;
  //hist_missing[2]={1, 8};//The 1st and 8th tissues doesn't exist in H3K4me1 and H3K27ac respectively, hardcoding here.Don't change the order of colume in rpm data 
  
  time_t start,end;
  int i, j;
  string line;
  string tmp;
  int line_num;
  Range::iterator it1;
  map<string,Segment*>::iterator it2;

  time(&start);
  /*Annotation file
   */
  ofstream outfs;
  outfs.open(output);//print header
  if (padj_bool)
  {  
    outfs << "chrom\tquer_start\tquery_end\tDRM_id\tgene_id\tgene_name\tmeth.type\tmeth.corr\tmeth.raw-pval\tmeth.Bonferroni\tmeth.Holm\tmeth.BH\tmeth.BY\tH3K4me1.type\tH3K4me1.corr\tH3K4me1.raw-pval\tH3K4me1.Bonferroni\tH3K4me1.Holm\tH3K4me1.BH\tH3K4me1.BY\tH3K27ac.type\tH3K27ac.corr\tH3K27ac.raw-pval\tH3K27ac.Bonferroni\tH3K27ac.Holm\tH3K27ac.BH\tH3K27ac.BY\n";
  }else if (padj_new)
  {
    outfs << "chrom\tquer_start\tquery_end\tDRM_id\tgene_id\tgene_name\ttype\tcorr\traw-pval\tBonferroni\tHolm\tBH\tBY\n";
  }
  else
  {
    outfs << "chrom\tquer_start\tquery_end\tDRM_id\tgene_id\tgene_name\tmeth.type\tmeth.corr\tmeth.raw-pval\tH3K4me1.type\tH3K4me1.corr\tH3K4me1.raw-pval\tH3K27ac.type\tH3K27ac.corr\tH3K27ac.raw-pval\n";
  }
  
  

  gff* pAnnot=new gff(string(annotgtf));
  pAnnot->read();

  /*RPM file
   */
  ifstream rpmio;
  rpmio.open(rpmfile);
  while(getline(rpmio,line))
  {
    //Format: gene_id.tss\trpm1..rpm20
    istringstream ss(line);
    vector<double> rpm;
    string gid;
    int cid=0;
    while(getline(ss,tmp,'\t'))
    {
      if (cid==0)
        gid=tmp;
      else
        rpm.push_back(atof(tmp.c_str()));
      cid++;
    }
    gencode_rpm.insert(pair<string, vector<double> >(gid, rpm));
  }
  rpmio.close();
  
  /*open meta information
  format:H3k27ac\t$x\t$linecnt\t$outf.gz
  */
  ifstream metafin;
  metafin.open(gzlist);
  while(getline(metafin,line))
  {
    istringstream ss(line);
    vector<string> cols;
    while(getline(ss,tmp,'\t'))
    {
      cols.push_back(tmp);
    }
    string measures=cols[0];
    string sample=cols[1];
    int totalreads=atoi(cols[2].c_str());
    string gzfile=cols[3];
    if (measures=="Methyl")
      updateGZInfo(meth_gz, sample, totalreads, gzfile);
    else if (measures=="H3K4me1")
      updateGZInfo(h3k4me1_gz, sample, totalreads, gzfile);
    else if (measures=="H3K27ac")
      updateGZInfo(h3k27ac_gz, sample, totalreads, gzfile);
  }  
  metafin.close();

  /*DRM regions and init methylation, H3K4me1 and H3K27ac levels
    Range *drm and input bed 
  */
  line_num=0;
  ifstream drmio;
  drmio.open(drmbed);
  while(getline(drmio,line))
  {
    line_num++;
    istringstream ss(line);
    vector<string> cols;
    while(getline(ss, tmp, '\t'))
    {
      cols.push_back(tmp);
    }
    //only chr,start,end, the id is line_num
    string chrom=cols[0];
    int start=atoi(cols[1].c_str())+1;//convert to 1 based
    int end=atoi(cols[2].c_str());
    //    char id[20];
    //sprintf(id, "%d", line_num);
    //init methylation, histone levels, hardcoding the size of vector<double>=58
    vector<double> drmlevs(58, 0.0);
    it1=drm.find(chrom);
    if (it1!=drm.end())
    {
      it1->second->chrom=chrom;
      it1->second->chromStart.push_back(start);
      it1->second->chromEnd.push_back(end);
      it1->second->id.push_back(cols[3]);
      it1->second->levels.push_back(drmlevs);
    }
    else
    {
      Fragment *sgm=new Fragment();
      sgm->chrom=chrom;
      sgm->chromStart.push_back(start);
      sgm->chromEnd.push_back(end);
      sgm->id.push_back(cols[3]);
      sgm->levels.push_back(drmlevs);
      drm.insert(pair<string, Fragment*> (chrom,sgm));
    }
  }
  drmio.close();
  time(&end);
  cout << "To read the input file and meta information cost "<< (int)((end-start)/60) << " Mins\n";
  
#ifdef DEBUG
  getVM();
#endif

  time(&start);
  sk_conf_t h3k4me1_type={1,0,1,2,0,0};  
  sk_conf_t h3k27ac_type={2,0,1,2,0,0};  
  sk_conf_t meth_type={0,0,1,2,0,0};
  GZInfoMap::iterator striter;
  unsigned  stridx,stridx0,strcnt;
  updateDRMcoverage(data_dir,meth_type, meth_gz, 0, drm);
  updateDRMcoverage(data_dir,h3k4me1_type, h3k4me1_gz, 20, drm);
  updateDRMcoverage(data_dir, h3k27ac_type, h3k27ac_gz, 39, drm);

  time(&end);
  cout << "To calculate the coverage of methylatioin, h3k4me1, and h3k27ac costs "<< (int)((end-start)/60) << " Mins\n"; 
#ifdef DEBUG
  getVM();
#endif

  time(&start);
  //load rawp-adjp table,begin column(bc) is floor(-log(p)),end column (ec) is ceiling(-log(p)),ec+1 column is the rawp
  tabSearch rawp2adjptbl(pgz, pidx, 1, 1);

  //var for new-adjustment
  vector<ADJP> MeAdjpNew,K4AdjpNew, K27AdjpNew;
  map<int, DRMTarPair> drmtarget;
  int idx_new=0;
  //vector<DRMTarPair> drmtarget2;
  map<string, Segment*> regtmp=pAnnot->rgn;
  for(it1=drm.begin(); it1!=drm.end();it1++)
  {
    string chrom=it1->first;
    Fragment* querysgm=it1->second;
    Segment* srcsgm=NULL;
    it2=regtmp.find(chrom);
    if (it2!=regtmp.end())
      srcsgm=it2->second;
    else
      continue; 
    vector<int>::size_type sz_sSt,sz_sEd,sz_qSt,sz_qEd,sz,szs;
    sz_sSt=srcsgm->chromStart.size();
    sz_sEd=srcsgm->chromEnd.size();
    sz_qSt=querysgm->chromStart.size();
    sz_qEd=querysgm->chromEnd.size();
    assert((sz_sSt==sz_sEd) && (sz_qSt==sz_qEd));
    sz=sz_qSt;//use a general one
    szs=sz_sSt;
    //vector<DRMTarPair>* target=new vector<DRMTarPair>();

    for (unsigned i=0; i< sz; i++)
    {
      int qst=querysgm->chromStart[i];
      int qed=querysgm->chromEnd[i];
      string qid=querysgm->id[i];
      int closest_dist=maxDist;
      vector<DRMTarPair> target;      
      int dist0=0;
      for (unsigned j=0; j<szs; j++)
      {
        int sst=srcsgm->chromStart[j];
        int sed=srcsgm->chromEnd[j];
        string strand=srcsgm->strand[j];
        int tss=sst;
        if (strand=="-")
          tss=sed;
        //stringstream ss or to_string
        char sid[15+33];
        int tmpid=0;
        string tmpcol;
        string gname;
        istringstream ss(srcsgm->id[j]);
        while(getline(ss, tmpcol, '\t'))
        {
          if (tmpid==0)sprintf(sid,"%s.%d",tmpcol.c_str(),tss);
          else if (tmpid==1) gname=tmpcol;
          tmpid++;
        }
        //string sid0=srcsgm->id[j].substr(0,15);
        //sid.append(".").append(to_string(tss));//gene_id.tss, no,here should
        DRMTarPair dpair;
        if (qed < sst)
        {
          if(sst - qed < maxDist)
          {
            dist0=sst-qed;
            closest_dist=(closest_dist < dist0)?dist0:closest_dist;
          }else //here no bracket
            break;
        }
        else if (qst > sed)
        {
          if (qst - sed < maxDist)
          {
            dist0=qst-sed;
            closest_dist=(closest_dist < dist0)?dist0:closest_dist;
          }else //here no bracket
            continue;
        }
        else //overlap 
        {
          dist0=0;
          closest_dist=dist0;
          break;
        }
        dpair.chrom=chrom;
        dpair.start=qst;
        dpair.end=qed;
        dpair.drmId=qid;
        dpair.geneId=string(sid);
        dpair.geneName=gname;//add gene name information
        dpair.strand=strand;
        dpair.dist=dist0;
        dpair.closestDist=closest_dist;
        target.push_back(dpair);
      }
      //only consider the DRMs
      if (closest_dist < minDist)
      {
        //target.clear();target.shrink_to_fit();//c++0x
        continue;
      }
      
      vector<double> levels=querysgm->levels[i];
      vector<double> meth(levels.begin(),levels.begin()+20);
      vector<double> h3k4me1(levels.begin()+20,levels.begin()+39);
      vector<double> h3k27ac(levels.begin()+39, levels.end());
      map<string, vector<double> >::iterator giter;
      vector<double> rpm;
      vector<DRMTarPair>::iterator it3;
      for (it3=target.begin(); it3!=target.end();it3++)
      {
        //gene_id
        string gid=it3->geneId;//here
        string gname=it3->geneName;
        giter=gencode_rpm.find(gid);
        if (giter!=gencode_rpm.end())
        {
          rpm=giter->second;
          vector<double> rpm4h3k4me1(rpm.begin()+1,rpm.end());//exclude the first tissue
          vector<double> rpm4h3k27ac;//exclude the 8th tissue(7 for 0-based)
          //rpm4h3k27ac.reserve(rpm.size()-1);
          rpm4h3k27ac.insert(rpm4h3k27ac.end(),rpm.begin(),rpm.begin()+7);
          rpm4h3k27ac.insert(rpm4h3k27ac.end(),rpm.begin()+8,rpm.end());

          vector<double> pvalue_meth(4,0.0),pvalue_h3k4me1(4,0.0),pvalue_h3k27ac(4,0.0);
          corrSigTest(meth,rpm,meth_corrtype, pvalue_meth);
          corrSigTest(h3k27ac,rpm4h3k27ac,h3k27ac_corrtype,pvalue_h3k27ac);
          corrSigTest(h3k4me1, rpm4h3k4me1, h3k4me1_corrtype, pvalue_h3k4me1);

          RAW2ADJ_conf_t meth_conf,h3k4me1_conf, h3k27ac_conf;
          genRawp2AdjpConf(pvalue_meth,0,meth_corrtype,meth_tailtype,meth_conf);
          genRawp2AdjpConf(pvalue_h3k4me1,1,h3k4me1_corrtype,h3k4me1_tailtype,h3k4me1_conf);
          genRawp2AdjpConf(pvalue_h3k27ac,2,h3k27ac_corrtype,h3k27ac_tailtype,h3k27ac_conf);
          
          if (padj_bool)
          {
            vector<double> meth_adjp(4,1.0),h3k4me1_adjp(4,1.0),h3k27ac_adjp(4,1.0);
            ti_iter_t iter;
            if (meth_conf.pval <= pv_cutoff)
            {
              iter=rawp2adjptbl.search(meth_conf.chrom, meth_conf.start, meth_conf.end);//
              rawp2adjptbl.coverage_adjustp(meth_conf.pval, meth_type, meth_adjp,iter);
            }
            if (h3k4me1_conf.pval <= pv_cutoff)
            {
              iter=rawp2adjptbl.search(h3k4me1_conf.chrom, h3k4me1_conf.start, h3k4me1_conf.end);
              rawp2adjptbl.coverage_adjustp(h3k4me1_conf.pval, h3k4me1_type, h3k4me1_adjp,iter);
            }
            if (h3k27ac_conf.pval <= pv_cutoff)
            {
              iter=rawp2adjptbl.search(h3k27ac_conf.chrom, h3k27ac_conf.start, h3k27ac_conf.end);
              rawp2adjptbl.coverage_adjustp(h3k27ac_conf.pval, h3k27ac_type, h3k27ac_adjp,iter);
            }
            
            int adjplog=0;
            for (int adjidx=0; adjidx <4; adjidx++)
            {
              if (h3k4me1_adjp[adjidx] <= pv_cutoff || meth_adjp[adjidx]<= pv_cutoff || h3k27ac_adjp[adjidx] <= pv_cutoff){
                adjplog=1;
                break;
              }
            }
            
            if (adjplog==1)
            {//output
              outfs << chrom << "\t" << qst << "\t" << qed << "\t"<<qid<<"\t" << gid << "\t" << gname << "\t";
              outfs << meth_conf.chrom << "\t" << pvalue_meth[0] << "\t" << meth_conf.pval << "\t" << meth_adjp[0] << "\t" << meth_adjp[1] << "\t" << meth_adjp[2] << "\t" << meth_adjp[3] << "\t";
              outfs << h3k4me1_conf.chrom << "\t" << pvalue_h3k4me1[0] << "\t" << h3k4me1_conf.pval << "\t" << h3k4me1_adjp[0] << "\t" << h3k4me1_adjp[1] << "\t" << h3k4me1_adjp[2] << "\t" << h3k4me1_adjp[3]<<"\t";
              outfs << h3k27ac_conf.chrom << "\t" << pvalue_h3k27ac[0] << "\t" << h3k27ac_conf.pval << "\t" << h3k27ac_adjp[0] << "\t" << h3k27ac_adjp[1] << "\t" << h3k27ac_adjp[2] << "\t" << h3k27ac_adjp[3]<<"\n";
            }
          }else if (padj_new==1)
          {
            if ( std::isnan(meth_conf.pval) && std::isnan(h3k4me1_conf.pval) && std::isnan(h3k27ac_conf.pval))
              continue;
            
            //drmtarget.insert(pair<int, DRMTarPair*>(idx_new, target.data()+(it3-target.begin())));
            DRMTarPair pair0;
            pair0.chrom=it3->chrom;
            pair0.start=it3->start;
            pair0.end=it3->end;
            pair0.drmId=it3->drmId;
            pair0.geneId=it3->geneId;
            pair0.geneName=it3->geneName;//add gene name information
            pair0.strand=it3->strand;
            pair0.dist=it3->dist;
            pair0.closestDist=it3->closestDist;
            drmtarget.insert(pair<int, DRMTarPair>(idx_new, pair0));
            
            if (!std::isnan(meth_conf.pval))
            {
              ADJP tmpadjp(idx_new, pvalue_meth[0], meth_conf.pval, meth_conf.chrom);
              MeAdjpNew.push_back(tmpadjp);
            }
            if (!std::isnan(h3k4me1_conf.pval))
            {
              ADJP tmpadjp(idx_new, pvalue_h3k4me1[0], h3k4me1_conf.pval, h3k4me1_conf.chrom);
              K4AdjpNew.push_back(tmpadjp);
            }

            if (!std::isnan(h3k27ac_conf.pval))
            {
              ADJP tmpadjp(idx_new, pvalue_h3k27ac[0], h3k27ac_conf.pval, h3k27ac_conf.chrom);
              K27AdjpNew.push_back(tmpadjp);
            }
            idx_new++;
          }else
          {
            if (meth_conf.pval <= pv_cutoff || h3k4me1_conf.pval <= pv_cutoff || h3k27ac_conf.pval <= pv_cutoff)
            {
              outfs << chrom << "\t" << qst << "\t" << qed << "\t"<<qid<<"\t" << gid << "\t" << gname << "\t";   
              outfs << meth_conf.chrom << "\t" << pvalue_meth[0] << "\t" << meth_conf.pval << "\t";
              outfs<< h3k4me1_conf.chrom << "\t" << pvalue_h3k4me1[0] << "\t" << h3k4me1_conf.pval << "\t" ;
              outfs << h3k27ac_conf.chrom << "\t" << pvalue_h3k27ac[0] << "\t" << h3k27ac_conf.pval << "\n";
            }
          }
        }
      }
    }
  }

  if (padj_new)
  {
    sort(MeAdjpNew.begin(), MeAdjpNew.end());
    sort(K4AdjpNew.begin(), K4AdjpNew.end());
    sort(K27AdjpNew.begin(), K27AdjpNew.end());
  
    double** MeAdjp=adjust(MeAdjpNew);
    double** K4Adjp=adjust(K4AdjpNew);
    double** K27Adjp=adjust(K27AdjpNew);
  
    //output
     
    int AdjpSize=MeAdjpNew.size();
    for (int i=0; i< AdjpSize; i++)
    {
      if (MeAdjp[i][0] > pv_cutoff && MeAdjp[i][1] > pv_cutoff 
          && MeAdjp[i][2] > pv_cutoff && MeAdjp[i][3] > pv_cutoff)
        continue;
      int tmpidx=MeAdjpNew[i].idx;
      DRMTarPair tmpdrm=drmtarget[tmpidx];

      outfs << tmpdrm.chrom <<"\t" << tmpdrm.start <<"\t" << tmpdrm.end << "\t";
      outfs << tmpdrm.drmId <<"\t" << tmpdrm.geneId << "\t" << tmpdrm.geneName <<"\t";
      outfs << MeAdjpNew[i].type <<"\t" << MeAdjpNew[i].corr<< "\t" << MeAdjpNew[i].rawp <<"\t" ;
      outfs << MeAdjp[i][0] << "\t" << MeAdjp[i][1] << "\t" << MeAdjp[i][2] <<"\t" << MeAdjp[i][3] <<"\n";
      //outfs.flush();
    }
    AdjpSize=K4AdjpNew.size();
    for (int i=0; i< AdjpSize; i++)
    {
      if (K4Adjp[i][0] > pv_cutoff && K4Adjp[i][1] > pv_cutoff && K4Adjp[i][2] > pv_cutoff && K4Adjp[i][3] > pv_cutoff)
        continue;
      int tmpidx=K4AdjpNew[i].idx;
      DRMTarPair tmpdrm=drmtarget[tmpidx];

      outfs << tmpdrm.chrom <<"\t" << tmpdrm.start <<"\t" << tmpdrm.end << "\t";
      outfs << tmpdrm.drmId <<"\t" << tmpdrm.geneId << "\t" << tmpdrm.geneName <<"\t";
      outfs << K4AdjpNew[i].type << "\t" << K4AdjpNew[i].corr<< "\t" << K4AdjpNew[i].rawp <<"\t" ;
      outfs << K4Adjp[i][0] << "\t" << K4Adjp[i][1] << "\t" << K4Adjp[i][2] <<"\t" << K4Adjp[i][3] <<"\n";
    }

    AdjpSize=K27AdjpNew.size();
    for (int i=0; i< AdjpSize; i++)
    {
      if (K27Adjp[i][0] > pv_cutoff && K27Adjp[i][1] > pv_cutoff && K27Adjp[i][2] > pv_cutoff && K27Adjp[i][3] > pv_cutoff)
        continue;
      int tmpidx=K27AdjpNew[i].idx;
      DRMTarPair tmpdrm=drmtarget[tmpidx];

      outfs << tmpdrm.chrom <<"\t" << tmpdrm.start <<"\t" << tmpdrm.end << "\t";
      outfs << tmpdrm.drmId <<"\t" << tmpdrm.geneId << "\t" << tmpdrm.geneName <<"\t";
      outfs <<K27AdjpNew[i].type << "\t" << K27AdjpNew[i].corr<< "\t" << K27AdjpNew[i].rawp <<"\t" ;
      outfs << K27Adjp[i][0] << "\t" << K27Adjp[i][1] << "\t" << K27Adjp[i][2] <<"\t" << K27Adjp[i][3] <<"\n";
    }
  }
  
  outfs.close();
  rawp2adjptbl.unloadidx();
  time(&end);
  cout << "To calculate correlatioins and statistical significance costs "<< (int)((end-start)/60) << " Mins\n";
  
#ifdef DEBUG
  getVM();
#endif
  return 0;
  
}


void updateDRMcoverage(char* data_dir,sk_conf_t conf,GZInfoMap gzinfo, int st_id, Range &drm)
{
  GZInfoMap::iterator striter;
  Range::iterator it1;
  unsigned stridx,stridx0,strcnt;
  strcnt=0;//meth 0-19, h3k4me1 20-38, h3k27ac 39-57
  for(striter=gzinfo.begin();striter!=gzinfo.end();striter++)
  { //start from each gzfile
    //sk_conf_t meth_conf={0,0,1,2,0,0};  
    string sample=striter->first;
    int sampleNum=striter->second->totalreads.size();
    vector<tabSearch*> ttss;
    for (stridx=0; stridx < sampleNum; stridx++)
    {
      string gzfile=string(data_dir).append(striter->second->gzfile[stridx]);
      string gzidx=gzfile;
      gzidx.append(".tbi");
      int totalreads=striter->second->totalreads[stridx]; 
      //TODO:test file's existance
      tabSearch* ts=new tabSearch(gzfile.c_str(),gzidx.c_str(),totalreads,1);
      ttss.push_back(ts);
    }
    
    for (it1=drm.begin(); it1!=drm.end();it1++)
    { //for each chrome
      string chrom=it1->first;
      unsigned drmsize=it1->second->chromStart.size();
      for (unsigned ssz=0; ssz < drmsize; ssz++)
      { //for each region
        double tmpval[sampleNum];
        double tmpval0=0.0;
        int averNum=sampleNum;
        
        int start=it1->second->chromStart[ssz];
        int end=it1->second->chromEnd[ssz];
        for (stridx0=0;stridx0< sampleNum; stridx0++)
        {//for each sample
          ti_iter_t iter=ttss[stridx0]->search(chrom.c_str(), start, end);
          tmpval[stridx0]=ttss[stridx0]->coverage(conf,iter);
          /*if (isnan(tmpval[stridx0]))
          { 
            averNum--;
            continue;
            }*/
          tmpval0+=(std::isnan(tmpval[stridx0]))?0.0:tmpval[stridx0];//here need to check whether nan
        }
        tmpval0=(averNum>0)?tmpval0/averNum:NAN;//math.h NAN()
        it1->second->levels[ssz][strcnt+st_id]=tmpval0;
      }
    }
  
    for (stridx0=0; stridx0< sampleNum;stridx0++)
    {
      ttss[stridx0]->unloadidx();
      delete ttss[stridx0];
    }
    //ttss.clear();
    //ttss.shrink_to_fit();
    strcnt++;
  }
  /*  
  for (it1=drm.begin(); it1!=drm.end();it1++)
  {//chrom
    
    for (unsigned ssz=0; ssz < it1->second->chromStart.size(); ssz++)
    {//region
      cout << conf.present <<"\t"<<  it1->first<<"\t";
    
      cout << it1->second->chromStart[ssz]<<"\t" << it1->second->chromEnd[ssz] <<"\t";
      cout << it1->second->id[ssz]<<"\t";
      strcnt=0;
      for(striter=gzinfo.begin();striter!=gzinfo.end();striter++)
      {
        cout << it1->second->levels[ssz][strcnt+st_id] << "\t";
        strcnt++;
      }
      cout << strcnt <<"\n";
    }
  }
  */
}



void updateGZInfo(GZInfoMap& gzinfo, string tissue, int totalreads, string gzfile)
{
  GZInfoMap::iterator it=gzinfo.find(tissue);
  if (it!=gzinfo.end())
  {
    it->second->sample=tissue;
    it->second->totalreads.push_back(totalreads);
    it->second->gzfile.push_back(gzfile);
  }else
  {
    GZInfo *gzi=new GZInfo();
    gzi->sample=tissue;
    gzi->totalreads.push_back(totalreads);
    gzi->gzfile.push_back(gzfile);
    gzinfo.insert(pair<string, GZInfo*>(tissue, gzi));
  }
}


void genRawp2AdjpConf(vector<double> pval,int dtype, int corrtype, int tailtype, RAW2ADJ_conf_t &conf)
{
  int idx=0;  
  switch(dtype)//meth
  {
  case 0:strcpy(conf.chrom,"Me");break;
  case 1:strcpy(conf.chrom,"K4");break;
  case 2:strcpy(conf.chrom,"K27");break;
    //default:strcpy(conf.chrom,"Me");break;
  }
  switch(corrtype)
  {//change idx=idx+0 for case0 and idx=idx+4 for case1
  case 0:strcat(conf.chrom,"_ps");break;
  case 1:strcat(conf.chrom,"_sp");break;
    //default:idx=0;strcat(conf.chrom,"_ps");break;
  }
  switch(tailtype)
  {
  case 0:idx=idx+2;strcat(conf.chrom,"L");break;
  case 1:idx=idx+3;strcat(conf.chrom,"R");break;
  case 2:idx=idx+1;strcat(conf.chrom,"B");break;
    //default:idx=idx+2;strcat(conf.chrom,"L");break;
  }
  int fval=int(floor(-log(pval[idx])));
  conf.start=(fval>=0)?fval:0;//extend the range so that we can get the rank on the boundary
  conf.end=int(ceil(-log(pval[idx])));//extend the range
  conf.pval=pval[idx];
}

void corrSigTest(vector<double> &a, vector<double> &b, int corrtype, vector<double> &pvalue)
{
  if (corrtype==0)
    pvalue[0]=gsl_pearson_corr(a,b,pvalue[1],pvalue[2],pvalue[3]);//both, left, right
  else if(corrtype==1)
    pvalue[0]=spearman(a,b,pvalue[1],pvalue[2],pvalue[3]);
}


void usage(char *title)
{                                                                                                     
  cerr << title << ": version " << VERSION << endl;                                
  cerr << "Usage: " << title << " <options>" << endl;                                                                
  cerr << "\tVersion: "<<VERSION << endl;
  cerr << '\t' << "-b [string:required]\t"<<"Bed file for DRM regions" << endl;
  cerr << '\t' << "-g [string:required]\t"<<"GTF file for gencode_v7 tss based annotation ../input_data/opt_g.gencode.v7.tss.gff" << endl;
  cerr << '\t' << "-r [string:required]\t"<<"RPM file in gff format ../input_data/opt_r.gencode.v7.tss.rpm.tsv" << endl;   
  cerr << '\t' << "-s [string:required]\t"<<"Correlation significance test file is: opt_s.bed.gz ../input_data/opt_s.bed.gz"<<endl;
  cerr << '\t' << "-d [string:required]\t"<<"Directory for compressed reads alignment file for methylation, H3K4me1, and H3K27ac ../bgzips/" << endl;                               
  cerr << '\t' << "-l [string:required]\t"<<"Meta information file for aligned data file ../input_data/opt_l.drm.meta" << endl;
  cerr << '\t' << "-m [int:10000]      \t"<<"Minimum distance between DRM and its closest gene" << endl;
  cerr << '\t' << "-M [int:1000000]    \t"<<"Maximum distance between DRM and genes," << endl;
  cerr << '\t' << "-p [int:0.05]       \t"<<"pvalue cutoff, should be not greater than 0.05, and if '-j' is set, it is used as the adjusted pvalue cutoff for any of four adjustment methods: Bonferroni, Holm, BH and BY. Tips: if you want to output all results, please set '-p' greater than 0.05" << endl; 
  cerr << '\t' << "-j                  \t"<<"Whether p-value is adjusted based on pre-calculated data. Default: not adjusted." << endl;
  cerr << '\t' << "-J                  \t"<<"To enable '-J' will disable '-j'.  Whether p-value is adjusted within user's data. Default: not adjusted." << endl;
  cerr << '\t' << "-N [int: 3000]      \t"<<"The maximum number of DRM regions can be suggested to run with this program. This tool only works for a small number of regions in a batch" << endl;
  cerr << '\t' << "-c [string:1,1,1]   \t"<<"Correlation methods used for methylation, H3K4me1 and H3K27ac, string delimeted by ',' with 0 and 1.[0, Pearson; 1,Spearman]" << endl;                                   
    cerr << '\t' << "-t [string:0,1,1]   \t"<<"Tailtype for correlation significance test for methylation, H3K4me1 and H3K27ac. [0,left-tail; 1,right-tail; 2,both-tails]. When using '-j' option, use default only. Other combinations are not available for adjustment with pre-calculated data. If you want to get the raw p-value or adjustment within user's data, please don't enable '-j' option! " << endl;                             
  cerr << '\t' << "-o [string:required]\t"<<"output file" << endl;

  //  cerr << '\t' << "-1 [string:optional]\t"<<"p-value adjust methods for methylation" << endl;                    
  //cerr << '\t' << "-2 [int]\t"<<"p-value adjust methods for H3K4me1" << endl;  
  //cerr << '\t' << "-3 [int]\t"<<"p-value adjust methods for H3K27ac" << endl;                              
  //cerr << '\t' << "-p [double]      pvalue[default: 0.01]" << endl;                                                                                                           
  cerr << '\t' << "-h                  \t"<<"print this help" << endl;
  cerr << " " << endl;                                                                                                                                                   
  return;
}
