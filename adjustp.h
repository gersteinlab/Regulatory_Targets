#include <iostream>
#include <string>
#include <vector>
#include <math.h>
#include "drmPair.h"

//DRMTarPair
//map<int, DRMTarPair>  
using namespace std;

class ADJP 
{
 public:
 ADJP(int idx0, double corr0, double rawp0,string type0):idx(idx0),corr(corr0),rawp(rawp0), type(type0){};
  
  bool operator<(const ADJP& adjp)
  {
    return this->rawp < adjp.rawp;  
  }
  string type;
  int idx;
  double corr, rawp;
};



//std::sort(vecAdjp.begin(),vecAdjp.end());
double sumInverseM(int n)
{
  double rtn=0.0;
  for (int i=1; i<=n; i++)
    rtn+=(1/i);
  return rtn;
}

//template<class T>
double** adjust(vector<ADJP> & sortedVec)
{
  int M=sortedVec.size();
  double A=sumInverseM(M);
  double** adjps=new double*[M];
  for (int i=0; i< M; i++)
  {
    adjps[i]=new double[4];
    for (int j=0; j< 4; j++)
    {
      adjps[i][j]=1.0;
    }    
  }
  
  adjps[0][0]=min(1.0, M*sortedVec[0].rawp);//Bonferroni
  adjps[0][1]=min(1.0, M*sortedVec[0].rawp); //Holm
  adjps[M-1][2]=sortedVec[M-1].rawp;//BH
  adjps[M-1][3]=min(A*sortedVec[M-1].rawp, 1.0); //BY

  for (int i=1; i< M; i++)
  {
    adjps[i][0]=min(1.0, M*sortedVec[i].rawp);
    adjps[i][1]=max(adjps[i-1][1], min(1.0, (M-i)*sortedVec[i].rawp));
  }
  
  for (int i= M-2; i>=0; i--)
  {
    adjps[i][2]=min(adjps[i+1][2], min((M/(i+1))*sortedVec[i].rawp, 1.0));
    adjps[i][3]=min(adjps[i+1][3], min( (M*A/(i+1))*sortedVec[i].rawp,1.0));
  }

  return adjps;
}


