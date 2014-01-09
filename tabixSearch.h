#include <stdio.h>
#include <vector>
#include "drmPair.h"
#include "tabix.h"


class tabSearch
{
 public:
  tabix_t* tab;
  const char* gz;
  const char* index;
  
  tabSearch(const char* gzfile, const char* gzidx,int readcnt, int loadnow);
  ~tabSearch(){};
  int idxload();
  ti_iter_t search(const char* chrom,int start, int end);
  void getlineIter();
  double coverage(sk_conf_t type,ti_iter_t iter);
  void unloadidx();
  void coverage_adjustp(double pval, sk_conf_t datatype,vector<double> &adjp, ti_iter_t iter);
  

 private:
  sk_conf_t ticonf;
  //ti_iter_t iter;//ti_iter_destroy(iter)
  char qChrome[256];//q for query
  int qStart;
  int qEnd;
  int totalreads;
  double qRawp;
  
  
  //const char* line;
  
};

