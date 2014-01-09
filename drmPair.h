#ifndef SK_DRMPAIR
#define SK_DRMPAIR
#include <iostream>
#include <string>

using namespace std;

typedef struct _drmtarget
{ 
  string chrom;
  int start;
  int end;
  string drmId;
  string geneId;
  string geneName;
  string strand;
  int dist;
  int closestDist;
} DRMTarPair;


typedef struct 
{
  int32_t present;//0 meth, 1, h3k4me1, 2,h3k27ac
  int32_t sc, bc, ec;// seq col., beg col. and end col.
  int32_t v1c, v2c; //values
}  sk_conf_t;

#endif
