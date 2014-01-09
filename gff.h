#include <string>
#include "region.h"
#ifndef SK_GFF_H
#define SK_GFF_H
using namespace std;

class gff:public Region
{
 public:
  
  gff(string fn);
  void read();
  virtual ~gff(){};
  void add(const char* chr,int start, const char* strand, double readcnt){};
  
 private:
  string fname;

};


#endif
