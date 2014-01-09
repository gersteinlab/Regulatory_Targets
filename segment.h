
#include <iostream>
#include <string>
#include <vector>
#include <map>
//can test multimap
using namespace std;

#ifndef SK_SEGMENT_H
#define SK_SEGMENT_H
class Segment
{
 public:
  string chrom;
  vector<int> chromStart;
  vector<int> chromEnd;
  vector<string> strand;
  vector<string> id;
  vector<double> value;
  vector<int> bin;
  vector<double> value2;
  string name;
};


#endif





  
