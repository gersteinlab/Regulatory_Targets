/*
stat.h statistics function
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>
-lstdc++ -lgsl -lgslcblas
 */
#include "stdafx.h"
#include <vector>
#include <string>
#include <stdlib.h>
#include <assert.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_sort_vector.h>
//#include "spearman/gsl_stats_spearman.h"
#include "ap.h"
#include "statistics.h"
#define MAX(a, b) ((a) > (b) ? (a):(b) )

using namespace std;

/*
Resolve a sequence of ties.

The input ranks array is expected to take the same value for all indices in
tiesTrace. The common value is recoded with the average of the indices. For
example, if ranks = <5,8,2,6,2,7,1,2> and tiesTrace = <2,4,7>, the result
will be <5,8,3,6,3,7,1,3>.

Source: http://commons.apache.org/math/apidocs/src-html/org/apache/commons/math/stat/ranking/NaturalRanking.html#line.312
*/
void
resolveTies (double * ranks, size_t * tiesTrace, const size_t n_ties)
{
  
  size_t i;
  
  
  // constant value of ranks over tiesTrace
  double c = ranks[tiesTrace[0]];
  
  
  // new rank (ie. the average of the current indices)
  double avg = (2*c + n_ties - 1) / 2;
  
  
  for(i=0; i<n_ties; ++i)
    ranks[tiesTrace[i]] = avg;
  
}


/*
Rank data using the natural ordering on doubles, ties being resolved by taking their average.

Source: http://commons.apache.org/math/apidocs/src-html/org/apache/commons/math/stat/ranking/NaturalRanking.html#line.190
*/
void
gsl_stats_spearman_rank (double * ranks, double * d, size_t * p, const double data[], const size_t stride, const size_t n)
{
  
  size_t i;
  
  
  // copy the input data and sort them
  for(i=0; i<n; ++i)
    d[i] = data[i];
  gsl_sort (d, 1, n);
  // get the index of the input data as if they were sorted
  gsl_sort_index (p, data, stride, n);  
  // walk the sorted array, filling output array using sorted positions, resolving ties as we go
  size_t pos = 1;
  ranks[p[0]] = pos;
  size_t n_ties = 1;
  size_t * tiesTrace = (size_t*) calloc (1, sizeof(size_t));
  tiesTrace[0] = p[0];
  for(i=1; i<n; ++i)
  {
    if(d[i] - d[i-1] > 0)
    {   
      pos = i + 1; 
      if(n_ties > 1)
        resolveTies(ranks, tiesTrace, n_ties);
      tiesTrace = (size_t*) realloc (tiesTrace, sizeof(size_t));
      n_ties = 1;
      tiesTrace[0] = p[i];
    }
    
    else
    {
      ++n_ties; 
      tiesTrace = (size_t*) realloc (tiesTrace, n_ties * sizeof(size_t));
      tiesTrace[n_ties-1] = p[i];
    }
    ranks[p[i]] = pos;
  }
  
  if(n_ties > 1)
    resolveTies(ranks, tiesTrace, n_ties);
  free (tiesTrace);
  
}
void
gsl_stats_spearman_alloc (double ** ranks1, double ** ranks2, double ** d, size_t ** p, size_t n)
{
  
  *ranks1 = (double*) calloc (n, sizeof(size_t));
  *ranks2 = (double*) calloc (n, sizeof(size_t));
  *d = (double *) calloc (n, sizeof(double));
  *p = (size_t*) calloc (n, sizeof(size_t));
}
void
gsl_stats_spearman_free (double ** ranks1, double ** ranks2, double ** d, size_t ** p)
{  
  free (*ranks1);
  free (*ranks2);
  free (*d);  
  free (*p);
}
/*
Calculate Spearman correlation by computing Pearson correlation on the ranked variables.
Could become gsl_stats_correlation_spearman()
*/
double
gsl_stats_spearman (const double data1[], const size_t stride1,
                    const double data2[], const size_t stride2,
                    const size_t n)
{
  double rs = 0.0;
  double * ranks1, * ranks2, * d;
  size_t * p;
  gsl_stats_spearman_alloc (&ranks1, &ranks2, &d, &p, n);
  gsl_stats_spearman_rank (ranks1, d, p, data1, stride1, n);
  gsl_stats_spearman_rank (ranks2, d, p, data2, stride2, n);
  rs = gsl_stats_correlation((double*) ranks1, 1, (double*) ranks2, 1, n );
  gsl_stats_spearman_free (&ranks1, &ranks2, &d, &p);  
  return rs; 
}
double gsl_pearson_corr(vector<double> &x, vector<double> &y, double &bothtails,double &lefttail, double &righttail)
{
  const size_t stride = 1;
  int n=x.size();
  assert(x.size()==y.size());
  /*
  string xstr("c(");
  string ystr("c(");  
  for (int i=0; i< n-1; i++)
  {
    xstr.append(to_string(x[i])).append(",");
    ystr.append(to_string(y[i])).append(",");
  }
  xstr.append(to_string(x[n-1])).append(")");
  ystr.append(to_string(y[n-1])).append(")");
  
  string cmd("echo 'cor(");
  cmd.append(xstr).append(",").append(ystr).append(")' | R --no-save -q");
  system(cmd.c_str());
  string cmd2("echo 'cor.test(");
  cmd2.append(xstr).append(",").append(ystr).append(")' | R --no-save -q");
  system(cmd2.c_str());
  */  
  gsl_vector_const_view gsl_x = gsl_vector_const_view_array( &x[0], x.size() );
  gsl_vector_const_view gsl_y = gsl_vector_const_view_array( &y[0], y.size() );
  double pearson = gsl_stats_correlation( (double*) gsl_x.vector.data, stride,
                                          (double*) gsl_y.vector.data, stride,
                                          n );
  //determine whether -nan
  //stat<-sqrt(n-2)* pearson / sqrt(1- pearson * pearson);
  //printf("gsl:%f",pearson);
  alglib::pearsoncorrelationsignificance(pearson, (alglib::ae_int_t)n, bothtails, lefttail, righttail);
  bothtails=MAX(0,bothtails);
  lefttail=MAX(0,lefttail);
  righttail=MAX(0,righttail);
  return pearson; 
}

double spearman(vector<double> &x, vector<double> &y,double &bothtails, double &lefttail, double& righttail)
{
  const size_t stride = 1;
  int n=x.size();
  assert(x.size()==y.size());
  
  gsl_vector_const_view gsl_x = gsl_vector_const_view_array( &x[0], x.size() );
  gsl_vector_const_view gsl_y = gsl_vector_const_view_array( &y[0], y.size() );
  double spearman = gsl_stats_spearman( (double*) gsl_x.vector.data, stride,
                                        (double*) gsl_y.vector.data, stride, n);
  alglib::spearmanrankcorrelationsignificance(
                                              spearman,
                                              (alglib::ae_int_t)n ,
                                              bothtails,
                                              lefttail,
                                              righttail);  
  bothtails=MAX(0,bothtails);
  lefttail=MAX(0,lefttail);
  righttail=MAX(0,righttail);
  
  return spearman;
}

