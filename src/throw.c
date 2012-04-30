
/* R */
#include "R.h" // nicht gut
#include "Rinternals.h" // nicht gut
#include "Rdefines.h"

/* C */
#include <string.h>
//#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include "help.h"

/* variables help.h */
#ifndef SORTMAX
#define SORTMAX 500
extern double values[SORTMAX];
extern double valuestmp[SORTMAX];
extern char str[100];
extern double values1[SORTMAX];
extern double values2[SORTMAX];
#endif


/* implementations */
//.c
int throw( double *eset, double *ngenes, double *ng1, double *ng2, double *min, double *result) 
{
    int nsamples = *ng1+*ng2;
    int tmpng1 = *ng1;
    int tmpng2 = *ng2;
    int i, j, currentgene, throws;
    double t = 0.0; // t-statistic value
    double df, pvalue;

    //how often to throw away
    throws=mini(*ng1,*ng2) - *min + 1;
    if (throws<1) throws=1;

    //char stmp[2000];

    // calc
    currentgene = 0; //gene variable
    for(i = 0; i < (*ngenes); i++)
    {
    	// values of gene X cover in two arrays, each for a group
	for(j = 0; j < *ng1; j++)
         {
             values1[j] = eset[currentgene];
             currentgene+=1;
         }
	for(j = 0; j < *ng2; j++)
         {
             values2[j] = eset[currentgene];
             currentgene+=1;
         }

         // sort vectors
         bsortdesc(values2, 0, *ng2); //if outlier in group 2
	bsortdesc(values1, 0, *ng1); //if outlier in group 1

	// calc every throw proc. group 2
         for(j = 0; j < throws; j++)
         {
              t = fabs(student_t(values1, values2, *ng1, *ng2));
              df = nsamples-j-2;
 	     pvalue = p_value(t,df);

              if (result[i] > pvalue || result[i]==0) result[i] = pvalue;

              //throw away biggest value of outlier group
              *ng2-=1;
         }

         // init/set original set size
         *ng1=tmpng1;
         *ng2=tmpng2;

         // calc every throw proc. group 1
         for(j = 0; j < throws; j++)
         {
              t = fabs(student_t(values1, values2, *ng1, *ng2));
              df = nsamples-j-2;
 	     pvalue = p_value(t,df);

	     if (result[i] > pvalue || result[i]==0) result[i] = pvalue;

              //throw away biggest value of outlier group
              *ng1-=1;
         }

         // init/set original set size
         *ng1=tmpng1;
         *ng2=tmpng2;

         //check if values only 0.0
         if (sum(values1, 0, *ng1)==0) result[i]=99;
	if (sum(values2, 0, *ng2)==0) result[i]=99;
    }
    return(1);
}
