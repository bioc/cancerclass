
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
double student_t(double *values1, double *values2, int ng1, int ng2);
void shift2to1(double *values1, double *values2, int ng1, int ng2);
void shift1to2(double *values1, double *values2, int ng1, int ng2);

//.c
int shift( double *eset, double *ngenes, double *ng1, double *ng2, double *min, double *result ) 
{
    int nsamples = *ng1+*ng2;
    int tmpng1 = *ng1;
    int tmpng2 = *ng2;
    int i, j, currentgene, shifts;
    double t = 0.0; // t-statistic value
    double df;
    double pvalue;

    //how often to throw away
    shifts=mini(*ng1,*ng2) - *min + 1;
    if (shifts<1) shifts=1;	
	
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
        
         bsortdesc(values2, 0, *ng2); //if outlier in group 2
         bsortdesc(values1, 0, *ng1); //if outlier in group 1


         for(j = 0; j < shifts; j++)
         {
	     // calc t-statistic
 	     t = fabs(student_t(values1, values2, *ng1, *ng2));
              df = nsamples-2;
 	     pvalue = p_value(t,df);
              if (result[i] > pvalue || result[i]==0) result[i] = pvalue;

              shift2to1(values1, values2, *ng1, *ng2);
	     *ng2-=1;
	     *ng1+=1;
         }

         // init/set original set size
         *ng1=tmpng1;
         *ng2=tmpng2;

         for(j = 0; j < shifts; j++)
         {
              // calc t-statistic
 	     t = fabs(student_t(values1, values2, *ng1, *ng2));
              df = nsamples-2;
 	     pvalue = p_value(t,df);
              if (result[i] > pvalue || result[i]==0) result[i] = pvalue;

              shift1to2(values1, values2, *ng1, *ng2);
    	     *ng2+=1;
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

void shift2to1(double *values1, double *values2, int ng1, int ng2)
{
    values1[ng1] = values2[ng2-1];
}
void shift1to2(double *values1, double *values2, int ng1, int ng2)
{
    values2[ng2] = values1[ng1-1];
}
