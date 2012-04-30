
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
int welch( double *eset, double *ngenes, double *ng1, double *ng2, double *result) 
{
    int i, j, currentgene;
    double t = 0.0; // t-statistic value
  
    // calc
    currentgene = 0; //gene variable
    for(i = 0; i < (*ngenes); i++)
    {
    	// values of gene X cover in two arrays, each for a group; prepare gene
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

         t = welch_t(values1, values2, *ng1, *ng2);
         result[i] = t;
    }
    return(1);
}

int pwelch( double *eset, double *ngenes, double *ng1, double *ng2, double *result) 
{
    int i, j, currentgene;
    double t = 0.0; // t-statistic value
    double df, pvalue;

    // calc
    currentgene = 0; //gene variable
    for(i = 0; i < (*ngenes); i++)
    {
    	// values of gene X cover in two arrays, each for a group; prepare gene
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

         //t = fabs(welch_t(values1, values2, *ng1, *ng2));
         t = welch_t(values1, values2, *ng1, *ng2);
         df = welch_df(values1, values2, *ng1, *ng2);
 	       pvalue = p_value(t,df);
         result[i] = pvalue;
    }
    return(1);
}

