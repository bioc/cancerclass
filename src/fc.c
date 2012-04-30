
/* R */
#include "R.h" 
#include "Rinternals.h" 
#include "Rdefines.h"

/* C */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "help.h"

/* variables */
#ifndef SORTMAX
#define SORTMAX 500
extern double values[SORTMAX];
extern double valuestmp[SORTMAX];
extern char str[100];
#endif


/* implementations */
int fc( double *eset, double *ngenes, double *ng1, double *ng2, double *result ) 
{ 
   int nsamples = *ng1+*ng2;
   int currentgene =0;
   int i=0;
   int j=0;
   double m1,m2 = 0.0;

   currentgene = 0;
   for(i = 0; i < (*ngenes); i++) 
   {
         for(j = 0; j < nsamples; j++) 
         {

             values[j] = eset[currentgene];
             currentgene+=1;
         }
         // means
         m1 = mean(values,0, *ng1);
         m2 = mean(values,*ng1, nsamples);
         result[i] = m1-m2;
    }
    return(1);
}

