
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
 

/* functions */

/* implementations */

int copa( double *eset, double *ngenes, double *ng1, double *ng2, double *r, double *result ) 
{
   int nsamples = *ng1+*ng2;	//ok
   int currentgene, i, j;
   double med,mad;
   double tmp = 0.0;   
   double q1 = 0.0;
   double q2 = 0.0;
   double q3 = 0.0;
   double q4 = 0.0;

   // all sample median , mad und q(r)-quantil
   currentgene = 0;
   for(i = 0; i < (*ngenes); i++)
   {
         for(j = 0; j < nsamples; j++)
         {

             values[j] = eset[currentgene];
             currentgene+=1;
         }
             
         q1 = quantile7( values, (*ng1), nsamples, *r );  //percentile from cancer samples; quantile(double r, int from, int to)         
         q2 = quantile7( values, (*ng1), nsamples, (1-*r) );  //percentile from cancer samples
         q3 = quantile7( values, 0, (*ng1), *r );
         q4 = quantile7( values, 0, (*ng1), (1-*r) );

         med = median(values, 0, nsamples ); //ERROR!!!         
         mad = madian(values, 0, nsamples );                  

         valuestmp[0] = fabs( (q1 - med)/mad );
         valuestmp[1] = fabs( (q2 - med)/mad );
         valuestmp[2] = fabs( (q3 - med)/mad );
         valuestmp[3] = fabs( (q4 - med)/mad );

         tmp = getmax(valuestmp, 0, 4);
         result[i] = tmp;
    }
    return(1);
}
