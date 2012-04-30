
/* R */
#include "R.h" // nicht gut
#include "Rinternals.h" // nicht gut
#include "Rdefines.h"

/* C */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
//#include <windows.h>
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

int os( double *eset, double *ngenes, double *ng1, double *ng2, double *r, double *result ) 
{
   int nsamples = *ng1+*ng2;
   int currentgene, i, j;
   double tmp = 0.0;
   double med, mad, iqr;
   double q75 = 0.0;
   double q25 = 0.0;

   // all sample median , mad und q(r)-quantil
   currentgene = 0;
   for(i = 0; i < (*ngenes); i++) 
   {
         for(j = 0; j < nsamples; j++) 
         {

             values[j] = eset[currentgene];
             currentgene+=1;
         }

         med = median(values, 0, nsamples );
         mad = madian(values, 0, nsamples );

         q75 = quantile7(values, 0, nsamples, *r ); 
         q25 = quantile7(values, 0, nsamples, (1-*r) );
         iqr = IQR(values, 0, nsamples);

         q75 = q75 + iqr;
         q25 = q25 - iqr;

         // add/get Subset of > .75quantile7
         valuestmp[0]=0.0;
         valuestmp[1]=0.0;
         valuestmp[2]=0.0;
         valuestmp[3]=0.0;
         for(j = 0; j < *ng1; j++)      
	{
	        if ( values[j] > q75 ) { valuestmp[0] += values[j]-med; }
                 if ( values[j] < q25 ) { valuestmp[1] += values[j]-med; }
         }

	for(j = *ng1; j < nsamples; j++)
	{
	        if ( values[j] > q75 ) { valuestmp[2] += values[j]-med; }
                 if ( values[j] < q25 ) { valuestmp[3] += values[j]-med; }
         }


        valuestmp[0] = fabs(valuestmp[0]); //mad;
		valuestmp[1] = fabs(valuestmp[1]); //mad;
        valuestmp[2] = fabs(valuestmp[2]); //mad;
        valuestmp[3] = fabs(valuestmp[3]); //mad;

        tmp = getmax(valuestmp, 0, 4);
        result[i] = tmp / mad;

    }
    return(1);
}
