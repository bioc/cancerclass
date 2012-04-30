
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
double calcmad(double *values, int ng1, int ng2, double med1, double med2);
/* implementations */

int ort( double *eset, double *ngenes, double *ng1, double *ng2, double *r, double *result ) 
{
   int nsamples = *ng1+*ng2;
   int currentgene, i, j;
   double tmp = 0.0;
   double med1 = 0.0;
   double med2 = 0.0;
   double mad = 0.0;
   double q1 = 0.0;
   double q2 = 0.0;
   double q3 = 0.0;
   double q4 = 0.0;
   double iqrX = 0.0;
   double iqrY = 0.0;

   // all sample median , mad und q(r)-quantil
   currentgene = 0;
   for(i = 0; i < (*ngenes); i++)
   {
         for(j = 0; j < nsamples; j++)
         {

             values[j] = eset[currentgene];
             currentgene+=1;
         }

         med1 = median(values, 0, *ng1 );
         med2 = median(values, (*ng1), nsamples );
         mad = calcmad(values, *ng1, *ng2, med1, med2);

         q1 = quantile7(values, 0, *ng1, *r ); //g1
         q2 = quantile7(values, 0, *ng1, (1-*r) ); //g1
         q3 = quantile7(values, *ng1, nsamples, *r ); //g2
         q4 = quantile7(values, *ng1, nsamples, (1-*r) ); //g2

         iqrX = IQR(values, 0, *ng1);
         iqrY = IQR(values, *ng1, nsamples);

         q1 = q1+iqrX;  //g1+iqr(g1)
         q2 = q2-iqrX;  //g1-iqr(g1)
         q3 = q3+iqrY;  //g2+iqr(g2)
         q4 = q4-iqrY;  //g2-iqr(g2)

         // get Subset of > .75quantile7
         valuestmp[0]=0.0;
         valuestmp[1]=0.0;
         valuestmp[2]=0.0;
         valuestmp[3]=0.0;
	for(j = *ng1; j < nsamples; j++)
	{
	        if ( values[j] > q1 ) valuestmp[0] += values[j]-med1;
                 if ( values[j] < q2 ) valuestmp[1] += values[j]-med1;
         }
         for(j = 0; j < *ng1; j++)       
	{
	        if ( values[j] > q3 ) valuestmp[2] += values[j]-med2;
                 if ( values[j] < q4 ) valuestmp[3] += values[j]-med2;
         }

         valuestmp[0] = fabs(valuestmp[0]);
	valuestmp[1] = fabs(valuestmp[1]);
         valuestmp[2] = fabs(valuestmp[2]);
         valuestmp[3] = fabs(valuestmp[3]);

         tmp = getmax(valuestmp, 0, 4);
         result[i] = tmp / mad;
    }
    return(1);
}

double calcmad(double *values, int ng1, int ng2, double med1, double med2)
{
   double res = 0.0;
   int nn = ng1+ng2;
   for(int i=0; i<nn; i++)
   {
      if (i < ng1) valuestmp[i] = fabs(values[i]-med1);
      else         valuestmp[i] = fabs(values[i]-med2);
   }

   res = median(valuestmp, 0, nn);
   return(res);
}
