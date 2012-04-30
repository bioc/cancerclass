
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

double ideal[SORTMAX] = { 0.0 };
 

/* functions */
double skalar(double *array, double *array2, int nsamples);
double norm(double *, int);
void center(double mittel,int nsamples);
void changeideal(double *array, int nsamples);

/* implementations */

int cor( double *eset, double *ngenes, double *ng1, double *ng2, double *result ) {

   int nsamples = *ng1+*ng2;
   int currentgene, i, j = 0;
   double mittel  = 0.0;
   double norm1 = 0.0;
   double norm2 = 0.0;
   double tmp1 = 0.0;
   double tmp2 = 0.0;

   for(i = 0; i < *ng1; i++)
   {
       ideal[i] = -1.0;
   }
   for(i = *ng1; i < nsamples; i++)
   {
       ideal[i] = 1.0;
   }


   currentgene = 0;
   for(i = 0; i < (*ngenes); i++)
   {
         for(j = 0; j < nsamples; j++)
         {

             values[j] = eset[currentgene];
             currentgene+=1;
         }

         mittel = mean(values, 0, nsamples);
         center(mittel, nsamples);

	norm1 = norm(ideal, nsamples);
	norm2 = norm(valuestmp, nsamples);

         tmp1 = skalar(ideal, valuestmp, nsamples) / (norm1 * norm2);
         changeideal(ideal,nsamples);
         tmp2 = skalar(ideal, valuestmp, nsamples) / (norm1 * norm2);

         if (tmp1 < tmp2) tmp1=tmp2;

         result[i]  = tmp1;
    }

    return(1);
}


/* -------------------------------------------------------------------------------------------------- */
void changeideal(double *array, int nsamples)
{
    int i;
    for(i = 0; i < nsamples; i++)
    {
	array[i] = array[i] * (-1);
    }
}

double norm(double *array, int nsamples)
{
    double result = 0.0;
    int i;
    for(i = 0; i < nsamples; i++)
    {
	result = result + (array[i]*array[i]);
    }
    return sqrt(result);
}

double skalar(double *array, double *array2, int nsamples)
{
    double result = 0.0;
    int i;
    for(i = 0; i < nsamples; i++)
    {
	result = result + (array[i]*array2[i]);
    }
    return (result);
}

void center(double mittel,int nsamples)
{
    int i;
    for(i = 0; i < nsamples; i++)
    {
	valuestmp[i] = values[i]-mittel;
    }
}
