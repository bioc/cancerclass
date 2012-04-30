
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
double calcrang(double *values, int from, int to);
void calcwilcoxstat(double *values, int from, int to);

/* implementations */

int wilcoxon( double *eset, double *ngenes, double *ng1, double *ng2, double *result ) {

   #define n 40000
   int nsamples = *ng1+*ng2;
   int currentgene, i, j;
   //double tmp = 0.0;
   int min1 = 0;
   int min2 = 0;
   int overall = 0;
   double rang1 = 0.0;
   double rang2 = 0.0;

   // all sample median , mad und q(r)-quantil
   currentgene = 0;
   for(i = 0; i < (*ngenes); i++) // für jedes gen mache:
   {
         for(j = 0; j < nsamples; j++) // alle samples in array zusammenfassen
         {

             values[j] = eset[currentgene];
             currentgene+=1;
         }

         //ok

	       min1 = (*ng1*(*ng1+1)) / 2;
	       min2 = (*ng2*(*ng2+1)) / 2;
         overall = ((nsamples+1)*nsamples)/2;  
         calcwilcoxstat(values, 0, nsamples);               
         rang1 = calcrang(valuestmp, 0, *ng1) - min1;
         rang2 = calcrang(valuestmp, *ng1, nsamples) - min2;

         if (rang1 < rang2) rang1 = rang2;       
         result[i] = rang1;
    }

    return(1);
}

int pwilcoxon( double *eset, double *ngenes, double *ng1, double *ng2, double *result ) {

   #define n 40000
   int nsamples = *ng1+*ng2;
   int currentgene, i, j;
   int min1 = 0;
   int min2 = 0;   
   int overall = 0;
   double rang1 = 0.0;
   double rang2 = 0.0;

   currentgene = 0;
   for(i = 0; i < (*ngenes); i++) 
   {
         for(j = 0; j < nsamples; j++) 
         {

             values[j] = eset[currentgene];
             currentgene+=1;
         }

         //ok

	       min1 = (*ng1*(*ng1+1)) / 2;
         min2 = (*ng2*(*ng2+1)) / 2;	       
         overall = ((nsamples+1)*nsamples)/2;
         calcwilcoxstat(values, 0, nsamples); 
         rang1 = calcrang(valuestmp, 0, *ng1) - min1; 
	       rang2 = calcrang(valuestmp, *ng1, nsamples) - min2;            

         if (rang1 < rang2) rang1 = p_wilcoxon(rang1,*ng1,*ng2);
         else rang1 = p_wilcoxon(rang2,*ng2,*ng1);                
         result[i] = rang1;        

    }

    return(1);
}

/* -------------------------------------------------------------------------------------------------- */
void calcwilcoxstat(double *array, int from, int to)
{
   int i=0;
   int j=0;
   int rang=1;
   for(i=from; i<to; i++)
   {
      rang=1;
      for(j=from; j<to; j++)
      {
      	if (array[i] > array[j]) rang +=1;
      }

      valuestmp[i] = rang;
   }
}

double calcrang(double *array, int from, int to)
{
   double tmp = 0.0;
   int i=0;
   for(i=from; i<to; i++)
   {
	   tmp += array[i];
   }
   return tmp;
}
