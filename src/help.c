
#include <stdio.h>
#include <math.h>
#include "help.h"
#include "Rmath.h"

#define SORTMAX 500


double values[SORTMAX]={0.0};
double valuestmp[SORTMAX]={0.0};
char str[100];
/* bshift/bthrow */
double values1[SORTMAX] = { 0.0 };
double values2[SORTMAX] = { 0.0 };


// min() for ints
int mini(double a, double b)
{
  if (a<b) return(a);
  else return(b);  
}

double testx(double tmp)
{
    double erg=0;
    erg=tmp*4;

    return(erg);
}        
  
int vergleich(int *i1, int *i2)
{ return  *i1-*i2; }

// compare_double
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
int compare_doubles (const void *xx, const void *yy)
{
  const double *y = *((double **)xx);
  const double *x = *((double **)yy);

   if (*x > *y)
      return 1;
   else if (*x < *y)
      return -1;
   else
      return 0;
}


// getmax
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
double getmax(double *array, int from, int to)
{
   double tmp = 0.0;
   int i=0;

   tmp = array[from];
   for(i=from; i<(to-1); i++)
   {
      if (tmp < array[i+1]) tmp = array[i+1];
   }
   return tmp;
}


// mean
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
double mean(double *values, int from, int to)
{
    double result = 0.0;
    int nn = to-from;

    result = sum(values, from, to);
    result = result/nn;

    return( result );
}


// median
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
double median(double *values, int from, int to)
{
    double result = 0.0;
    int nsamples = to-from;
    int i,j = 0;

    /* sort */
    for (i=from; i<to; i++)
    {
	valuestmp[j] = values[i];
	j += 1;
    }

    sorttmp(nsamples);

    if (nsamples%2==1) { 
             result = valuestmp[((nsamples+1)/2)-1];
    } else {    
             result = 0.5 * (valuestmp[(nsamples/2)-1] + valuestmp[ nsamples/2 ] );
    }

    return( result );

}


// MAD median average deviation
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
double madian(double *values, int from, int to)
{
    double result = 0.0;
    double med = 0.0;
    int i;

    med = median(values, from, to);
    for (i=from; i<to; i++) {
        valuestmp[i] = fabs( values[i]-med );
    }

    sorttmp(to);
    med = median(valuestmp, from, to);

    result = 1.4826 * med;

    return( result );
}


// quantile type 1
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
double quantile(double *values, int from, int to, double r)
{
    double result = 0.0;
    double tmp = 0.0;
    int i,j = 0;
    double n=(double)(to-from);

    /* sort */
    for (i=from; i<to; i++) {
        valuestmp[j] = values[i];
        j += 1;
    }

    sorttmp(to-from);

    /* get quantil */
    for (i=0; i<(to-from); i++) {   

        tmp = ( (i+1) ) / n;
        if (tmp>=r) {
            result = valuestmp[i]; 
            break;
        }
    }

    return( result );
}

// quantile type 7 .R
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
double quantile7(double *values, int from, int to, double r)
{
    double result = 0.0;   
    double ktmp,rr;
    int k,i,j = 0;
    double n=(double)(to-from);
   
    /* sort */
    for (i=from; i<to; i++) {
        valuestmp[j] = values[i];
        j += 1;
    }  

    // sort valuestmp values
    sorttmp(to-from); 

    // continuous sample quantile
    ktmp = ((n-1)*r)+1;
    k = floor(ktmp);   
    rr = ktmp - k;

    result = valuestmp[k-1] + rr*( valuestmp[k] - valuestmp[k-1] );
 
    return( result );
}

// sum
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
double sum(double *values, int from, int to)
{
    int i=0;
    double result = 0.0;
    for (i=from; i<to; i++) {
        result += values[i];
    }

    return( result );
}


// IQR inter quantile range
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
double IQR(double *values, int from, int to)
{
    double result;
    double tmp1 = 0.0;
    double tmp2 = 0.0;

    tmp1 = quantile7(values, from, to, 0.75);      
    tmp2 = quantile7(values, from, to, 0.25);
    result = tmp1 - tmp2;

    return( result );

}


// student_t
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
double student_t(double *values1, double *values2, int ng1, int ng2)
{
    double mean1, mean2, var1, var2 = 0.0;
    double t=0.0;

    mean1=mean(values1,0,ng1);
    mean2=mean(values2,0,ng2);

    var1 = var(values1, 0, ng1);
    var2 = var(values2, 0, ng2);

    t = (mean1-mean2) / sqrt( (var2+var1) / ((ng1+ng2)-2) ); 
    t = t * sqrt((double)(ng1*ng2)/(double)(ng1+ng2)); 

    return(t);
}


// welch_t
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
double welch_t(double *values1, double *values2, int ng1, int ng2)
{
    double mean1, mean2, var1, var2 = 0.0;
    double t=0.0;

    mean1=mean(values1,0,ng1);
    mean2=mean(values2,0,ng2);

    var1 = var(values1, 0, ng1);
    var2 = var(values2, 0, ng2);

    t = (mean1-mean2) / sqrt( ((var1/(ng1-1))/ng1) + ((var2/(ng2-1))/ng2) );      
    return(t);
}


// welch_df
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
double welch_df(double *values1, double *values2, int ng1, int ng2)
{
    double mean1, mean2, var1, var2 = 0.0;
    double wdf=0.0;

    mean1=mean(values1,0,ng1);
    mean2=mean(values2,0,ng2);

    var1 = var(values1, 0, ng1);
    var2 = var(values2, 0, ng2);
    
    var1=var1/(ng1-1);
    var2=var2/(ng2-1);

    wdf = ( pow(var1/ng1 + var2/ng2, 2) ) / ( pow(var1/ng1,2)/(ng1-1) + pow(var2/ng2,2)/(ng2-1) );

    return wdf;
}


// pvalue welch/student
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
double p_value(double stat, double df)
{
    double pv=0.0;
    stat = fabs(stat);  
    pv = 2 * pt(stat, df, 0, 0); 
    return pv;
}


// pvalue wilcoxon
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
double p_wilcoxon(double stat, double m, double n)
{
    double pv=0.0;
    stat = fabs(stat); 
    pv = 2 * pwilcox(stat, m, n, 1, 0); 
    return pv;
}


// sort
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*void sort(double *values, int nsamples) 
{ 
    R_qsort(values, 1, nsamples); 
}
*/

// sort tmp
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void sorttmp(int nsamples) 
{
    R_qsort(valuestmp, 1, nsamples);    
}


// var (variance)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
double var(double *array, int from, int to) //from,to stimmt
{
   double var=0.0;
   double m=0.0;
   int count=0;

   m=mean(array, from, to); //stimmt

   for(int j = from; j < to; j++)
   {
       var = var + ( (array[j]-m) * (array[j]-m) );
       count+=1;
   }

   var=var; ///(to-from-1); R=var(N-1) => /(count-1)
   return(var);
}


// sort bsortdesc: descending order
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void bsortdesc(double *array, int from, int to)
{
    int i,k,anzahl=0;
    double temp=0.0;
    anzahl=from-to;

    for (k=0; k<(to-1); k++)                         // Äußere Schleife
    {
       for (i=0; i<(to-1); i++)                      // innere Schleife
       {
	  if(array[i]<array[i+1])                     // Abfrage, ob der aktuelle Wert größer als der nächste ist
	  {
	         temp = array[i];                        // dann wird der aktuelle in die temp-Variable geschrieben
	         array[i] = array[i+1];                  // die Werte werden getauscht
	         array[i+1] = temp;                      // der Wert der temp-Variable wir an den nächsten zugewießen
	  }

       }
    }
}


// sort bsort: ascending order
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void bsort(double *array, int from, int to)
{
    int i,k,anzahl=0;
    double temp=0.0;
    anzahl=from-to;

    for (k=0; k<(to-1); k++)                         // Äußere Schleife
    {
       for (i=0; i<(to-1); i++)                      // innere Schleife
       {
	  if(array[i]>array[i+1])                     // Abfrage, ob der aktuelle Wert größer als der nächste ist
	  {
	         temp = array[i];                        // dann wird der aktuelle in die temp-Variable geschrieben
	         array[i] = array[i+1];                  // die Werte werden getauscht
	         array[i+1] = temp;                      // der Wert der temp-Variable wir an den nächsten zugewießen
	  }

       }
    }
}

