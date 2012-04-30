

#include <stdio.h>
#include <R_ext/Rdynload.h>
#include "help.h"
#include "copa.h"
#include "ort.h"
#include "os.h"

#include "cor.h"
#include "fc.h"
#include "student.h"
#include "welch.h"
#include "wilcox.h"

#include "shift.h"
#include "throw.h"

static const
R_CMethodDef cMethods[] = {
	{"copa", (DL_FUNC) &copa, 6},
	{"ort", (DL_FUNC) &ort, 6},
	{"os", (DL_FUNC) &os, 6},
	{"cor", (DL_FUNC) &cor, 5},
	{"fc", (DL_FUNC) &fc, 5},
	{"student", (DL_FUNC) &student, 5},
	{"welch", (DL_FUNC) &welch, 5},
	{"wilcoxon", (DL_FUNC) &wilcoxon, 5},
	{"shift", (DL_FUNC) &shift, 6},
	{"throw", (DL_FUNC) &throw, 6},
	{NULL, NULL, 0}
};

void R_init_statistics(DllInfo *info)
{
	R_registerRoutines(info,
		cMethods, NULL, NULL, NULL);
}

int statistics( int *fx, double *eset, double *ngenes, double *ng1, double *ng2, double *r, double *result )
{
  int q=0;
  int i=0;
  int n=*ng1+*ng2;
  double tmp_result[(int)*ngenes];

  switch (*fx)
  {
    case 1: // CORRELATION
    	q=cor(eset, ngenes, ng1, ng2, result);
    	break;
    case 2: // WELCH.test
    	q=pwelch(eset, ngenes, ng1, ng2, result);
        //for(i=0;i<*ngenes;i++) result[i]=1-result[i];
    	break;
    case 3: // WELCH
    	q=welch(eset, ngenes, ng1, ng2, result);
        //for(i=0;i<*ngenes;i++) result[i]=1-result[i];
    	break;    	
    case 4: // STUDENT.test
	    q=pstudent(eset, ngenes, ng1, ng2, result);
    	break;
    case 5: // STUDENT
	    q=student(eset, ngenes, ng1, ng2, result);
    	break;
    case 6: // WILCOXON.test = p-values
	    q=pwilcoxon(eset, ngenes, ng1, ng2, result);
    	break;
    case 7: // WILCOXON
	    q=wilcoxon(eset, ngenes, ng1, ng2, result);
    	break;    	
    case 8: // FOLDCHANGE (foldchange)
    	q=fc(eset, ngenes, ng1, ng2, result);
    	break;
    case 9: // FOLDCHANGE (fc)
    	q=fc(eset, ngenes, ng1, ng2, result);
    	break;
    case 10: //COPA
    	q=copa(eset, ngenes, ng1, ng2, r, result);
    	break;    	
    case 11: // ORT
    	q=ort(eset, ngenes, ng1, ng2, r, result);
    	break;    	
    case 12: // OS
	    q=os(eset, ngenes, ng1, ng2, r, result);
    	break;

    case 13: // SHIFT
	    q=shift(eset, ngenes, ng1, ng2, r, result);
          // CHANGE SIGN OF NUMBER OF ESET    	
          for(i=0;i<(*ngenes*n);i++) eset[i]=eset[i]*(-1);
          // SAVE IN tmp_result CURRENT results
          for(i=0;i<(*ngenes);i++) tmp_result[i]=result[i];
          
      q=shift(eset, ngenes, ng1, ng2, r, result);
          for(i=0;i<(*ngenes);i++)
          { 
            if (result[i]>tmp_result[i]) result[i]=tmp_result[i];
            result[i]=1-result[i];
          }  
                
    	break;
    case 14: //THROW
    	q=throw(eset, ngenes, ng1, ng2, r, result);
          // CHANGE SIGN OF NUMBER OF ESET    	
          for(i=0;i<(*ngenes*n);i++) eset[i]=eset[i]*(-1);
          // SAVE IN tmp_result CURRENT results
          for(i=0;i<(*ngenes);i++) tmp_result[i]=result[i];
          
      q=throw(eset, ngenes, ng1, ng2, r, result);
        // CHECK WHICH PVALUE IS MOST LESS
          for(i=0;i<(*ngenes);i++)
          { 
            if (result[i]>tmp_result[i]) result[i]=tmp_result[i];
            result[i]=1-result[i];
          }  
                
    	break;

    default: q=-1;
  }

  return q;
}


