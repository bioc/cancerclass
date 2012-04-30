\name{nvalidation-class}
\docType{class}
\alias{nvalidation-class}

\title{Class "nvalidation"   }
\description{	
  This class of objects is returned by the function \code{nvalidation} 
  and represents the validation 
  of a nearest centroid predictor in a random validation protocol.
  Objects of this class have methods for the function \code{plot}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("nvalidation", ...)}.
 describe objects here 
}
\section{Slots}{
	 \describe{
    \item{\code{ngenes}:}{Numerical vector containing the numbers features that are used for classification.}    
    \item{\code{method}:}{Character string specifying the feature selection method. Possible values are "cor", "student.test", "welch.test", "wilcoxon.test", "foldchange", "copa", "os", "ort", "shift", "throw".}
    \item{\code{dist}:}{Character string specifiying the method for calculation of the distance between test samples and the centroids. Possible values are "euclidean", "angle", "cor", "center".}
    \item{\code{ntrain}:}{The number of samples in the training set. }
    \item{\code{nrep}:}{The number of repetitions for each training set size.}
    \item{\code{hparam}:}{Hyperparameter needed for some of the feature selection methods. For methods copa, ors and os: Quantile (e.g. 0.75, 0.9, 0.95) used in order to detect outliers. For methods shift and throw: the minimum number of samples in each class after applying shift or throw. }  
    \item{\code{misclass}:}{A list containing the total, the class1 and the class2 misclassification rates. }
    \item{\code{nselected}:}{Contains information, how often each of the genes was selected for a predictor. }
    \item{\code{samples}:}{Numeric matrix containing the classification rates for each of the samples. }
    \item{\code{classifier}:}{Numerical or character vector containing the class membership of the samples. }
    \item{\code{fdata}:}{Feature annotations inherited from the ExpressionSet. }
  }
}
\section{Methods}{
  There are three \code{plot} methods for visualization of the validation results.
}
 
\seealso{
  \code{\link{nvalidate}}, 
  \code{\linkS4class{nvalidation}} 
}
\examples{
showClass("nvalidation")
}
\keyword{classes}
