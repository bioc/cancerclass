\name{validation-class}
\docType{class}
\alias{validation-class}

\title{ Class "validation"   }
\description{	 
  An object of this class is returned by the function \code{validate} and represents the 
  validation of a nearest centroid predictor in a random validation protocol. 
}
\section{Objects from the Class}{
  Objects can of this class be created by \code{new("validation", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{ngenes}:}{Numerical vector containing the numbers features that are used for classification.}                                           
    \item{\code{method}:}{Character string specifying the feature selection method. Possible values are "cor", "student.test", "welch.test", "wilcoxon.test", "foldchange", "copa", "os", "ort", "shift", "throw".}
    \item{\code{dist}:}{Character string specifiying the method for calculation of the distance between test samples and the centroids. Possible values are "euclidean", "angle", "cor", "center". } 
    \item{\code{ntrain}:}{A numeric vector containing the numbers of samples in the training sets. }
    \item{\code{nrep}:}{The number of repetitions for each training set size. }
    \item{\code{hparam}:}{Hyperparameter needed for some of the feature selection methods. For methods copa, ors and os: Quantile (e.g. 0.75, 0.9, 0.95) used in order to detect outliers. For methods shift and throw: the minimum number of samples in each class after applying shift or throw. }  
    \item{\code{misclass}:}{A list containing  the total, the class1 and the class2 misclassification rates. }
    \item{\code{nselected}:}{Contains information, how often each of the genes was selected for a predictor. }
    \item{\code{samples}:}{Numeric matrix containing the classification rates for each of the samples. }
    \item{\code{classifier}:}{Numerical or character vector containing the class membership of the samples. }
    \item{\code{fdata}:}{Numerical or character vector containing the class membership of the samples. }
  }
}
\section{Methods}{
  There are three \code{plot} methods for visualization of the validation results.
}

\seealso{
  \code{\link{validate}}, 
  \code{\linkS4class{validation}} 
}
\examples{
showClass("validation")
}

\keyword{classes}
