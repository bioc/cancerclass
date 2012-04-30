\name{fit}
\alias{fit}

\title{Fitting of a predictor}
\description{
  Fits a predictor to a training data set.
}

\usage{
  fit(eset, class="class", method = "welch.test", hparam = 0.75)  
}

\arguments{
  \item{eset}{Bioconductor ExpressionSet}
  \item{class}{Character vector specifying classes related to the samples.} 
  \item{method}{Character string specifying the feature selection method. Possible values are "cor", "student.test", "welch.test", "wilcoxon.test", "foldchange", "copa", "os", "ort", "shift", "throw".}
  \item{hparam}{Hyperparamter needed for the feature selection methods: Confidential Interval for copa, ort, os (e.g. 0.75, 0.95). Minimum number of samples in each class after applying shift/throw (only necessary for the feature selection methods: throw, shift).}
}
\details{
The matrix \code{eset} contains the expression signatures of the patients in the columns. The vector \code{class} contains the class membership of each sample or patient.  Only two-class problems are supported. The colnames of \code{eset} are matched to the names of \code{classifier} (if both exist).

The hyperparameter \code{hparam} describes the minimum number of samples in each class after  applying shift/throw. For \code{copa} the hyperparameter is quanilte for the definition of outliers. Typical values are 0.75 (default), 0.9, 0.95.

A nearest centroid predictor is constructed by calculating the average level ofeach feature in each of the two classes of the trainig data set.
}

\value{
  A \code{predictor} object, see \code{predictor.object} for details.
}

\examples{
### see: help(GOLUB);
}

\seealso{ 
  \code{\linkS4class{predictor}}
}

\keyword{ classif }
