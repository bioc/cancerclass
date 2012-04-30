\name{loo}
\alias{loo}
\title{Leave-one-out cross-validation}
\description{
  Fitting and validation of a predictor in a leave-one-out protocol. 
}
\usage{
  loo(eset, class="class", method = "welch.test", ngenes=50, dist="cor", hparam = 0.75, positive="") 
}

\arguments{
  \item{eset}{Bioconductor ExpressionSet} 
  \item{class}{String specifying the column in pData(eset) that contains the class information.} 
  \item{method}{Specifying the feature selection method. Possible values are "cor", "student.test", "welch.test", "wilcoxon.test", "foldchange", "copa", "os", "ort", "shift", "throw".}
  \item{ngenes}{Number of features used for classification.}
  \item{dist}{Metric for distance calculation}
  \item{hparam}{Hyperparameter needed for some of the feature selection methods. For methods copa, ors and os: Quantile (e.g. 0.75, 0.9, 0.95) used in order to detect outliers. For methods shift and throw: the minimum number of samples in each class after applying shift or throw.}
  \item{positive}{One of the two classes. Membership to this class is considered as positive.    Needed in order to calculate sensitivity and specificity of the validation.}
}

\details{
A leave-one-out cross-validation is performend by calling fit and predict in a loop.
}

\value{
  A \code{pvalidation} object, see \code{pvalidation.object} for details.
}

\examples{
### see: help(GOLUB);
}

\keyword{ classif }
