%R CMD check -no-examples cancerclass
%R CMD build --binary cancerclass
\name{cancerclass-package}
\alias{cancerclass}
\docType{package}

\title{
Development and validation of diagnostic tests from high-dimensional molecular data
}

\description{
This package implements classification and validation methods for high-dimensional applications,  such as gene expression data. The classification protocol starts with a feature selection step and continues with nearest-centroid classification. The accurarcy of the predictor can be evaluated using training and test set validation, leave-one-out cross-validation or in a multiple random validation protocol [1]. Methods for calculation and visualization of continuous prediction scores allow to balance sensitivity and specificity and define a cutoff value according to clinical requirements.
}

\details{
\tabular{ll}{
Package: \tab cancerclass\cr
Type: \tab Package\cr
Version: \tab 1.5.1 \cr
Date: \tab 2013-09-04\cr
License: \tab GPL (>=2)
}}

\author{
Jan Budczies \email{jan.budczies@charite.de},
Daniel Kosztyla \email{danielkossi@hotmail.com}
}

\references{
[1] Michiels S, Koscielny S, Hill C (2005), \emph{Prediction of cancer outcome
with microarrays: a multiple random validation strategy}, Lancet 365:488-492.
}

\seealso{ 
  \code{\link{fit}},
  \code{\link{GOLUB1}},
  \code{\link{loo}},
  \code{\link{nvalidate}},
  \code{\link{nvalidation-class}},
  \code{\link{plot}},
  \code{\link{plot,nvalidation-method}},
  \code{\link{plot,prediction-method}},
  \code{\link{plot,predictor-method}},
  \code{\link{plot,validation-method}},
  \code{\link{plot3d}},
  \code{\link{plot3d,nvalidation-method}},
  \code{\link{plot3d,validation-method}},
  \code{\link{predict}},
  \code{\link{prediction-class}},
  \code{\link{predictor-class}},
  \code{\link{summary}},
  \code{\link{validate}},
  \code{\link{validation-class}},
  \code{\link{cancerclass-internal}}, 
  \code{\link{ilogit}}, 
  \code{\link{calc.roc}}, 
  \code{\link{calc.auc}}, 
  \code{\link{get.d}}, 
  \code{\link{get.d2}}, 
  \code{\link{get.prop}}, 
  \code{\link{get.ntrain}}, 
  \code{\link{prepare}}, 
  \code{\link{filter}} 
}  

\keyword{ package }
