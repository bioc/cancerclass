\name{predict-methods}
\docType{methods}
\alias{predict}
\alias{predict-methods}
\alias{predict,predictor-method}

\title{  Predict Method for 'predictor' Class}

\description{
  Assessment of the performance of a predictor in a test data set. 
}

\usage{
  predict(object="predictor", ...)
}

\arguments{
  \item{object}{Object of class \code{predictor}.}     
  \item{\dots}{Further arguments described below.}
}

\section{Methods}{
\describe{
  \item{object}{Object of class \code{predictor}.}
  \item{eset}{Test set, stored in a Bioconductor ExpressionSet object.}   
  \item{positive}{String refering to one of the two classes. Sensitivity and specificity calculations are carried out with respect to this class.}   
  \item{class}{Specification of the column in \code{pData(eset)} that contains the class information. The default value is class="class".}   
  \item{ngenes}{Number of features used for classification. The default value is ngenes=50.}   
  \item{dist}{Character string specifiying the method for calculation of the distance between test samples and the centroids. Possible values are "euclidean", "angle", "center" and the default "cor".}     
}}

\section{Details}{
The test samples are classified to the class with the nearest centroid. 
For methods are available for calculation of the distance between 
test samples and the centroids: Euclidean distance, euclidean distance after centering, angle and Pearson correlation. Calculation of distances is executed using the internal function \code{get.d}.
}

\value{
  A \code{prediction} object, see \code{prediction.object} for details. 
  Objects of this class have a method for the function \code{plot}.
}

\seealso{
  \code{\link{fit}} 
  \code{\linkS4class{prediction}} 
  \code{\linkS4class{predictor}} 
}

\examples{
### see: help(GOLUB);
}

\keyword{ methods }
