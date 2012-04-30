\name{plot3d}
\docType{methods}
\alias{plot3d}
\alias{plot3d,nvalidation-method}
\alias{plot3d,validation-method}

\title{ Plot3d method for 'validtion and 'nvalidation' classes }

\description{
  \itemize{
  \item nvalidation: \cr
  For each gene, the frequency that it is selected is plotted in dependence of the number of genes in the predictor. 
  \item validation: \cr
  For each gene, the frequency that it is selected is plotted in dependence of the training set size. 
  }
}

\usage{
  plot3d(object,...)
}

\arguments{
  \item{object}{Object of class nvalidation or validation}
  \item{\dots}{Further arguments directly passed to \code{plot}.}
}


\seealso{ 
  %\code{\link{plot3d,validation-method}}, 
  %\code{\link{plot3d,nvalidation-method}},
  \code{\link{plot3d}}, 
  \code{\link{persp}}  
}

\keyword{methods}
\keyword{classif}
