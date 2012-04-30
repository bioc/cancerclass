\name{summary,prediction-methods}
\docType{methods}
\alias{summary}
%\alias{summary-methods}
%\alias{summary,ANY-method}
\alias{summary,prediction-method}

\title{  Summary Method for 'prediction' Class }

\description{
  Assessment of the performance of a predictor on a test data set.
}

\usage{
  summary(object="prediction", ...)
}

\arguments{
  \item{object}{Object of class \code{prediction}.}
  \item{\dots}{Further arguments desribed below.}
}

\section{Methods}{
\describe{
  \item{object}{Object of class \code{prediction}.}   
  \item{positive}{Character string specifying one of the two classes. Sensitivity and specificity calculations are carried out with respect to this class.}
}}

\seealso{
  \code{\linkS4class{prediction}} 
}

\examples{         
### see: help(GOLUB);
}

\keyword{methods}


