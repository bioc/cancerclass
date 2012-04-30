\name{predictor-class}
\docType{class}
\alias{predictor-class}

\title{Class "predictor"}
\description{	
An object of this class is returned by the function \code{predict(...)} 
and represents a nearest centroid predictor learned on a training data set.   
}
\section{Objects from the Class}{
  Objects can be created by the call \code{new("predictor", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{predictor}:}{Object of class \code{"predictor"} }
    \item{\code{cl}:}{Character string or number specifying one of the two classes. }
    \item{\code{method}:}{Character string specifying the feature selection method that was used for predictor construction: Possible values are "cor", "student.test", "welch.test", "wilcoxon.test", "foldchange", "copa", "os", "ort", "shift", "throw".  }   
    \item{\code{hparam}:}{Hyperparameter needed for some of the feature selection methods. For methods copa, ors and os: Quantile (e.g. 0.75, 0.9, 0.95) used in order to detect outliers. For methods shift and throw: the minimum number of samples in each class after applying shift or throw. }
    \item{\code{type}:}{Type of calculation. Loo (leave one out) or prediction (simple prediction). }
    \item{\code{fdata}:}{Feature annotations inherited from the training ExpressionSet.}
  }
}

\section{Methods}{
  \code{plot}{A plot method is only avaiable for the score "zeta" that is a linear combination of features.} 
}

\seealso{
	See Also as \code{\link{fit}}, 
	or \code{\linkS4class{predictor}} for links to other classes
}
\examples{
showClass("predictor")
}
\keyword{classes}
