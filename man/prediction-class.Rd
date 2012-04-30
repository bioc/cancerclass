\name{prediction-class}
\docType{class}
\alias{prediction-class}

\title{Class "prediction"}
\description{	
  An object of this class is returned by the function \code{fit.predictor} 
  and contains the results of a validation of a predictor on a test set.
}
\section{Objects from the Class}{
	 Objects of this class can be created by \code{new("prediction", ...)}. 
}
\section{Slots}{
	 \describe{	 
    \item{\code{type}:}{String specifying the type of validation: "traintest" for validation in a test set or "loo" for leave-one-out cross-validation. }    
    \item{\code{predictor}:}{The predictor that is validated, object of class \code{predictor}. }
    \item{\code{method}:}{Character string specifying the feature selection method. Possible values are "cor", "student.test", "welch.test", "wilcoxon.test", "foldchange", "copa", "os", "ort", "shift", "throw".}  
    \item{\code{ngenes}:}{Numerical vector containing the numbers features that are used for classification. }    
    \item{\code{dist}:}{Character string specifiying the method for calculation of the distance between test samples and the centroids. Possible values are "euclidean", "angle", "cor", "center". }
    \item{\code{prediction}:}{Matrix containing the prediction results. Each row represents a sample of the test data set. The first two columns contain the actual and the predicted class membership of the sample. Columns three, four and five contain the scores "z", "zeta" and "ratio". The score "z" is defined as difference of the distance of the sample from the class1 centroid minus the distance of the sample from the class2 centroid divided by the distance of the two centroids. For calculation of the score "zeta", the sample is orthogonally projected to the straight line through the centroids. Then the difference is calculated between the distance from the class1 and the class2 centroid. The score "ratio" is defined as the logarithm of the ratio of the distance of the sample from the class1 centroid divided by the correponding distance from the class2 centroid. For all three prediction scores, the sample is predicted to belong to class1, if the prediction score is negative, while the sample is predicted to belong to class2, if the prediction score is positive.}   
    \item{\code{cl}:}{The column in pData(eset) that was used to define the class membership. } 
    \item{\code{positive}:}{Character string or number specifying one of the two classes. This information is used for sensitivity and specificity calculations. }   
  }
}

\seealso{
	\code{\link{predict}}, \code{\link{plot}}, \code{\link{fit}} 
}
\examples{
showClass("prediction")
}
\keyword{classes}
