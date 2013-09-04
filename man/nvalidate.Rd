\name{nvalidate}
\alias{nvalidate}

\title{Classification in a multiple random validation protocol in pependence of the number of features used for predictor construction
}

\description{
  Feature selection and class prediction in a multiple random validation protocol as it was introduced in [1]. Misclassifications rates are calculated for predictors that include different numbers of features.   
}

\usage{
nvalidate(eset, class="class", ngenes = c(5, 10, 20, 50, 100, 200, 500, 1000), method = "welch.test", dist="cor", ntrain ="balanced", nrep = 200, hparam = 0.75)
}

\arguments{
  \item{eset}{Bioconductor ExpressionSet}
  \item{class}{Specification of the column in pData(eset) that contains the class information. }     
  \item{ngenes}{Numerical vector specifying the numbers features that are used for classification. }  
  \item{method}{Character string specifying the feature selection method. Possible values are "cor", "student.test", "welch.test", "wilcoxon.test", "foldchange", "copa", "os", "ort", "shift", "throw".}
  \item{dist}{Character string specifying the method for calculation of the distance between test samples and the centroids. Possible values are "euclidean", "angle", "cor", "center". }
  \item{ntrain}{One of the strings "balanced" or "prevalence" or a numeric vector specifying the number of  samples of class1 and the number of samples of class2 in the training sets. }
  \item{nrep}{The number of repetitions for each training set size.}
  \item{hparam}{Hyperparameter needed for some of the feature selection methods: Quantile used for the methods for copa, ort, os (e.g. 0.75, 0.95). Minimum number of samples in each class after applying shift/throw (only necessary for the feature selection methods: throw, shift). }
}

\details{
The matrix \code{exprs(eset)} contains the expression signatures of the patients in the columns.  The character vector \code{pData(eset)[[class]]} contains the class membership of each sample or patient. Only tow-class problems are supported.

The hyperparameter \code{hparam} describes the minimum number of 
samples in each class after applying \code{shift/throw}. For \code{copa}, \code{ort} and \code{os} the hyperparameter specifies the quantile that has to be exceeded in order to consider a sample as an outlier. Typical values are 0.75 (default), 0.9, 0.95.

Validation is implemented in a multiple random validation protocol [1]. For each training set size, \code{nrep} training sets are randomly drawn from the patients. Features are selected and the centroid is calculated for each of the two classes in feature space. The test samples are classified to the class with the nearest centroid. 

Four methods are available for calculation of the distance between test samples and the centroids: euclidean distance, euclidean distance after centering, angle and Pearson correlation. Calculation of distances is executed using the internal function \code{get.d}.

Feature selection, classification and validation are for predictors that include \code{ngenes} features. 
}

\value{
  A \code{nvalidation} object, see \code{nvalidation.object} for details. 
  Objects of this class have a method for the function \code{plot}.
}

\references{
[1] Michiels S, Koscielny S, Hill C (2005), \emph{Prediction of cancer outcome
with microarrays: a multiple random validation strategy}, Lancet 365:488-92. 
}

\seealso{
  \code{\linkS4class{nvalidation}} 
}

\examples{
### see: help(GOLUB);
}

\keyword{ classif }
