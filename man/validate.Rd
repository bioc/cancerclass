\name{validate}
\alias{validate}

\title{
  Classification in a Multiple Random Validation Protocol in Dependence of the Training Set Size
}

\description{
  Feature selection and class prediction in a multiple random validation protocol.
  Misclassifications rates are calculated for different sizes of the training set.
}

\usage{
validate(eset, class="class", ngenes = 50, method = "welch.test", dist="cor", ntrain = "balanced", nrep = 200, hparam = 0.75)
}

\arguments{
  \item{eset}{Bioconductor ExpressionSet}
  \item{class}{Specification of the column in \code{pData(eset)} that contains the class information. }     
  \item{ngenes}{Numerical vector specifying the numbers features that are used for classification. }  
  \item{dist}{Character string specifiying the method for calculation of the distance between test samples and the centroids. Possible values are "euclidean", "angle", "cor", "center".}
  \item{method}{Character string specifying the feature selection method. Possible values are "cor", "student.test", "welch.test", "wilcoxon.test", "foldchange", "copa", "os", "ort", "shift", "throw".}  
  \item{ntrain}{
  One of the strings "balanced" or "prevalence" or a numeric matrix that contains the numbers of training samples of the first class in the in first row and the numbers of training samples of the second class in the second row. } 
  \item{nrep}{The number of repeated training-test splits for each training set size. }
  \item{hparam}{Hyperparameter needed for some of the feature selection methods. For methods copa, ors and os: Quantile (e.g. 0.75, 0.9, 0.95) used in order to detect outliers. For methods shift and throw: the minimum number of samples in each class after applying shift or throw. }
}

\details{
The matrix \code{exprs(eset)} contains the expression signatures of the patients in the columns.  
The character vector \code{pData(eset)[[class]]} contains the class membership of each sample or patient. Only two-class problems are supported.

The hyperparameter \code{hparam} describes the minimum number of samples in each class after applying shift/throw. 
For \code{copa}, \code{ort} and \code{os} the hyperparameter specifies the quantile that has to be exceeded in order to consider a sample as an outlier. Typical values are 0.75 (default), 0.9, 0.95.
 
Validation is implemented in a multiple random validation protocol [1]. For each training set size, \code{nrep} training sets are randomly drawn from the patients. Features are selected and the centroid is calculated for each of the two classes in feature space. The test samples are classified to the class with the nearest centroid. 

Four methods are available for calculation of the distance between test samples and the centroids: euclidean distance, centered euclidean distance, angle and Pearson correlation. Calculation of distances is executed using the internal function \code{get.d}.

The parameter \code{ntrain} should be equal to one of the strings "balanced" or "prevalence" or a numeric matrix with two rows. For \code{ntrain = "balanced"}, a balanced layout is used, i.e. half of the training set is chosen from each of the two classes. For \code{ntrain = "prevalence"} the training sets are balanced according to the prevalence of the two classes in the entire data set. Further, the user can manually specify the sizes of the training sets. 
}

\value{
A \code{validation} object, see \code{validation.object} for details. 
Objects of this class have a method for the function \code{plot}. 
}

\references{
[1] Michiels S, Koscielny S, Hill C (2005), \emph{Prediction of cancer outcome
with microarrays: a multiple random validation strategy}, Lancet 365:488-92. 
}

\seealso{
  \code{\linkS4class{validation}}
}

\examples{
### see: help(GOLUB); 
}

\keyword{ classif }
