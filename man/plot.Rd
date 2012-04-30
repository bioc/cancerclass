\name{plot}
\docType{methods}
\alias{plot}
\alias{plot,validation-method}
\alias{plot,nvalidation-method}
\alias{plot,prediction-method}
\alias{plot,predictor-method}


\title{  Plot Method for 'validation, nvalidation, prediction, predictor' Classes }

\description{
  \itemize{
  \item class nvalidation: \cr
  Plot of misclassification rates in dependence of the number of features that were used for classification. 
  Total, class1 and class2 misclassification rates including confidence intervals can be plotted separately.
  \item class prediction: \cr
  Plot methods for continuous predictions scores. Prediction scores are obtained by validation 
  of a predictor in a test set. Four methods for assessment and visualization of predictor performance can be selected by the parameter \code{type}. 
  \item class predictor: \cr
  Among the three continuous prediction scores (\code{z}, \code{zeta} and \code{ratio}), \code{zeta} has the special property to be a linear combination of gene expression values. The plot method works only for the prediction score \code{zeta} and visualizes the contribution of each gene to the score.
  \item class validation: \cr
  Plots the misclassification rate in dependence of the training set size. Total, class1 and class2 misclassification rates can be plotted including confidence intervals.
  }
}

\usage{
  plot(x, y, ...)
}

\arguments{
  \item{x}{Object of class \code{nvalidation}.}
  \item{y}{missing}
  \item{\dots}{Further arguments directly passed to \code{plot}.}
}

%validation
\section{Methods}{
\describe{
  \item{x}{Object of class \code{validation}}
  \item{y}{missing}
  \item{type}{Three different kinds of plots can be generated: a xy-plot showing the misclassification rate in dependence of the training set size (type="xy"), a barplot showing the misclassification rates for each of the samples (type="samples"), a barplot showing how often (in \%) a gene is included in a predictor (type="genes"). }
  \item{method}{A character vector specifying the types of misclassification rates to be plotted. Possible types are the names of the classes and all for the total misclassification rate. }
  \item{anno}{Only relevant if type="genes": annotation of the features by array probes (anno="probe") or gene symbols (anno="symbol"). }
  \item{sig}{Vector of numerical values corresponding to the method vector. The numerical values are equal to the confidence level, if equal to NULL, the corresponding confidence interval is not plotted. }
  \item{xlog}{A logical value. If TRUE, a logarithmic scale is used for the x-axis. }
  \item{pos}{Position of legend specified by a keyword from the list "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center". }
  \item{ntrain}{Only relevant if type="samples" or type="genes". Either results for predictors trained in training sets of different sample sizes (ntrain="all"), or results for predictors trained in training sets of the sample size specified by ntrain. }
  \item{min.percent}{Only relevant if type="genes". Theshold for selection of the genes that are plotted. }
  \item{n}{Only relevant if type="genes": Number of top genes that are plotted. }
  \item{col}{Only relevant if type="samples". Color of the boxes for samples of class1 and of class2. } 
  \item{ylim}{Range of the y-axis. } 
  \item{cex.names}{Only relevant if type="samples" or type="genes". Scaling factor for the labels of the x-axis. }
  \item{\dots}{Further arguments directly passed to \code{plot}. }
}}

%nvalidation
\section{Methods}{
\describe{
  \item{x}{Object of class \code{nvalidation}}
  \item{y}{missing}  
  \item{type}{Three different kinds of plots can be generated: a xy-plot showing the misclassification rate in dependence of the training set size (type="xy"), a barplot showing the misclassification rates for each of the samples (type="samples"), a barplot showing how often (in \%) a gene is included in a predictor (type="genes"). }
  \item{method}{A character vector specifying the types of misclassification rates to be plotted. Possible types are the names of the classes and all for the total misclassification rate. }
  \item{anno}{Only relevant if type="genes": annotation of the features by array probes or gene symbols. }
  \item{sig}{Vector of numerical values corresponding to the method vector. The numerical values are equal to the confidence level, if equal to NULL, the corresponding confidence interval is not plotted. }
  \item{xlog}{A logical value. If TRUE, a logarithmic scale is used for the x-axis. }
  \item{pos}{Position of legend specified by a keyword from the list "bottomright", "bottom", "bottomleft", "left", "topleft", "top", "topright", "right" and "center". }
  \item{ngenes}{Only relevant if type="samples" or type="genes". Either results for predictors including different number of genes (ngenes="all"), or results for predictors including the number of genes specified by the numeric value ngenes. }
  \item{min.percent}{Only relevant if type="genes". Theshold for selection of the genes that are plotted. }
  \item{n}{Only relevant if type="genes": Number of top genes that are plotted. }
  \item{col}{Only relevant if type="samples". Color of the boxes for samples of class1 and of class2. }
  \item{ylim}{Range of the y-axis. }
  \item{cex.names}{Only relevant if type="samples" or type="genes". Scaling factor for the labels of the x-axis. }
  \item{\dots}{Further arguments directly passed to \code{plot}.}
}}

%prediction
\section{Methods}{
\describe{
  \item{x}{Object of class \code{prediction}}
  \item{y}{missing}        
  \item{type}{Four different kinds of plots can be generated: a histogram showing the distribution of the prediction score in class1 and class2 (type="histogram"), a xy-plot showing sensitivity, specificity, positive prediction value (PPV) and negative prediction value (NPV) in dependence on cutoffs for the prediction score (type="curves"), an ROC curve including calculation of the area under the curve (type="roc"), a barplot showing the prediction score for each of the samples (type="samples"). }
  \item{score}{Specification of the prediction scores used for the plot: score="z", "zeta" or "ratio".  If type="roc" a comparative analysis of two or three scores can be done. }
  \item{breaks.dist}{Distance of breaks. }
  \item{ci}{Only relevant, if type="roc". The method to calculate confidence intervals for sensitivity and specificity. Possibly values are  "exact", "ac", "asymptotic", "wilson", "prop.test", "bayes", "logit", "cloglog", "probit", see R package binom for details. }
  \item{col}{Only relevant if type="samples". Numerical or character vector of length two specifying the color of symbols for correct and wrong classifications. }
  \item{curves}{Only relevant if type="curves". A vector of strings specifying the curves that are included into the plot. Can include "sensitivity", "specificity", "PPV" and "NPV". }
  \item{col.curves}{Only relevant if type="curves". A vector of strings corresponding to curves. Specifies the colors of the plot curves. }
  \item{lty}{Only relevant if type="roc". Numerical or character vector corresponding to the vector score. Specifies the line types used for the ROC plot. }
  \item{npoints}{Only relevant if type="logistic". Number of points to be plotted. }
  \item{alpha}{Only relevant if type="logistic". The probability of class membership is estimated by logistic regression. The parameter alpha specifies the confidence level for the confidence interval of this probability. }
  \item{main}{Title of the plot. }
  \item{cex.names}{Only relevant if type="samples". Scaling factor for the labels of the x-axis. }
  \item{\dots}{Further arguments directly passed to \code{plot}.}
}}

%predictor
\section{Methods}{
\describe{
  \item{x}{Object of class \code{validation}}
  \item{y}{missing}
  \item{type}{Currently only type="genes" is supported. }
  \item{ngenes}{Number of genes in the predictor. }
  \item{dist}{Character string specifying the method for calculation of the distance between test samples and the centroids. Possible values are "euclidean", "angle", "cor", "center". }
  \item{anno}{Annotation of the features by array probes (anno="probe") or gene symbols (anno="symbol"). }
  \item{ylab}{Label of the y-axis. }
  \item{main}{Title of the plot. }  
  \item{\dots}{Further arguments directly passed to \code{plot}.}
}}

\seealso{ 
  %\code{\linkS4class{plot,validation-method}}, 
  %\code{\linkS4class{plot,nvalidation-method}},
  %\code{\linkS4class{plot,predictor-method}},
  %\code{\linkS4class{plot,prediction-method}},
  %\code{\link{plot}}, 
  \code{\link{validate}}, 
  \code{\linkS4class{validation}} 
}        

\examples{
### see: help(GOLUB);
}

\keyword{methods}
\keyword{classif}
