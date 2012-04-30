\name{GOLUB}
\docType{data}

\alias{GOLUB}
\alias{GOLUB1}

\title{GOLUB DATA}

\description{ 
  Gene expression data from the leukemia
  microarray study of Golub et al. [1]. Dataset GOLUB has a dimention of
  7129 genes in 72 tumors samples. Dataset GOLUB1 has a dimention of 
  3571 genes in 72 tumors samples. This dataset is filtered and preprocessed
  as described in [2].
}

\usage{
  data(GOLUB)
  data(GOLUB1)
}

\value{
 Data and annotations are organized in a ExtressenSet of the package Biobase.
 \item{GOLUB}{ExpressionSet (7129 genes in 72 tumors) }
 \item{GOLUB1}{ExpressionSet (3571 genes in 72 tumors) }
}

\references{
  [1] Golub TR et al (1999), \emph{Molecular Classification of cancer: class Discovery and Class Prediction by gene expression monitoring}, Science 286:531-7.\cr
  [2] Dudoit S, Fridlyand J (2002), \emph{A prediction-based resampling method for estimating the number of clusters in a dataset}, Genome Biol. 3(7):RESEARCH0036.
}

\examples{
### nvalidate
data(GOLUB1)
nval <- nvalidate(GOLUB1[1:200, ])
# Use only the first 200 genes for speed-up of the calculations
plot(nval, type="xy")
plot(nval, type="genes")
plot(nval, type="samples")

### validate
data(GOLUB1)
val <- validate(GOLUB1[1:200, ])
# Use only the first 200 genes for speed-up of the calculations
plot(val, type="xy")
plot(val, type="genes")
plot(val, type="samples")

### fit und predict
data(GOLUB1)
train <- GOLUB1[, 1:38]
test <- GOLUB1[, 39:72]
predictor <- fit(train, method="welch.test")
prediction <- predict(predictor, test, positive="AML", ngenes=50, dist="cor")
plot(prediction, type="histogram", score="zeta")
plot(prediction, type="curves", score="zeta")
plot(prediction, type="roc", score="zeta")
summary(prediction)

### loo
data(GOLUB1)
cv <- loo(GOLUB1, positive="AML", ngenes=10, method="welch.test", dist="cor")
plot(cv, type="histogram", score="zeta")
plot(cv, type="samples", score="zeta")
plot(cv, type="curves", score="zeta")
plot(cv, type="roc", score="zeta")
}

\keyword{datasets}
