
#.First.lib <- function(lib, pkg) {
.onLoad <- function(lib, pkg) {
if( !require(methods) ) stop("we require methods for package Foo")
where <- match(paste("package:", pkg, sep=""), search())
.initFoo(where)
}

.initFoo <- function(where) {
  #assign("FILTER", c("cor","welch.test","welch","student.test","student","wilcoxon.test","wilcoxon","foldchange","fc","copa","ort","os","shift","throw") );       
  #assign("FILTER.DECREASING", c(cor=TRUE,welch.test=TRUE,welch=TRUE,student.test=TRUE,student=TRUE,wilcoxon.test=TRUE,wilcoxon=TRUE,foldchange=TRUE,fc=TRUE,copa=TRUE,ort=TRUE,os=TRUE,shift=TRUE,throw=TRUE) );
  assign("DPI",800) #not found
  DPI<-800          #not found 
} #init end

get.d <- function(X, y, method="cor") {
  f <- paste("d", method, sep=".")
  L <- list()
  L[[2]] <- y
  X=as.matrix(X);

  #apply(X, 2, function(x){ L[[1]] <- x; do.call(f, L) }) #for each column of X
  ###cor
  if (method=="cor")
  {
    z=apply(X, 2, function(x,y) 
    {
      x <- x - mean(x)
      x <- x / sqrt(sum(x*x))
      y <- y - mean(y)
      y <- y / sqrt(sum(y*y))
      z <- sqrt(sum((x-y)^2))  
    }, y=y)
  }  
  ###eucledean
  if (method=="euclidean")
  { 
    z=apply(X, 2, function(x,y) 
    {   
      z <- sqrt(sum((x-y)^2))      
    }, y=y)
  }   
  ###angle        
  if (method=="angle")
  { 
    z=apply(X, 2, function(x,y) 
    {   
      x <- x / sqrt(sum(x*x))
      y <- y / sqrt(sum(y*y))
      z <- sqrt(sum((x-y)^2))      
    }, y=y)
  }  
  ###center      
  if (method=="center")
  { 
    z=apply(X, 2, function(x,y) 
    {   
      x <- x - mean(x)
      y <- y - mean(y)
      z <- sqrt(sum((x-y)^2))    
    }, y=y)
  }   
           
  return(z)
}
  

get.d2 <- function(X, x1, x2, method = "cor") {
    f <- paste("d2", method, sep = ".")
    L <- list()
    L[[2]] <- x1
    L[[3]] <- x2
    X = as.matrix(X)
    if (method == "cor") {
        z = apply(X, 2, function(x, x1, x2) {
            x1 <- x1 - mean(x1)
            x1 <- x1/sqrt(sum(x1 * x1))
            x2 <- x2 - mean(x2)
            x2 <- x2/sqrt(sum(x2 * x2))
            x <- x - mean(x)
            x <- x/sqrt(sum(x * x))
            c = 0.5 * (x1 + x2)
            d = x1 - x2
            d_norm = (sqrt(sum(d^2)))^2
            dist1 = (x - c) %*% d
            z = -2 * dist1/d_norm
        }, x1 = x1, x2 = x2)
    }
    if (method == "euclidean") {
        z = apply(X, 2, function(x, x1, x2) {
            c = 0.5 * (x1 + x2)
            d = x1 - x2
            d_norm = sqrt(sum(d^2))^2
            dist1 = (x - c) %*% d
            z = -2 * dist1/d_norm
        }, x1 = x1, x2 = x2)
    }
    if (method == "angle") {
        z = apply(X, 2, function(x, x1, x2) {
            x1 <- x1/sqrt(sum(x1 * x1))
            x2 <- x2/sqrt(sum(x2 * x2))
            x <- x/sqrt(sum(x * x))
            c = 0.5 * (x1 + x2)
            d = x1 - x2
            d_norm = (sqrt(sum(d^2)))^2
            dist1 = (x - c) %*% d
            z = -2 * dist1/d_norm
        }, x1 = x1, x2 = x2)
    }
    if (method == "center") {
        z = apply(X, 2, function(x, x1, x2) {
            x1 <- x1 - mean(x1)
            x2 <- x2 - mean(x2)
            x <- x - mean(x)
            c = 0.5 * (x1 + x2)
            d = x1 - x2
            d_norm = (sqrt(sum(d^2)))^2
            dist1 = (x - c) %*% d
            z = -2 * dist1/d_norm
        }, x1 = x1, x2 = x2)
    }
    return(z)
}
 
 

ilogit <- function(x) {
  y <- 1 / (1 + exp(-x))
  return(y)
}

calc.roc <- function(x, y, ci="wilson") {
  index <- intersect(which(!is.na(x)), which(!is.na(y)))
  x <- x[index]
  y <- y[index]
  n <- length(index)
  index.0 <- which(y == 0)
  index.1 <- which(y == 1)
  a <- mean(x[index.1])
  b <- mean(x[index.0])
  if (a > b) index <- order(-x)
  else index <- order(x)
  x <- x[index]
  y <- y[index]
  S <- matrix(ncol=5, nrow=n+1)
  colnames(S) <- c("threshold", "sensitivity", "specificity", "PPV", "NPV")
  rownames(S) <- 0:n
  for (i in 1:(n-1)) {
    y.P <- y[1:i]
    y.N <- y[(i+1):n]
    n.TP <- length(which(y.P == 1))
    n.FP <- length(which(y.P == 0))
    n.TN <- length(which(y.N == 0))
    n.FN <- length(which(y.N == 1))
    S[i+1, "threshold"] <- (x[i] + x[i+1])/2
    S[i+1, "sensitivity"] <- n.TP / (n.TP + n.FN)
    S[i+1, "specificity"] <- n.TN / (n.TN + n.FP)
    S[i+1, "PPV"] <- n.TP / (n.TP + n.FP)
    S[i+1, "NPV"] <- n.TN / (n.TN + n.FN)
  }
  if (a > b) {  
    S[1, ] <- c(x[1]+0.01, 0, 1, 0, 1)
    S[n+1, ] <- c(x[n]-0.01, 1, 0, 1, 0)
  }
  else {
    S[1, ] <- c(x[1]-0.01, 0, 1, 0, 1)
    S[n+1, ] <- c(x[n]+0.01, 1, 0, 1, 0)
  }
  index <- which(S[, "threshold"] %in% x)
  if (length(index) > 0) S <- S[-index, , drop=FALSE]
  S[, 2:5] <- 100 * S[, 2:5]
  m <- c(length(which(y == 1)), length(which(y == 0)))
  names(m) <- c("sensitivity", "specificity") 
  if ((length(ci) == 1) && (ci %in% c("exact", "ac", "asymptotic", "wilson", "prop.test", "bayes", "logit", "cloglog", "probit"))) {
    for (x in c("sensitivity", "specificity")) {
      CI <- binom.confint(S[, x]/100 * m[x], m[x], methods=ci)  
      S <- cbind(S, CI[, c("lower", "upper")] * 100)
      colnames(S)[(ncol(S)-1):ncol(S)] <- paste(x, c("lower", "upper"), sep="_")
    }
  }
  return(S)
} 

calc.auc <- function(S) {
  sen <- S[, "sensitivity"]
  spe <- S[, "specificity"]
  if (max(sen) <= 1) {
    sen <- sen * 100 
    spe <- spe * 100 
  } 
  n <- length(sen)
  x <- 100-spe
  y <- sen
  z <- y - x
  a <- sqrt(2)*x + z/sqrt(2)
  b <- z/sqrt(2) 
  a <- a/100
  b <- b/100
  auc <- try(integrate(function(x){approx(a, b, x)$y}, 0, sqrt(2), subdivisions=500))
  if (length(auc) == 1) auc <- 0
  else auc <- auc$value + 0.5
  return(auc)
}

get.prop <- function(n, m, ci="exact", conf.level=.9) {
  x <- binom.confint(n, m+n, method=ci, conf.level=conf.level)
  y <- vector()
  for (a in c("mean", "lower", "upper")) y[a] <- round(x[a]*100)
  result <- paste(y["mean"], "% (", y["lower"], "%-", y["upper"], "%)", sep="")
  return(result)
}

get.ntrain <- function(n1, n2, percent, n.min=5) {
  if (percent > 50) {
    n <- n1
    n1 <- n2
    n2 <- n
    percent <- 100-percent
  }
  else n <- 0
  m1 <- 5
  m <- round(100/percent * n.min)
  m2 <- m-m1
  result <- NULL
  while (m1 < n1 && m2 < n2) {
    result <- cbind(result, c(m1, m2))
    m1 <- m1 + 1
    m <- round(100/percent * m1)
    m2 <- m-m1
  }
  if (n > 0) result <- result[2:1, ]
  rownames(result) <- c("n1", "n2")
  colnames(result) <- result["n1", ] + result["n2", ]
  return(result)
}

prepare <- function(eset, class="class") {
  nsamples <- ncol(eset)
  ngenes <- nrow(eset)
  index <- which(!is.na(pData(eset)[[class]]))
  nmissing <- nsamples - length(index)
  if (nmissing > 0) {
    message("WARNING: ", class, " is not defined for all samples!")
    message("Removing ", nmissing, " samples from the data set.")
    eset <- eset[, index]
    nsamples <- ncol(eset)
  }
  X <- exprs(eset)
  n <- apply(X, 1, function(x){length(which(!is.finite(x)))})
  index <- which(n < nsamples)  
  nmissing <- ngenes - length(index)
  if (nmissing > 0) {
    message("WARNING: Data set contains missing genes!")
    message("Removing ", nmissing, " genes from the data set.")
    eset <- eset[index, ]
    ngenes <- nrow(eset)
  }
  X <- exprs(eset)
  A <- which(!is.finite(X), arr.ind=TRUE)
  if (nrow(A) > 0) {
    message("WARNING: Data set contains missing values!")
    message("Imputing ", nrow(A), " missing values.")
    for (i in 1:nrow(A)) {
      m <- mean(X[A[i, 1], ], na.rm=TRUE)
      X[A[i, 1], A[i, 2]] <- m
    }
  }
  exprs(eset) <- X
  return(eset)
}

get.lm <- function(v, ngenes=10, dist="cor") {
  x1 <- v@predictor[1:ngenes, 2]
  x2 <- v@predictor[1:ngenes, 3]
  if (dist == "cor") {
    x1 <- x1 - mean(x1)
    x1 <- x1/sqrt(sum(x1 * x1))
    x2 <- x2 - mean(x2)
    x2 <- x2/sqrt(sum(x2 * x2))
  }
  if (dist == "angle") {
    x1 <- x1/sqrt(sum(x1 * x1))
    x2 <- x2/sqrt(sum(x2 * x2))
  }
  if (dist == "center") {
    x1 <- x1 - mean(x1)
    x2 <- x2 - mean(x2)
  }
  d <- x1 - x2
  dnorm <- sum(d^2)
  x <- -2/dnorm * d
  names(x) <- names(x1)
  x["intercept"] <- (sum(x1^2) - sum(x2^2)) / dnorm
  model <- ""
  for (g in names(x)) {
    if (g == "intercept") model <- paste(model, round(x[g], 2)) 
    else model <- paste(model, round(x[g], 2), g, "+")
  } 
  model <- sub("^ ", "", model)
  result <- list()
  result$x <- x
  result$model <- model
  result$x1 <- x1
  result$x2 <- x2
  return(result)
}

filter <- function(method)
{
  FILTER = c("cor", "welch.test", "welch", "student.test", 
        "student", "wilcox.test", "wilcox", "foldchange", "fc", 
        "copa", "ort", "os", "shift", "throw")
  FILTER.DECREASING = c(cor = TRUE, welch.test = FALSE, welch = TRUE, 
        student.test = FALSE, student = TRUE, wilcox.test = FALSE, 
        wilcox = TRUE, foldchange = TRUE, fc = TRUE, copa = TRUE, 
        ort = TRUE, os = TRUE, shift = TRUE, throw = TRUE)
  fx = which(FILTER == method)

  FILTER=list(FILTER=FILTER,DECREASING=FILTER.DECREASING, fx=fx)

  return(FILTER);
}


