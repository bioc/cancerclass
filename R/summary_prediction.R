
if (!isGeneric("summary")) {
  setGeneric("summary", function(object="prediction", ...) standardGeneric("summary"), package="cancerclass");  
}

setMethod("summary", signature(object="prediction"),    
  
  function(object = "prediction", positive=sort(object@prediction[, "class_membership"])[1], score=NULL, ...) 
  {
    index <- union(which(!is.na(object@prediction[, "class_membership"])), which(!is.na(object@prediction[, "class_predicted"])))
    object@prediction <- object@prediction[index, ]
    
    result <- list()
    result$class <- positive
    label <- sort(as.character(unique(object@prediction[, "class_membership"])))
    if (label[1] != positive) label <- label[2:1]
    label.predicted <- paste(label, "predicted", sep = "_")
    
    result$crosstab <- matrix(nrow=2, ncol=2)
    rownames(result$crosstab) <- paste(label, "predicted", sep="_")
    colnames(result$crosstab) <- label
    membership <- list() 
    predicted <- list()
  
  for (i in 1:2) {
    membership[[i]] <- which(object@prediction[, "class_membership"] == label[i])
    predicted[[i]] <- which(object@prediction[, "class_predicted"] == label[i])
  }   
  for (i in 1:2) for (j in 1:2) result$crosstab[i, j] <- length(intersect(predicted[[i]], membership[[j]]))

    #result$p
    result$p <- fisher.test(result$crosstab, alternative="greater")$p.value
    enrich <- rbind(result$crosstab[1, ], result$crosstab[1, ]  + result$crosstab[2, ])
    
    #result$enrich
    result$enrich <- vector() 
    result$enrich[1] <- get.prop(enrich[1, 1], enrich[1, 2])
    result$enrich[2] <- get.prop(enrich[2, 1], enrich[2, 2])
    result$enrich[3] <- signif(fisher.test(enrich, alternative="greater")$p.value, 2)
    names(result$enrich) <- c("PPV", "prevalence", "p")
    
    #result$tab
    result$tab <- matrix(ncol=3, nrow=4)
    colnames(result$tab) <- c(label, "PV")
    rownames(result$tab) <- c(label.predicted, "accuracy", "enrichment")    
    print(result$crosstab)
    result$tab[1:2, 1:2] <- result$crosstab
    result$tab[3, 1] <- get.prop(result$crosstab[1, 1], result$crosstab[2, 1])
    result$tab[3, 2] <- get.prop(result$crosstab[2, 2], result$crosstab[1, 2])
    result$tab[1, 3] <- get.prop(result$crosstab[1, 1], result$crosstab[1, 2])
    result$tab[2, 3] <- get.prop(result$crosstab[2, 2], result$crosstab[2, 1])
    for (i in 1:2) result$tab[i, 3] <- get.prop(result$crosstab[i, 1], result$crosstab[i, 2])
    result$tab[3, 3] <- get.prop(result$crosstab[1, 1] + result$crosstab[2, 2], result$crosstab[1, 2] + result$crosstab[2, 1])
    result$tab[4, ] <- result$enrich
    result$tab[3, 1] <- paste("sens. =", result$tab[3, 1])
    result$tab[3, 2] <- paste("spec. =", result$tab[3, 2])
    result$tab[1, 3] <- paste("PPV =", result$tab[1, 3])
    result$tab[2, 3] <- paste("NPV =", result$tab[2, 3])
    result$tab[3, 3] <- paste("acc. =", result$tab[3, 3])
    result$tab[4, 1] <- paste("PPV =", result$tab[4, 1])
    result$tab[4, 2] <- paste("preval. =", result$tab[4, 2])
    result$tab[4, 3] <- paste("p = ", result$tab[4, 3], " (p_2x2 = ", signif(result$p, 2), ")", sep="")
    
  ###
  result$sensitivity <- vector()  
  for (i in 1:2) result$sensitivity[label[i]] <- result$crosstab[i, i] / sum(result$crosstab[, i])

  result$PPV <- vector()  
  for (i in 1:2) result$PPV[label[i]] <- result$crosstab[i, i] / sum(result$crosstab[i, ])

  result$misclassification <- 1 - result$sensitivity
  result$misclassification["all"] <- (result$crosstab[1, 2] + result$crosstab[2, 1]) / sum(result$crosstab)
   
  result$AUC <- vector()
  result$S <- list()
  for (z in c("z", "zeta", "ratio")) {
    x <- as.numeric(object@prediction[, z])
    y <- rep(0, nrow(object@prediction))
    y[which(object@prediction[, "class_membership"] == positive)] <- 1
    S <- calc.roc(x, y)
    result$S[[z]] <- S
    result$AUC[z] <- calc.auc(S)
  }
  
  #this <- .Alias(object);
  #slot(this,"summary")<-result
  #slot(object,"summary")<-result
  #message(str(object@summary))
  #assign(object@summary, result);
  return(result)
}
)


#2.
#setGeneric("plot", function(x="prediction", y="missing", type = "xy", method = names(x@misclass), sig = c(0.95, NA, NA), xlog = FALSE, pos = "topleft", ntrain = "all", min.percent=5, n=nrow(x@nselected), col = c(grey(0.7), grey(0.3)), ylim = c(0, 100), cex.names=0.5, ...) standardGeneric("plot"), package="cancerclass");
#setMethod("plot", signature(x="prediction", y="missing", type="ANY", method="ANY", sig="numeric", xlog="ANY", pos="ANY", ntrain="numeric", min.percent="numeric", n="numeric", col="numeric", ylim="numeric", cex.names="numeric"),
#type, method, sig, xlog, pos, ntrain, min.percent, n, col, ylim, cex.names
#type = "xy", method = names(x@misclass), sig = c(0.95, NA, NA), xlog = FALSE, pos = "topleft", ntrain = "all", min.percent=5, n=nrow(x@nselected), col = c(grey(0.7), grey(0.3)), ylim = c(0, 100), cex.names=0.5, 
#setGeneric("plot", function(x, y,...) standardGeneric("plot"));

