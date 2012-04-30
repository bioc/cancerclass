
if (!isGeneric("plot")) {
  setGeneric("plot", function(x="predictor", y="ANY", ...) standardGeneric("plot"), package="cancerclass");  
}
 
setMethod("plot", signature(x="predictor", y="ANY"),  
    function (x, y, type="genes", ngenes=10, dist="cor", anno="symbol", ylab="contribution to zeta", main=NULL, ...)
    {
      model <- get.lm(x, ngenes=ngenes, dist=dist)
      contrib <- sort(model$x[1:ngenes])
      if (!(anno %in% colnames(x@fdata))) anno <- "probe"
      if (anno != "probe") names(contrib) <- x@fdata[names(contrib), anno]
      barplot(contrib, las=2, lwd=2, ylab=ylab, main=main, ...)
    }
)