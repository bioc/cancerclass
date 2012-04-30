if (!isGeneric("plot3d")) {
  setGeneric("plot3d", function(object, ...) standardGeneric("plot3d"));  
}

setMethod("plot3d", "validation",
    function (object, ...) 
    {
 g=object@nselected #rows:genes; cols:training set size
  h=apply(g,1, function(x) {sum(x)} )

  ### 3D PLOT
  y=c(1:ncol(g)) #training set size
  x=c(1:nrow(g)) #genes
  z=g
  rep=paste("validation repeat 1:",object@nrep);
  persp(x, y, z, xlab="genes", ylab="number of training set size", zlab=rep,
  theta = 30, phi = 30, expand = 0.5, col = "lightblue", main="absolute frequency of genes in signature")
    }
)
