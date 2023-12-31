
if (!isGeneric("plot")) { 
  setGeneric("plot", function(x="validation", y="ANY", ...) standardGeneric("plot"), package="cancerclass");  
}

setMethod("plot", signature(x="validation", y="ANY"), 
    function (x, y, type = "xy", method = names(x@misclass), 
        anno = "symbol", sig = c(0.95, NA, NA), xlog = FALSE, 
        pos = "topleft", ntrain = "all", min.percent = 5, n = nrow(x@nselected), 
        col = c(grey(0.7), grey(0.3)), ylim = c(0, 100), cex.names = 0.5, col.curves=c("black", "red", "green"),
        ...) 
    {
        if (type == "genes") {
            if (ntrain == "all") X <- apply(x@nselected, 1, mean)
            else X <- x@nselected[, as.character(ntrain)]
            if (!(anno %in% colnames(x@fdata))) anno <- "probe"
            if (anno != "probe") names(X) <- x@fdata[,anno]
            X <- sort(X/x@nrep * 100, decreasing = TRUE)
            #print(X[1:20])
            n <- min(length(which(X >= min.percent)), n)
            X <- X[1:n]
            if (is.null(ylim)) ylim <- c(min(X), max(X))
            par(las = 2)
            barplot(X, ylim = ylim, xlab="genes", ylab = "predictors including gene (in %)", cex.names = cex.names, ...)
        }
        if (type == "samples") {
            if (ntrain == "all") 
                X <- apply(x@samples, 1, mean)
            else X <- x@samples[, as.character(ntrain)]
            X <- X * 100
            classes <- unique(as.character(x@classifier))
            index <- list()
            color <- vector()
            if (is.null(names(col))) 
                names(col) <- classes
            for (c in classes) {
                index[[c]] <- which(x@classifier == c)
                index[[c]] <- index[[c]][order(X[index[[c]]])]
            }
            X <- c(X[index[[1]]], X[index[[2]]])
            color <- c(rep(col[classes[1]], length(index[[1]])), rep(col[classes[2]], length(index[[2]])))
            if (is.null(ylim)) 
                ylim <- c(min(X), max(X))
            par(las = 2)
            barplot(X, ylim = ylim, xlab="samples", ylab = "missclassifications (in %)", 
                col = color, cex.names = cex.names, ...)
            legend("topleft", legend = classes, pch = 15, col = col, 
                inset = 0.02)
        }
        if (type == "xy") {
            nsig <- length(sig)
            if (nsig < 3) 
                sig = c(sig, rep(NA, 3 - nsig))
            method = as.vector(method)
            ix = x@ntrain
            if (length(ix) <= 10) 
                line.type = "b"
            else line.type = "l"
		col.curves=c("black", "red", "green")
            nmethod <- length(method)		
            for (i in nmethod:1) {
                name = method[i]
                error = x@misclass[[name]]
                m = apply(error, 2, mean)
                if (!is.na(sig[i])) {
                  q05 = apply(error, 2, function(x) quantile(x, 
                    1 - sig[i]))
                  q95 = apply(error, 2, function(x) quantile(x, 
                    sig[i]))
                }
                if (i == nmethod) {
                  if (xlog == "TRUE") 
                    plot(ix, m, ylim = c(0, 1), type = "n", log = "x", 
                      xlab = "Training set size", ylab = "Proportion of misclassifications", 
                      lwd = 2, ...)
                  else plot(ix, m, ylim = c(0, 1), type = "n", 
                    xlab = "Training set size", ylab = "Proportion of misclassifications", 
                    lwd = 2, ...)
                  abline(h = 0.5, col = "gray50", lty = 2)
                }
                if (!is.na(sig[i]) && sig[i] > 0) {
                  lines(ix, q05, col = i, type = "l", lty = 3, 
                    lwd = 2)
                  lines(ix, q95, col = i, type = "l", lty = 3, 
                    lwd = 2)
                }
                lines(ix, m, col = col.curves[i], type = line.type, lwd = 2)
            }
            legend.text <- method
            if (nmethod == 1) {
                if (!is.na(sig[1]) && sig[1] > 0) {
                  legend.text <- c(method, paste(legend.text, 
                    ", ", sig[1] * 100, "% CI", sep = ""))
                  legend(pos, legend = legend.text, col = 1, lty = c(1, 3), inset = 0.02, cex = 0.8)
                }
                else legend(pos, lenged = legend.text, col = 1, lty = 1, inset = 0.02, cex = 0.8)
            }
            else {
                for (i in 1:nmethod) if (!is.na(sig[i]) && sig[i] > 
                  0) 
                  legend.text[i] <- paste(legend.text[i], " with ", 
                    sig[i] * 100, "% CI", sep = "")
                legend(pos, legend = legend.text, col=col.curves, lty = 1, inset = 0.02, cex = 0.8)
            }
        }
    }
)