
if (!isGeneric("plot")) {
  #ok: setGeneric("plot", function(x="prediction", y="ANY", ...) standardGeneric("plot"), package="cancerclass");
  setGeneric("plot", function(x="prediction", y="ANY", ...) standardGeneric("plot"), package="cancerclass");  
}

setMethod("plot", signature(x="prediction", y="ANY"), 
    function (x, y, type = c("histogram", "curves", 
        "roc", "logistic"), score = c("zeta", "z", "ratio"), 
        breaks.dist = 0.1, ci = "wilson", pch = c(15, 16), col = c("green", 
            "red"), curves = c("sensitivity", "specificity", 
            "PPV", "NPV"), col.curves = c("blue", "cyan", "darkgreen", 
            "green"), lty = 1:3, npoints = 100, alpha = 0.05, 
        main = NULL, cex.names = 0.5, ...) 
    {
        positive = x@positive
        type <- type[1]
        class1 = positive
        class2 = setdiff(unique(x@prediction[, "class_membership"]), 
            positive)
        if (is.null(main)) 
            main <- paste(class1, "vs.", class2)
        type <- type[1]
        pvalidation = x@prediction
        if (type == "histogram") {
            negative <- setdiff(pvalidation[, "class_membership"], 
                positive)
            score <- score[1]
            x <- as.numeric(pvalidation[, score])
            y <- rep(0, nrow(pvalidation))
            y[which(pvalidation[, "class_membership"] == positive)] <- 1
            index0 <- which(y == 0)
            index1 <- which(y == 1)
            p <- t.test(x[index0], x[index1])$p.value
            par(lwd = 2)
            breaks <- seq(round(min(x), 1) - 0.1, round(max(x), 
                1) + 0.1, breaks.dist)
            counts0 <- hist(x[index0], plot = FALSE, breaks = breaks)$counts
            counts1 <- hist(x[index1], plot = FALSE, breaks = breaks)$counts
            counts <- c(counts0, counts1)
            hist(x[index0], col = col[1], xlim = range(x), ylim = c(0, 
                max(counts) * 1.4), breaks = breaks, main = main, 
                angle = 135, density = 10, xlab = score, cex.main=1, lwd = 2)
            hist(x[index1], add = TRUE, col = col[2], breaks = breaks, 
                density = 10, lwd = 2)
            legend("topleft", legend = c(negative, positive), 
                col = col, lwd = 2, box.lwd = 2, inset = 0.02, cex=0.8)
            legend("topright", legend = paste("p =", signif(p, 
                2)), lty = 0, bty = "n", inset = 0.02, cex=1)
        }
        if (type == "samples") {
            score <- score[1]
            index <- order(as.numeric(pvalidation[, score]))
            pvalidation <- pvalidation[index, ]
            scores <- as.numeric(pvalidation[, score])
            names(scores) <- rownames(pvalidation)
            colors <- rep("red", length(scores))
            index <- which(pvalidation[, "class_membership"] == 
                pvalidation[, "class_predicted"])
            if (length(index) > 0) 
                colors[index] <- "green"
            par(las = 2)
            barplot(scores, ylab = score, col = colors, cex.names = cex.names, 
                ...)
            legend("topleft", legend = c("right prediction", 
                "wrong prediction"), pch = 15, col = c("green", 
                "red"), inset = 0.02)
        }
        if (type == "curves") {
            score <- score[1]
            x <- as.numeric(pvalidation[, score])
            y <- rep(0, nrow(pvalidation))
            y[which(pvalidation[, "class_membership"] == positive)] <- 1
            S <- calc.roc(x, y, ci = "none")
            ncurves <- length(curves)
            plot(S[, "threshold"], S[, "sensitivity"], type = "n", 
                xlab = score, ylab = "percentage", ylim = c(0, 
                  100), main = main, cex.main = 1, lwd = 2, ...)
            for (i in 1:ncurves) lines(S[, "threshold"], S[, 
                curves[i]], col = col.curves[i], lwd = 2)
            legend("bottomleft", legend = curves, col = col.curves, 
                lty = 1, inset = 0.02, cex = 0.8)
        }
        if (type == "roc") {
            par(lwd = 2)
            plot(0, 0, type = "n", xlim = c(0, 100), ylim = c(0, 
                100), xlab = "1 - specificity", ylab = "sensitivity", 
                main = main, cex.main = 1, lwd = 1, ...)
            abline(0, 1, col = "gray50", lty = 2)
            label <- vector()
            nscore <- length(score)
            for (i in 1:nscore) {
                z <- score[i]
                x <- as.numeric(pvalidation[, z])
                y <- rep(0, nrow(pvalidation))
                y[which(pvalidation[, "class_membership"] == 
                  positive)] <- 1
                S <- calc.roc(x, y, ci = ci)
                auc <- calc.auc(S)
                if (auc < 0.99) 
                  auc = round(auc, 2)
                else auc = 1 - signif(1 - auc, 1)
                label[i] <- paste(z, ", AUC = ", auc, sep = "")
                if (nscore == 1 && ncol(S) > 5) {
                  for (j in 1:nrow(S)) {
                    s <- S[j, ]
                    par(lwd = 1)
                    lines(c(100 - s["specificity_lower"], 100 - 
                      s["specificity_upper"]), c(s["sensitivity"], 
                      s["sensitivity"]), col = col[1], lwd = 1)
                    lines(c(100 - s["specificity"], 100 - s["specificity"]), 
                      c(s["sensitivity_lower"], s["sensitivity_upper"]), 
                      col = col[2], lwd = 1)
                  }
                }
                par(lwd = 2)
                lines(100 - S[, "specificity"], S[, "sensitivity"], 
                  lty = lty[i], lwd = 2)
            }
            legend("bottomright", label, lty = lty, inset = 0.02, cex=0.8)
        }
        if (type == "logistic") {
            score <- score[1]
            x <- as.numeric(pvalidation[, score])
            y <- rep(0, nrow(pvalidation))
            y[which(pvalidation[, "class_membership"] == positive)] <- 1
            fit <- glm(y ~ x, family = binomial)
            print(summary(fit))
            delta.deviance <- fit$null.deviance - fit$deviance
            p <- 1 - pchisq(delta.deviance, df = 1)
            res <- predict(fit, se.fit = TRUE)
            y <- ilogit(res$fit)
            z <- qnorm(1 - alpha/2)
            y.l <- ilogit(res$fit - z * res$se.fit)
            y.u <- ilogit(res$fit + z * res$se.fit)
            r.x <- max(x) - min(x)
            xx <- seq(min(x), max(x), r.x/npoints)
            X <- as.data.frame(xx)
            colnames(X) <- "x"
            res.new <- predict(fit, newdata = X, se.fit = TRUE)
            yy <- ilogit(res.new$fit)
            yy.l <- ilogit(res.new$fit - z * res.new$se.fit)
            yy.u <- ilogit(res.new$fit + z * res.new$se.fit)
            index <- list()
            label <- sort(unique(pvalidation[, "class_membership"]))
            if (label[1] == positive) {
                label[1] <- label[2]
                label[2] <- class1
            }
            membership <- list()
            predicted <- list()
            for (i in 1:2) {
                membership[[i]] <- which(pvalidation[, "class_membership"] == 
                  label[i])
                predicted[[i]] <- which(pvalidation[, "class_predicted"] == 
                  label[i])
            }
            for (i in 1:2) {
                index[[i]] <- list()
                index[[i]][["TRUE"]] <- intersect(membership[[i]], 
                  predicted[[i]])
                index[[i]][["FALSE"]] <- intersect(membership[[i]], 
                  predicted[[setdiff(1:2, i)]])
            }
            model <- paste("y = ", signif(fit$coefficients, 3)["x"], 
                score, " + ", signif(fit$coefficients["(Intercept)"], 
                  3))
            main <- paste(model, " (p = ", signif(p, 2), ")", 
                sep = "")
            ylab <- paste("probability", positive)
            plot(x, y, type = "n", xlab = score, ylab = ylab, 
                ylim = c(-0.05, 1.05), main = main, cex.main = 0.8, 
                lwd = 2, ...)
            lines(xx, yy, lwd = 2)
            lines(xx, yy.l, lty = "dashed", lwd = 2)
            lines(xx, yy.u, lty = "dashed", lwd = 2)
            for (i in 1:2) {
                if (i == 1) 
                  yy <- -0.05
                if (i == 2) 
                  yy <- 1.05
                points(x[membership[[i]]], rep(yy, length(membership[[i]])), 
                  pch = "|")
            }
        }
    }
)