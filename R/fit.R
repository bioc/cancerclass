fit <- function (eset, class="class", method = "welch.test", hparam = 0.75) 
{         
	
### prepare and preprocess data
    eset = prepare(eset)      
    fdata = pData(featureData(eset))
    classifier = as.character(pData(eset)[, class])
    data = exprs(eset);
	
### check parameter for validity
    FILTER=filter(method)   
    if (length(FILTER$fx) == 0) { 
        message("Please select a correct method:\n", methods ,"\n")
        return(0)
    }  
    if (length(setdiff(unique(classifier), "NA")) != 2) {
        message("Please select a classifier including 2 class labels. \n")
        return(0)
    }
    if (length(colnames(data)) == 0) {
        message("No colnames found in data matrix...")
        colnames(data) = 1:ncol(data)
        names(classifier) = 1:length(classifier)
    } else {
        message("Colnames found in data matrix...")
        samples = is.element(colnames(data), names(classifier)[which(!is.na(classifier))])
        if (length(which(samples)) != 0) {
            data = data[, samples]
            classifier = classifier[samples]
        }
        is_na = which(is.na(classifier))
        not_na = which(!is.na(classifier))
        if (length(is_na) > 0) {
            message("Filtering Na values of data matrix and classifier...")
            classifier = classifier[not_na]
            data = data[, not_na]
        }
        if (length(table(classifier)) != 2) {
            message("Classifier contains no 2 labels intersected with colnames of data matrix. \n")
            return(0)
        }
    }    
    cl1 = sort(unique(classifier))[1]
    cl2 = setdiff(unique(classifier), cl1)
    ng1 = which(classifier == cl1)
    ng2 = which(classifier == cl2)
    if (length(which(is.na(classifier) == TRUE)) > 0) {
        message("The dataset contains NA values! Please prepare dataset correctly.")
        return(0)
    }
    
### calculate predictor
    message("Read and prepare data...")
    message("  Info (dataset): ", nrow(data), "x", ncol(data))
    message("  #", cl1, "/#", cl2, ": ", length(ng1), "/", length(ng2))
    message("  Method: ", method)
    message("\nStart analyse...")   
    training = data[, c(ng1, ng2)]
    n1 = as.double(length(ng1))
    n2 = as.double(length(ng2))	
### call C function for calculation of statistics and p-values
    stat = .C("statistics", (FILTER$fx), as.double(as.matrix(t(training))), as.double(nrow(training)), as.double(n1), as.double(n2), as.double(hparam), result = vector(length = nrow(training), mode = "double"), PACKAGE = "cancerclass")$result
    stat = abs(stat)
    myo = sort(stat, decreasing=FILTER$DECREASING[method], index.return = T)
    best = myo$ix
    res = matrix(nrow = nrow(eset), ncol = 3)

### prepare result data (class predictor)
    colnames(res) = c(method, cl1, cl2)
    rownames(res) = rownames(data)[best]
    res[, 1] = myo$x
    res[, 2] <- rowMeans(data[best, ng1])
    res[, 3] <- rowMeans(data[best, ng2])
    predictor=new("predictor")
    slot(predictor, "predictor") = res    
    slot(predictor, "cl") = class
    slot(predictor, "method") = method
    slot(predictor, "fdata") = fdata
    slot(predictor, "hparam") = hparam 
    message("Finished!")
    return(predictor)
}
