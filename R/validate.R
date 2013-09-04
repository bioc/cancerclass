validate <- function (eset, class = "class", ngenes = 50, method = "welch.test", 
    dist = "cor", ntrain = "balanced", nrep = 200, hparam = 0.75) 
{

### prepare and preprocess data
    eset = prepare(eset)
    fdata = pData(featureData(eset))
    classifier = as.character(pData(eset)[, class])
    eset = exprs(eset)
    FILTER=filter(method)

### check parameters for validity
    if (length(FILTER$fx) == 0) {
        message("Please select a correct method:\n", methods, "\n")
        return(0)
    }
    data = as.matrix(eset)
    topgenes = ngenes
    rep = nrep
    classifier = as.vector(classifier)
    if (length(setdiff(unique(classifier), NA)) != 2) {
        message("Please select a classifier including 2 class labels. \n")
        return(0)
    }
    if (length(colnames(data)) == 0) {
        message("No colnames found in data matrix...")
        colnames(data) = 1:ncol(data)
        names(classifier) = 1:length(classifier)
    }
    else {
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
    cl1 = unique(classifier)[1]
    cl2 = unique(classifier)[2]
    ng1 = which(classifier == cl1)
    ng2 = which(classifier == cl2)
    nng1 <- length(ng1)
    nng2 <- length(ng2)
    if (ntrain %in% c("balanced", "prevalence")) {
        if (ntrain == "balanced") 
            balance <- 50
        if (ntrain == "prevalence") 
            balance <- nng1/(nng1 + nng2) * 100
        ntrain <- get.ntrain(nng1, nng2, balance)
        rownames(ntrain) <- c(cl1, cl2)
    }
    else balance <- ""
    if (!is.numeric(ntrain) || (!is.matrix(ntrain)) || (nrow(ntrain) != 2)) {
        message("ntrain must be a numeric matrix with two rows or one of the strings \"balanced\" or \"prevalence\".")
        return(0)
    }
	if (length(setdiff(c(cl1, cl2), rownames(ntrain))) > 0) {
		message(paste("The rownames of ntrain should refer to the classes \"", cl1, "\" and \"", cl2, "\".", sep=""))
		return(0)
	}
    if (length(which(is.na(data) == TRUE)) > 0) {
        message("The dataset contains NA values! Please prepare dataset correctly.")
        return(0)
    }
    nrun <- ncol(ntrain)
	
### init result variables
    message("Read and prepare data...")
    message("  Info (dataset): ", nrow(data),"x",ncol(data) )
    message("  #", cl1, "/#", cl2, ": ", length(ng1), "/", length(ng2), sep = "")
    message("  Method: ", method, sep = "")
    message("\nStart analyse...")
    values_mean = matrix(0, nrow = rep, ncol = nrun)
    values_class1 = matrix(0, nrow = rep, ncol = nrun)
    values_class2 = matrix(0, nrow = rep, ncol = nrun)
    sel_genelist = matrix(0, ncol = nrun, nrow = nrow(data))
    misclass = matrix(0, nrow = ncol(data), ncol = nrun)
    selection = matrix(0, nrow = ncol(data), ncol = nrun)
    all_genes = 1:nrow(data)
    number_of_samples = nrow(data) 
	
### multiple random validation protocol
    for (count in 1:nrun) {   
        n1 <- ntrain[cl1, count]	 
        n2 <- ntrain[cl2, count]
        n <- n1 + n2
        message(count, " of ", nrun, " (training set size: ", n1, "+", n2, ")")		
        genelist = rep(0, nrow(data))
### for each random set/repeat 
        for (j in 1:rep) {
            tmpstat = rep(0, nrow(data))
            tr.g = c(sample(ng1, n1), sample(ng2, n2))
            te.g = setdiff(c(ng1, ng2), tr.g)
            training = data[, tr.g]
            test = data[, te.g]

### collect list of discriminating genes			
            if (topgenes < number_of_samples) {
### call C function for calculation of statistics and p-values
               stat = .C("statistics", FILTER$fx, as.double(as.matrix(t(training))), 
               as.double(nrow(training)), as.double(n1), as.double(n2), 
               as.double(hparam), result = vector(length = nrow(training), 
               mode = "double"), PACKAGE = "cancerclass")$result
               stat = abs(stat)
               myo = sort(stat, decreasing = FILTER$DECREASING[method], index.return = T)
               best = myo$ix[1:topgenes]
            }
            else {
               best = all_genes
            }
            genelist[best] = genelist[best] + 1		  
            class1_samples = training[, 1:n1]		
            class2_samples = training[, (n1 + 1):n]
### calculate centroids for both classes			           
	    class1_mean <- rowMeans(class1_samples[best, ])
	    class2_mean <- rowMeans(class2_samples[best, ])			
	    test_samples = test[best, ]
### calculate distance of classes by the selected distance method
            class1 = get.d(test_samples, class1_mean, method = dist)
            class2 = get.d(test_samples, class2_mean, method = dist)
### calculate class prediction, calculate sensitivity and specificity
            tmp = class2 < class1
            class2_test = te.g[tmp]
            class1_test = te.g[!tmp]
            class1_real = setdiff(ng1, tr.g)
            class2_real = setdiff(ng2, tr.g)
            mean1 = 1 - (length(intersect(class1_real, class1_test))/(length(class1_real)))
            mean2 = 1 - (length(intersect(class2_real, class2_test))/(length(class2_real)))
            mean_all = ((length(class1_real) * mean1) + (length(class2_real) *  mean2))/length(c(class1_real, class2_real))
            values_mean[j, count] = mean_all
            values_class1[j, count] = mean1
            values_class2[j, count] = mean2
            i_cl = c(intersect(class1_real, class1_test), intersect(class2_real, class2_test))
            i_real = c(class1_real, class2_real)
            misclass[i_cl, count] = misclass[i_cl, count] + 1
            selection[i_real, count] = selection[i_real, count] + 1
        }
### save list of discriminating genes
        sel_genelist[, count] = genelist
    }
	
### prepare result data (class validation)
    n <- ntrain[1, ] + ntrain[2, ]
    misclass = 1 - (misclass/selection)
    colnames(misclass) = n
    colnames(values_mean) = n
    colnames(values_class1) = n
    colnames(values_class2) = n
    rownames(misclass) = colnames(data)
    rownames(values_mean) = 1:nrep
    rownames(values_class1) = 1:nrep
    rownames(values_class2) = 1:nrep
    error = list(values_mean, values_class1, values_class2)
    names(error) = c("all", cl1, cl2)
    rownames(sel_genelist) = rownames(data)
    colnames(sel_genelist) = n
    validation = new("validation")
    slot(validation, "ngenes") = ngenes
    slot(validation, "method") = method
    slot(validation, "dist") = dist
    slot(validation, "ntrain") = n
    slot(validation, "nrep") = nrep
    slot(validation, "hparam") = hparam
    slot(validation, "misclass") = error
    slot(validation, "nselected") = sel_genelist
    slot(validation, "samples") = misclass
    slot(validation, "classifier") = classifier
    slot(validation, "fdata") = fdata
    message("Finished!")
    return(validation)
}