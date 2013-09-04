nvalidate <- function (eset, class = "class", ngenes = c(5, 10, 20, 50, 100, 
    200, 500, 1000), method = "welch.test", dist = "cor", ntrain = "balanced", 
    nrep = 200, hparam = 0.75) 
{

### prepare and preprocess data
    ngenes <- intersect(2:nrow(eset), ngenes)   
    eset = prepare(eset)
    fdata = pData(featureData(eset))
    classifier = as.character(pData(eset)[, class])
    eset = exprs(eset)
	
### check parameters for validity
    FILTER=filter(method)   
    if (length(FILTER$fx) == 0) {
        message("Please select a correct method:\n", methods, "\n")
        return(0)
    }
    data = as.matrix(eset)
    geneix = ngenes
    rep = nrep
    classifier = as.vector(classifier)
    if (length(setdiff(unique(classifier), NA)) != 2) {
        message("Please select a classifier including 2 class labels. \n")
        return(0)
    }
    if (length(colnames(data)) == 0) {
        message("No colnames found...")
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
    if (length(ntrain) == 1 && ntrain == "balanced") {
        n <- min(nng1, nng2)
        ntrain <- c(round(2/3 * n), round(2/3 * n))
        names(ntrain) <- c(cl1, cl2)
    }
    if (length(ntrain) == 1 && ntrain == "prevalence") {
        ntrain <- c(round(2/3 * nng1), round(2/3 * nng2))
        names(ntrain) <- c(cl1, cl2)
    }
    if (length(ntrain) == 1 && is.numeric(ntrain)) {
        ntrain2 <- round(ntrain/2)
        ntrain <- c(ntrain2, ntrain2)
        names(ntrain) <- c(cl1, cl2)
    }
    if ((!is.numeric(ntrain)) || (length(ntrain) != 2)) {
        message("Warning: ntrain must be a vector of two numerical numbers or one of the strings \"balanced\" or \"prevalence\". \n")
        return(0)
    }
    if (max(ntrain) >= ncol(data)) {
        message("The ntrain vector contains a bigger trainingset size than possible. Max: ", ncol(data) - 1, "\n")
        return(0)
    }
    if (min(ntrain) < 3) {
        message("Warning: Make sure that the training set contains at least three patients of each class.\n")
    }
    if (max(geneix) > nrow(data)) {
        message("The ngenes vector contains more genes than available in dataset. \n")
        return(0)
    }
    if (length(which(is.na(data) == TRUE)) > 0) {
        message("The dataset contains NA values! Please prepare dataset correctly.")
        return(0)
    }
	
### init result variables
    n1 <- ntrain[cl1]
    n2 <- ntrain[cl2]
    n <- n1 + n2
    message("Read and prepare data...")
    message("  Info (dataset): ", nrow(data),"x",ncol(data) )
    message("  #", cl1, "/#", cl2, ": ", length(ng1), "/", length(ng2), sep = "")
    message("  Training set size: ", n1, "+", n2)
    message("  Method: ", method)
    message("\nStart analyse...")
    count = 0
    values_mean = matrix(0, nrow = rep, ncol = length(geneix))
    values_class1 = matrix(0, nrow = rep, ncol = length(geneix))
    values_class2 = matrix(0, nrow = rep, ncol = length(geneix))
    sel_genelist = matrix(0, ncol = length(geneix), nrow = nrow(data))
    misclass = matrix(0, nrow = ncol(data), ncol = length(geneix))
    selection = matrix(0, nrow = ncol(data), ncol = length(geneix))
    all_genes = c(1:nrow(data))
    number_of_samples = nrow(data)
	
### multiple random validation protocol
    for (topgenes in geneix) {
        count = count + 1
        ideal = c(rep(-1, n), rep(1, n))
        message(count, " of ", length(geneix), " (top genes: ", topgenes, ")")
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
                stat = .C("statistics", (FILTER$fx), as.double(as.matrix(t(training))), 
                  as.double(nrow(training)), as.double(n1), as.double(n2), 
                  as.double(hparam), result = vector(length = nrow(training), 
                    mode = "double"), PACKAGE = "cancerclass")$result
                stat = abs(stat)
                best <- order(stat, decreasing = FILTER$DECREASING[method])[1:topgenes]
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
            tmp = class2 < class1
            class2_test = te.g[tmp]
            class1_test = te.g[!tmp]
### calculate class prediction, calculate sensitivity and specificity			
            class1_real = setdiff(ng1, tr.g)
            class2_real = setdiff(ng2, tr.g)
            mean1 = 1 - (length(intersect(class1_real, class1_test))/(length(class1_real)))
            mean2 = 1 - (length(intersect(class2_real, class2_test))/(length(class2_real)))
            mean_all = ((length(class1_real) * mean1) + (length(class2_real) * 
                mean2))/length(c(class1_real, class2_real))
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
	
### prepare result data (class nvalidation)
    misclass = 1 - (misclass/selection)
    colnames(misclass) = geneix
    colnames(values_mean) = geneix
    colnames(values_class1) = geneix
    colnames(values_class2) = geneix
    rownames(misclass) = colnames(data)
    rownames(values_mean) = c(1:nrep)
    rownames(values_class1) = c(1:nrep)
    rownames(values_class2) = c(1:nrep)
    error = list(values_mean, values_class1, values_class2)
    names(error) = c("all", cl1, cl2)
    rownames(sel_genelist) = rownames(data)
    colnames(sel_genelist) = ngenes
    validation = new("nvalidation")
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