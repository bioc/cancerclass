
if (!isGeneric("predict")) {
  setGeneric("predict", function(object="predictor", ...) standardGeneric("predict"), package="cancerclass");  
}

setMethod("predict", signature(object="predictor"),    
  function(object, eset, positive="", class="class", ngenes = 50, dist = "cor", ...) 
  { 
   
### prepare data
    classifier = as.character(pData(eset)[,class])
    data = exprs(eset)
    predictor = object@predictor[1:ngenes, ]

### check parameters for validity
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
    }    

    if (length(positive)>0 & positive!=" " & positive!="")
    {    
      cl1=setdiff(colnames(object@predictor)[2:3],positive)
      cl2=positive
    } else
    {   
      cl1 = colnames(object@predictor)[2]
      cl2 = colnames(object@predictor)[3]
    } 
    ng1 = which(classifier == cl1)
    ng2 = which(classifier == cl2)
    if (length(which(is.na(classifier) == TRUE)) > 0) {
        message("The dataset contains NA values! Please prepare dataset correctly.")
        return(0)
    }
	
### calculate prediction
    message("Read and prepare data...")
    message("  Info (dataset): ", nrow(data),"x",ncol(data) )
    message("  #", cl1, "/#", cl2, ": ", length(ng1), "/", length(ng2), sep = "")
    message("\nStart analyse...")
    data_test = data[rownames(predictor), ]
### calculate continuous prediction scores
    dist1 = get.d(data_test, predictor[, cl1], method = dist)
    dist2 = get.d(data_test, predictor[, cl2], method = dist)
    dist3 = get.d(predictor[, cl1], predictor[, cl2], method = dist)
    z = (dist1 - dist2)/dist3
    tmp = which(dist1 < dist2)
    tmp0 = which(dist1 == dist2)
    e_cl1 = 1 - (length(intersect(tmp, ng1))/length(ng1))
    e_cl2 = (length(intersect(tmp, ng2))/length(ng2))
    zeta = get.d2(data_test, predictor[, cl2], predictor[, cl1], method = dist)
    ratio = log2(dist1 / dist2)
### calculate predicted class
    predicted_class = rep(cl2, length(classifier))
    predicted_class[tmp] = cl1
    predicted_class[tmp0] = NA
    true_class = classifier
    result = cbind(true_class, predicted_class, z, zeta, ratio)
    colnames(result) = c("class_membership", "class_predicted", "z", "zeta", "ratio")
    rownames(result) = colnames(data_test)
    
### prepare result data (class prediction)
    prediction=new("prediction");
    slot(prediction,"prediction")=result        
    slot(prediction,"type")="traintest"
    slot(prediction,"predictor")=slot(object,"predictor")
    slot(prediction,"cl")= slot(object,"cl") #class
    slot(prediction,"method")=slot(object,"method")
    slot(prediction,"ngenes")=ngenes
    slot(prediction,"dist")=dist       
    slot(prediction,"positive")=cl2       
    message("Finished!")
    return(prediction)
  }
)