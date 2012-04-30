loo <- function (eset, class="class", method = "welch.test", ngenes=50, dist="cor", hparam = 0.75, positive="") 
{  
  eset=prepare(eset);
  classifier=as.character(pData(eset)$class)
  #eset=exprs(eset) 
  #ix=order(classifier) # sortierung der für z.b. wilcoxon maximal wichtig
  #eset=eset[,ix]  
 
    n <- length(classifier)
    validation.loo <- NULL
    
    for (i in 1:n) {
        message("\nSample: ", i)
        message("################################")

        predictor1 <- fit(eset[, -i], class=class, method = method, hparam = hparam)      
        prediction <- predict(predictor1, eset[, i, drop = FALSE], positive=positive, class=class, ngenes=ngenes, dist=dist)       
        #print(prediction)
                
        validation.loo <- rbind(validation.loo, prediction@prediction[1,])
        rownames(validation.loo)[nrow(validation.loo)] <- rownames(prediction@prediction)
    }
    colnames(validation.loo) <- colnames(prediction@prediction)
    #attr(validation.loo, "class") <- "prediction"
    
    #prediction=new("prediction", loo=validation.loo)
    #class=sort(unique(classifier))[1]
    positive=prediction@positive
    
    prediction=new("prediction");            
    slot(prediction,"type")="loo"
    #slot(prediction,"predictor")="";
    slot(prediction,"cl")=class
    slot(prediction,"method")=method
    slot(prediction,"ngenes")=ngenes
    slot(prediction,"dist")=dist 
    slot(prediction,"prediction")=validation.loo  
    slot(prediction,"positive")=positive
        
    ###SUMMARY
    #...
        
    return(prediction)
}
