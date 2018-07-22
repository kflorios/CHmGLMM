#' Estimates the modelFit
#'
#' It acts on the Data data.frame, on the number of items Q and on the number of individuals
#'@param Data a data.frame with the data. 1st column id, 2nd column time, remaining Q columns are the y 0/1 values (Q items)
#'@param Q the number of items. Set this to four.
#'@param n number of individuals
#'@return Returns the model fitted
#'@export
estimateModelFit2 <- function(Data,Q,n) {
  
  require(lme4)
  
  modelFit <- try({
    pairs <- combn(Q, 2)
    P <- ncol(pairs)
    DD <- do.call(rbind, list(Data[1:3], Data[1:3]))
    DD$outcome <- gl(2, nrow(Data))
    Models <- vector("list", P)
    Sigma2YiRI   <- vector("list", P)
    Sigma2YiRS   <- vector("list", P)    
    Sigma2YiRIRS <- vector("list", P)
    Sigma2YjRI <- vector("list", P)    
    Sigma2YjRS <- vector("list", P)
    Sigma2YjRIRS <- vector("list", P)
    extraParam <- vector("list", P)
    
    ModelsOne <- vector("list", Q)
    Sigma2YiRIone <- vector("list", Q)
    Sigma2YiRSone <- vector("list", Q)
    Sigma2YiRIRSone <- vector("list", Q)
    extraParamOne   <- vector("list", Q)
    
    
    for (i in 1:P) {
      prs <- pairs[, i]
      yi <- Data[[paste("y", prs[1], sep = "")]]
      yj <- Data[[paste("y", prs[2], sep = "")]]
      DD$y <- c(yi, yj)
      # Fit model
      #Models[[i]] <- glmer(y ~ 0 + outcome + outcome:(time + group + time:group) + (0 + outcome + outcome:time | id),  #Florios
      #                     data = DD, family = binomial,verbose=T,control=glmerControl(optCtr=list(maxfun=3000)))               #Florios
      #Models[[i]] <- glmer(y ~ 0 + outcome + outcome:(time + group + time:group) + (0 + outcome + outcome:time | id), 
      #                     data = DD, family = binomial, nAGQ=3, verbose=T, control=list(maxIter=500, maxFN=1500))     #Florios
      ##Models[[i]] <- glmer(y ~ 0 + outcome + outcome:(time) + ( 0 + outcome + outcome:time | id ), 
      ##                     data = DD, family = binomial(logit), nAGQ=3, verbose=T, control=list(maxIter=500, maxFN=1500))     #Florios
      Models[[i]] <- glmer(y ~ 0 + outcome + outcome:(time) + ( 0 + outcome + outcome:time  | id ), 
                           data = DD, family = binomial(logit), nAGQ=1, verbose=T, control=glmerControl(optCtr=list(maxfun=1000)))     #Florios
      
      cp<-i
      #1
      #Florios K-L pair of items
      
      
      ###Sigma2YiRI[[cp]]   <-   VarCorr(Models[[i]])$id[1,1]
      ###Sigma2YiRS[[cp]]   <-   VarCorr(Models[[i]])$id[3,3]
      ###Sigma2YiRIRS[[cp]] <-   VarCorr(Models[[i]])$id[3,1]
      ###Sigma2YjRI[[cp]]   <-   VarCorr(Models[[i]])$id[2,2]
      ###Sigma2YjRS[[cp]]   <-   VarCorr(Models[[i]])$id[4,4]
      ###Sigma2YjRIRS[[cp]] <-   VarCorr(Models[[i]])$id[4,2]
      
      extraParam[[cp]]   <- VarCorr(Models[[i]])$id
      extraParam[[cp]]   <- nearPD(extraParam[[cp]])
      # ...
      #6
      
    }

    #ModelsOne
    
    DDD <- do.call(rbind, list(Data[1:3]))
    #DDD$outcome <- gl(1, nrow(Data))
    
    indic <- 1:Q
    for (i in 1:Q) {
      yi <- Data[[paste("y", indic[i], sep = "")]]
      DDD$y <- c(yi)  
      ModelsOne[[i]] <- glmer(y ~ time + ( time  | id ), 
                              data = DDD, family = binomial(logit), nAGQ=1, verbose=T, control=glmerControl(optCtr=list(maxfun=1000)))     #Florios
      
      ###Sigma2YiRIone[[i]]   <- VarCorr(ModelsOne[[i]])$id[1,1]
      ###Sigma2YiRSone[[i]]   <- VarCorr(ModelsOne[[i]])$id[2,2]
      ###Sigma2YiRIRSone[[i]] <- VarCorr(ModelsOne[[i]])$id[1,2]
      
      extraParamOne[[i]] <- VarCorr(ModelsOne[[i]])$id
      extraParamOne[[i]] <- nearPD(extraParamOne[[i]])
      
    }
    
    # Calculate Results
    
    DataRaw <- Data
    
    aveThetas2(Models,ModelsOne,Data,5,n,Q,extraParam,extraParamOne,1)
  }, TRUE)
  
  return(modelFit)
  
}