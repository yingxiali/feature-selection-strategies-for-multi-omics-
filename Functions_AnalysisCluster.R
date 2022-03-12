# This function performs one iteration of the five times repeated 
# 5-fold cross-validation on a specific data set using a specific
# method for the analysis of the multi-blocks case.
#
#
# It takes the whole number 'iter', which corresponds to the iter-th line 
# of 'scenariogrid', which contains the necessary information
# on the iter-th setting.



evaluatesetting <- function(iter) {
  
  # Load some packages
  library(FSelector)
  #library(rJava)
  library(mRMRe)
  library(e1071)
  library(ranger)
  library(measures)#brier.score
  
#### 1  
  # Initiate lists in which the results
  # will be stored:
  
  runningtime <- list()
  
  auc_svm <- list()
  brier_svm <- list()
  accuracy_svm <- list()
  
  auc_rf <- list()
  brier_rf <- list()
  accuracy_rf <- list()
  
  
#### 2  
  # Obtain information for the iter-th setting:
  #iter = 1
  dat <- scenariogrid$dat[iter]
  seed <- scenariogrid$seed[iter]
  featuremethod <- scenariogrid$featuremethod[iter]
  nselvars <- scenariogrid$nselvars[iter]
  selectseparately <- scenariogrid$selectseparately[iter]
  alwaysincludclin <- scenariogrid$alwaysincludclin[iter]
  
  cvind <- scenariogrid$cvind[iter]
  cvfoldind <- scenariogrid$cvfoldind[iter]
  
#### 3  
  # Load data set:
  
  load(paste("./VarselCompStudy/Additional_file/Data/", dat, sep=""))
  
#### 4  
  # Make covariate matrix and target variable:
  
  blocknames <- c("clindata", 
                  "mirnadata",
                  "mutationdata",
                  "cnvdata",
                  "rnadata")
  
  blocknamesav <-  c("")
  
  if(class(try(ncol(clindata)))!="try-error")
    blocknamesav <- c(blocknamesav, "clindata")
  if(class(try(ncol(mirnadata)))!="try-error")
    blocknamesav <- c(blocknamesav, "mirnadata")
  if(class(try(ncol(mirnadata)))!="try-error")
    blocknamesav <- c(blocknamesav, "mutationdata")
  if(class(try(ncol(cnvdata)))!="try-error")
    blocknamesav <- c(blocknamesav, "cnvdata")
  if(class(try(ncol(rnadata)))!="try-error")
    blocknamesav <- c(blocknamesav, "rnadata")  
   #block[[1]]
   #block[[2]]
  
  blocknamesav <- blocknamesav[-1]
  
  eval(parse(text=paste("X <- cbind(", paste(blocknamesav, collapse=", "), ")", sep="")))
  
  eval(parse(text=paste("block <- rep(1:length(blocknamesav), times=c(", paste(paste("ncol(", blocknamesav, ")", sep=""), collapse=", "), "))", sep="")))
  
  block <- lapply(1:length(blocknamesav), function(x) which(block==x))
 
  #y <- cbind(targetvar$time, targetvar$status)
  y <- TP53_mutation

#### 5 
  # Divide data set into training and test data:

  ncv <- 5
  set.seed(seed)

  cvdiv <- makeCVdiv(n=nrow(X), ncv=ncv)
  
  Xtrain <- X[cvdiv!=cvfoldind,] # dataframe
  ytrain <- y[cvdiv!=cvfoldind,]
  
  Xtest <- X[cvdiv==cvfoldind,] # vector
  ytest <- y[cvdiv==cvfoldind,]
  
#### 6  
  #feature selection
  
  timestart <-Sys.time()
  
  if(featuremethod=="t-test") {
    
    index <- t_index(Xtrain, ytrain, block, nselvars, selectseparately, alwaysincludclin)
    
  }else if(featuremethod=="chi-squared"){
    
    index <- chi_index(Xtrain, ytrain, block, nselvars, selectseparately, alwaysincludclin)
    
  }else if(featuremethod=="infor"){
  
    index <- infor_index(Xtrain, ytrain, block, nselvars, selectseparately, alwaysincludclin)
    
  }else if(featuremethod=="relief"){
    
    index <- relief_index(Xtrain, ytrain, block, nselvars, selectseparately, alwaysincludclin)
    
  }else if(featuremethod=="mrmr"){
    
    index <- mrmr_index(Xtrain, ytrain, block, nselvars, selectseparately, alwaysincludclin)
    
  }else if(featuremethod=="rf"){
    
    index <- rf_index(Xtrain, ytrain, block, nselvars, selectseparately, alwaysincludclin)
    
  }
  
  timeend <-Sys.time()
  runningtime[[1]] <-timeend-timestart
  
#### 7 
  #classifiers
  data_train <- cbind(ytrain, Xtrain[,index])
  data_test <- cbind(ytest, Xtest[,index])
  #rf
  fit <- ranger(ytrain ~ ., data = data_train,
                  num.threads = 1,
                #gamma = g,cost=c,
                probability=1)
  pred <- predict(fit, data = data_test, probability=1)
  yh  <- pred[["predictions"]]
  positive = 1
  negative = 0
  auc_rf[[1]] <- AUC(yh[,"1"], data_test$ytest, negative, positive)
  brier_rf[[1]] <- Brier(yh[,"1"], data_test$ytest, negative, positive)
  response = as.factor(as.numeric(yh[,"1"] > 0.5))
  accuracy_rf[[1]] <- ACC(data_test$ytest, response)
    
  #SVM
  fit <- svm(ytrain ~ ., data = data_train,
             type = 'C-classification', kernel = 'radial',
             #gamma = g,cost=c,
             probability=1)
  pred <- predict(fit, newdata = data_test, probability=1)
  yh <- attr(pred, "probabilities")
  positive = 1
  negative = 0
  auc_svm[[1]] <- AUC(yh[,"1"], data_test$ytest, negative, positive)
  brier_svm[[1]] <- Brier(yh[,"1"], data_test$ytest, negative, positive)
  response = as.factor(as.numeric(yh[,"1"] > 0.5))
  accuracy_svm[[1]] <- ACC(data_test$ytest, response)
  
  
  
  
#### 8
  # Combine results in list:
  
  res <- list(runningtime = runningtime,
              
              auc_svm = auc_svm,
              brier_svm = brier_svm,
              accuracy_svm = accuracy_svm,
              
              auc_rf = auc_rf,
              brier_rf = brier_rf,
              accuracy_rf = accuracy_rf,
             
              settingind=iter)
  
#### 9   
  # Save results in folder:
    save(res, file=paste("./VarselCompStudy/Additional_file/Results/res", iter, ".Rda", sep=""))

    
    #return(res)
    
}






# Function to generate the splittings for cross-validation:

# Input parameters:

# n    - number of observations in the data set
# ncv  - number of folds to use

makeCVdiv <- function(n, ncv) {
  
  nperfold <- ceiling(n/ncv)
  partition <- sample(rep(1:ncv, nperfold), n)
  
  partition
  
}

# Load data set:

#load("C:/Users/cellar/Desktop/R_try/data/LUAD.RData")
#nselvars = 10


#### 1
# t-test

t_index_inner <- function(Xtrain, ytrain, nselvars){
  t_Xtrain <- sapply(Xtrain, function(x) { t.test(x ~  ytrain)$p.value })#no_lable
  index <-  order(t_Xtrain, decreasing = 1)[1:nselvars]
}

t_index <- function(Xtrain, ytrain, block, nselvars, selectseparately, alwaysincludclin){
  
  if(selectseparately == "0" && alwaysincludclin == "0"){         # 1. naive + non-clinic
    clinindex <- block[[1]]
    Xtrain_noclin <- Xtrain[,-clinindex]
    index <- t_index_inner(Xtrain_noclin, ytrain, nselvars)
    index <- sapply(index, function(x) { x+length(clinindex) })
    return(index)
  }else if (selectseparately == "0" && alwaysincludclin == "1"){    # 2. naive + clinic
    clinindex <- block[[1]]
    Xtrain_noclin <- Xtrain[,-clinindex]
    index <- t_index_inner(Xtrain_noclin, ytrain, nselvars)
    index <- sapply(index, function(x) { x+length(clinindex) })
    index <- c(clinindex,index)
    return(index)
  } else if (selectseparately == "1" && alwaysincludclin == "0"){   # 3. group + non-clinic
    clinindex <- block[[1]]
    sumlength <- ncol(Xtrain[,-clinindex])
    #mirna
    Xtrain_mirna <- Xtrain[,block[[2]]]
    index_mirna <- t_index_inner(Xtrain_mirna, ytrain, ceiling(nselvars/sumlength*length(Xtrain_mirna)))
    mirna <- Xtrain_mirna[index_mirna]
    #mutation
    Xtrain_mutation <- Xtrain[,block[[3]]]
    index_mutation <- t_index_inner(Xtrain_mutation, ytrain, round(nselvars/sumlength*length(Xtrain_mutation)))
    mutation <- Xtrain_mutation[index_mutation]
    #cnv
    Xtrain_cnv <- Xtrain[,block[[4]]]
    index_cnv <- t_index_inner(Xtrain_cnv, ytrain, round(nselvars/sumlength*length(Xtrain_cnv)))
    cnv <- Xtrain_cnv[index_cnv]
    #rna
    Xtrain_rna <- Xtrain[,block[[5]]]
    index_rna <- t_index_inner(Xtrain_rna, ytrain, round(nselvars/sumlength*length(Xtrain_rna)))
    rna <- Xtrain_rna[index_rna]
    # combime
    data_fs <- cbind.data.frame(mirna, mutation, cnv, rna)
    index <- which(colnames(Xtrain) %in% colnames(data_fs))
    #index <- sapply(index, function(x) { x+length(clinindex) })
    return(index)
    
    
  } else if (selectseparately == "1" && alwaysincludclin == "1"){    # 4. group + clinic
    clinindex <- block[[1]]
    sumlength <- ncol(Xtrain[,-clinindex])
    #mirna
    Xtrain_mirna <- Xtrain[,block[[2]]]
    index_mirna <- t_index_inner(Xtrain_mirna, ytrain, ceiling(nselvars/sumlength*length(Xtrain_mirna)))
    mirna <- Xtrain_mirna[index_mirna]
    #mutation
    Xtrain_mutation <- Xtrain[,block[[3]]]
    index_mutation <- t_index_inner(Xtrain_mutation, ytrain, round(nselvars/sumlength*length(Xtrain_mutation)))
    mutation <- Xtrain_mutation[index_mutation]
    #cnv
    Xtrain_cnv <- Xtrain[,block[[4]]]
    index_cnv <- t_index_inner(Xtrain_cnv, ytrain, round(nselvars/sumlength*length(Xtrain_cnv)))
    cnv <- Xtrain_cnv[index_cnv]
    #rna
    Xtrain_rna <- Xtrain[,block[[5]]]
    index_rna <- t_index_inner(Xtrain_rna, ytrain, round(nselvars/sumlength*length(Xtrain_rna)))
    rna <- Xtrain_rna[index_rna]
    # combime
    data_fs <- cbind.data.frame(mirna, mutation, cnv,  rna)
    index <- which(colnames(Xtrain) %in% colnames(data_fs))
    index <- c(clinindex,index)
    return(index)
  }
}

#### 2
#chi-squared

chi_index_inner <- function(Xtrain, ytrain, nselvars){
  chi_Xtrain <- sapply(Xtrain, function(x) { chisq.test(table(x, ytrain))$p.value })#no_lable
  index <- order(chi_Xtrain, decreasing = 1)[1:nselvars]
  #chi_Xtrain_sub <- Xtrain[ order(chi_Xtrain, decreasing = 1)[1:nselvars]] #n is here
  #index <- which(colnames(Xtrain) %in% colnames(chi_Xtrain_sub))
}

chi_index <- function(Xtrain, ytrain, block, nselvars, selectseparately, alwaysincludclin){
  
  if(selectseparately == "0" && alwaysincludclin == "0"){         # 1. naive + non-clinic
    clinindex <- block[[1]]
    Xtrain_noclin <- Xtrain[,-clinindex]
    index <- chi_index_inner(Xtrain_noclin, ytrain, nselvars)
    index <- sapply(index, function(x) { x+length(clinindex) })
    return(index)
  }else if (selectseparately == "0" && alwaysincludclin == "1"){    # 2. naive + clinic
    clinindex <- block[[1]]
    Xtrain_noclin <- Xtrain[,-clinindex]
    index <- chi_index_inner(Xtrain_noclin, ytrain, nselvars)
    index <- sapply(index, function(x) { x+length(clinindex) })
    index <- c(clinindex,index)
    return(index)
  } else if (selectseparately == "1" && alwaysincludclin == "0"){   # 3. group + non-clinic
    clinindex <- block[[1]]
    sumlength <- ncol(Xtrain[,-clinindex])
    #mirna
    Xtrain_mirna <- Xtrain[,block[[2]]]
    index_mirna <- chi_index_inner(Xtrain_mirna, ytrain, ceiling(nselvars/sumlength*length(Xtrain_mirna)))
    mirna <- Xtrain_mirna[index_mirna]
    #mutation
    Xtrain_mutation <- Xtrain[,block[[3]]]
    index_mutation <- chi_index_inner(Xtrain_mutation, ytrain, round(nselvars/sumlength*length(Xtrain_mutation)))
    mutation <- Xtrain_mutation[index_mutation]
    #cnv
    Xtrain_cnv <- Xtrain[,block[[4]]]
    index_cnv <- chi_index_inner(Xtrain_cnv, ytrain, round(nselvars/sumlength*length(Xtrain_cnv)))
    cnv <- Xtrain_cnv[index_cnv]
    #rna
    Xtrain_rna <- Xtrain[,block[[5]]]
    index_rna <- chi_index_inner(Xtrain_rna, ytrain, round(nselvars/sumlength*length(Xtrain_rna)))
    rna <- Xtrain_rna[index_rna]
    # combime
    data_fs <- cbind.data.frame(mirna, mutation, cnv, rna)
    index <- which(colnames(Xtrain) %in% colnames(data_fs))
    #index <- sapply(index, function(x) { x+length(clinindex) })
    return(index)
    
    
  } else if (selectseparately == "1" && alwaysincludclin == "1"){    # 4. group + clinic
    clinindex <- block[[1]]
    sumlength <- ncol(Xtrain[,-clinindex])
    #mirna
    Xtrain_mirna <- Xtrain[,block[[2]]]
    index_mirna <- chi_index_inner(Xtrain_mirna, ytrain, ceiling(nselvars/sumlength*length(Xtrain_mirna)))
    mirna <- Xtrain_mirna[index_mirna]
    #mutation
    Xtrain_mutation <- Xtrain[,block[[3]]]
    index_mutation <- chi_index_inner(Xtrain_mutation, ytrain, round(nselvars/sumlength*length(Xtrain_mutation)))
    mutation <- Xtrain_mutation[index_mutation]
    #cnv
    Xtrain_cnv <- Xtrain[,block[[4]]]
    index_cnv <- chi_index_inner(Xtrain_cnv, ytrain, round(nselvars/sumlength*length(Xtrain_cnv)))
    cnv <- Xtrain_cnv[index_cnv]
    #rna
    Xtrain_rna <- Xtrain[,block[[5]]]
    index_rna <- chi_index_inner(Xtrain_rna, ytrain, round(nselvars/sumlength*length(Xtrain_rna)))
    rna <- Xtrain_rna[index_rna]
    # combime
    data_fs <- cbind.data.frame(mirna, mutation, cnv, rna)
    index <- which(colnames(Xtrain) %in% colnames(data_fs))
    index <- c(clinindex,index)
    return(index)
  }
}

#### 3
#infor-gain

infor_index_inner <- function(Xtrain, ytrain, nselvars){
  
  #t-test was used to get top 10% features
  t_Xtrain <- sapply(Xtrain, function(x) { t.test(x ~  ytrain)$p.value })#no_lable
  t_Xtrain <- Xtrain[ order(t_Xtrain, decreasing = 1)[1:round(ncol(Xtrain )/10)]] #n is here
  
  #information.gain
  weights <- information.gain(ytrain~., t_Xtrain)
  subset <- cutoff.k(weights, nselvars)  #n is here, subset is colname
  #r_data<- t_data[,which(colnames(t_data) %in% subset)]#return
  index <- which(colnames(Xtrain) %in% subset)
}


infor_index <- function(Xtrain, ytrain, block, nselvars, selectseparately, alwaysincludclin){
  
  if(selectseparately == "0" && alwaysincludclin == "0"){         # 1. naive + non-clinic
    clinindex <- block[[1]]
    Xtrain_noclin <- Xtrain[,-clinindex]
    index <- infor_index_inner(Xtrain_noclin, ytrain, nselvars)
    index <- sapply(index, function(x) { x+length(clinindex) })
    return(index)
  }else if (selectseparately == "0" && alwaysincludclin == "1"){    # 2. naive + clinic
    clinindex <- block[[1]]
    Xtrain_noclin <- Xtrain[,-clinindex]
    index <- infor_index_inner(Xtrain_noclin, ytrain, nselvars)
    index <- sapply(index, function(x) { x+length(clinindex) })
    index <- c(clinindex,index)
    return(index)
  } else if (selectseparately == "1" && alwaysincludclin == "0"){   # 3. group + non-clinic
    clinindex <- block[[1]]
    sumlength <- ncol(Xtrain[,-clinindex])
    #mirna
    Xtrain_mirna <- Xtrain[,block[[2]]]
    #information.gain
    weights <- information.gain(ytrain~., Xtrain_mirna)
    subset <- cutoff.k(weights, ceiling(nselvars/sumlength*length(Xtrain_mirna)))  #n is here, subset is colname
    #r_data<- t_data[,which(colnames(t_data) %in% subset)]#return
    index_mirna <- which(colnames(Xtrain_mirna) %in% subset)
    mirna <- Xtrain_mirna[index_mirna]
    #mutation
    Xtrain_mutation <- Xtrain[,block[[3]]]
    index_mutation <- infor_index_inner(Xtrain_mutation, ytrain, round(nselvars/sumlength*length(Xtrain_mutation)))
    mutation <- Xtrain_mutation[index_mutation]
    #cnv
    Xtrain_cnv <- Xtrain[,block[[4]]]
    index_cnv <- infor_index_inner(Xtrain_cnv, ytrain, round(nselvars/sumlength*length(Xtrain_cnv)))
    cnv <- Xtrain_cnv[index_cnv]
    #rna
    Xtrain_rna <- Xtrain[,block[[5]]]
    index_rna <- infor_index_inner(Xtrain_rna, ytrain, round(nselvars/sumlength*length(Xtrain_rna)))
    rna <- Xtrain_rna[index_rna]
    # combime
    data_fs <- cbind.data.frame(mirna, mutation, cnv, rna)
    index <- which(colnames(Xtrain) %in% colnames(data_fs))
    #index <- sapply(index, function(x) { x+length(clinindex) })
    return(index)
    
  } else if (selectseparately == "1" && alwaysincludclin == "1"){    # 4. group + clinic
    clinindex <- block[[1]]
    sumlength <- ncol(Xtrain[,-clinindex])
    #mirna
    Xtrain_mirna <- Xtrain[,block[[2]]]
    #information.gain
    weights <- information.gain(ytrain~., Xtrain_mirna)
    subset <- cutoff.k(weights, ceiling(nselvars/sumlength*length(Xtrain_mirna)))  #n is here, subset is colname
    #r_data<- t_data[,which(colnames(t_data) %in% subset)]#return
    index_mirna <- which(colnames(Xtrain_mirna) %in% subset)
    mirna <- Xtrain_mirna[index_mirna]
    #mutation
    Xtrain_mutation <- Xtrain[,block[[3]]]
    index_mutation <- infor_index_inner(Xtrain_mutation, ytrain, round(nselvars/sumlength*length(Xtrain_mutation)))
    mutation <- Xtrain_mutation[index_mutation]
    #cnv
    Xtrain_cnv <- Xtrain[,block[[4]]]
    index_cnv <- infor_index_inner(Xtrain_cnv, ytrain, round(nselvars/sumlength*length(Xtrain_cnv)))
    cnv <- Xtrain_cnv[index_cnv]
    #rna
    Xtrain_rna <- Xtrain[,block[[5]]]
    index_rna <- infor_index_inner(Xtrain_rna, ytrain, round(nselvars/sumlength*length(Xtrain_rna)))
    rna <- Xtrain_rna[index_rna]
    # combime
    data_fs <- cbind.data.frame(mirna, mutation, cnv, rna)
    index <- which(colnames(Xtrain) %in% colnames(data_fs))
    index <- c(clinindex,index)
    return(index)
  }
}

#### 4
#relief

relief_index_inner <- function(Xtrain, ytrain, nselvars){
  
  #t-test was used to get top 10% features
  t_Xtrain <- sapply(Xtrain, function(x) { t.test(x ~  ytrain)$p.value })#no_lable
  t_Xtrain <- Xtrain[ order(t_Xtrain, decreasing = 1)[1:round(ncol(Xtrain )/10)]] #n is here
  
  #relief
  data <- cbind(ytrain, t_Xtrain)
  weights<- relief(ytrain~., data, neighbours.count = 5, sample.size = 10)
  #weights <- information.gain(ytrain~., t_Xtrain)
  subset <- cutoff.k(weights, nselvars)  #n is here, subset is colname
  index <- which(colnames(Xtrain) %in% subset)
  
}

relief_index <- function(Xtrain, ytrain, block, nselvars, selectseparately, alwaysincludclin){
  
  if(selectseparately == "0" && alwaysincludclin == "0"){         # 1. naive + non-clinic
    clinindex <- block[[1]]
    Xtrain_noclin <- Xtrain[,-clinindex]
    index <- relief_index_inner(Xtrain_noclin, ytrain, nselvars)
    index <- sapply(index, function(x) { x+length(clinindex) })
    return(index)
  }else if (selectseparately == "0" && alwaysincludclin == "1"){    # 2. naive + clinic
    clinindex <- block[[1]]
    Xtrain_noclin <- Xtrain[,-clinindex]
    index <- relief_index_inner(Xtrain_noclin, ytrain, nselvars)
    index <- sapply(index, function(x) { x+length(clinindex) })
    index <- c(clinindex,index)
    return(index)
  } else if (selectseparately == "1" && alwaysincludclin == "0"){   # 3. group + non-clinic
    clinindex <- block[[1]]
    sumlength <- ncol(Xtrain[,-clinindex])
    #mirna
    Xtrain_mirna <- Xtrain[,block[[2]]]
    #relief
    data <- cbind(ytrain,  Xtrain_mirna)
    weights<- relief(ytrain~., data, neighbours.count = 5, sample.size = 10)
    #weights<- relief(ytrain~., Xtrain_mirna, neighbours.count = 5, sample.size = 10)
    subset <- cutoff.k(weights, ceiling(nselvars/sumlength*length(Xtrain_mirna)))  #n is here, subset is colname
    index_mirna <- which(colnames(Xtrain_mirna) %in% subset)
    mirna <- Xtrain_mirna[index_mirna]
    #mutation
    Xtrain_mutation <- Xtrain[,block[[3]]]
    index_mutation <- relief_index_inner(Xtrain_mutation, ytrain, round(nselvars/sumlength*length(Xtrain_mutation)))
    mutation <- Xtrain_mutation[index_mutation]
    #cnv
    Xtrain_cnv <- Xtrain[,block[[4]]]
    index_cnv <- relief_index_inner(Xtrain_cnv, ytrain, round(nselvars/sumlength*length(Xtrain_cnv)))
    cnv <- Xtrain_cnv[index_cnv]
    #rna
    Xtrain_rna <- Xtrain[,block[[5]]]
    index_rna <- relief_index_inner(Xtrain_rna, ytrain, round(nselvars/sumlength*length(Xtrain_rna)))
    rna <- Xtrain_rna[index_rna]
    # combime
    data_fs <- cbind.data.frame(mirna, mutation, cnv, rna)
    index <- which(colnames(Xtrain) %in% colnames(data_fs))
    #index <- sapply(index, function(x) { x+length(clinindex) })
    return(index)
    
  } else if (selectseparately == "1" && alwaysincludclin == "1"){    # 4. group + clinic
    clinindex <- block[[1]]
    sumlength <- ncol(Xtrain[,-clinindex])
    #mirna
    Xtrain_mirna <- Xtrain[,block[[2]]]
    #relief
    data <- cbind(ytrain,  Xtrain_mirna)
    weights<- relief(ytrain~., data, neighbours.count = 5, sample.size = 10)
    #weights<- relief(ytrain~., Xtrain_mirna, neighbours.count = 5, sample.size = 10)
    subset <- cutoff.k(weights, ceiling(nselvars/sumlength*length(Xtrain_mirna)))  #n is here, subset is colname
    index_mirna <- which(colnames(Xtrain_mirna) %in% subset)
    mirna <- Xtrain_mirna[index_mirna]
    #mutation
    Xtrain_mutation <- Xtrain[,block[[3]]]
    index_mutation <- relief_index_inner(Xtrain_mutation, ytrain, round(nselvars/sumlength*length(Xtrain_mutation)))
    mutation <- Xtrain_mutation[index_mutation]
    #cnv
    Xtrain_cnv <- Xtrain[,block[[4]]]
    index_cnv <- relief_index_inner(Xtrain_cnv, ytrain, round(nselvars/sumlength*length(Xtrain_cnv)))
    cnv <- Xtrain_cnv[index_cnv]
    #rna
    Xtrain_rna <- Xtrain[,block[[5]]]
    index_rna <- relief_index_inner(Xtrain_rna, ytrain, round(nselvars/sumlength*length(Xtrain_rna)))
    rna <- Xtrain_rna[index_rna]
    # combime
    data_fs <- cbind.data.frame(mirna, mutation, cnv,  rna)
    index <- which(colnames(Xtrain) %in% colnames(data_fs))
    index <- c(clinindex,index)
    return(index)
  }
}



#### 5
#mrmr

mrmr_index_inner <- function(Xtrain, ytrain, nselvars){
  
  #t-test was used to get top 10% features
  t_Xtrain <- sapply(Xtrain, function(x) { t.test(x ~  ytrain)$p.value })#no_lable
  t_Xtrain <- Xtrain[ order(t_Xtrain, decreasing = 1)[1:round(ncol(Xtrain )/10)]] #n is here
  
  #mrmr
  mrmr_feature<-t_Xtrain
  mrmr_feature$y <-ytrain
  target_indices = which(names(mrmr_feature)=='y')
  for (m in which(sapply(mrmr_feature, class)!="numeric")){
    mrmr_feature[,m]=as.numeric(mrmr_feature[,m])
  }
  Data <- mRMR.data(data = data.frame(mrmr_feature))
  mrmr=mRMR.ensemble(data = Data, target_indices = target_indices, 
                     feature_count = nselvars, solution_count = 1) #n is here
  index=mrmr@filters[[as.character(mrmr@target_indices)]]#no lable, index is a list.
  mrmr_data <- mrmr_feature[as.numeric(index)]
  index <- which(colnames(Xtrain) %in% colnames(mrmr_data))
}

mrmr_index <- function(Xtrain, ytrain, block, nselvars, selectseparately, alwaysincludclin){
  
  if(selectseparately == "0" && alwaysincludclin == "0"){         # 1. naive + non-clinic
    clinindex <- block[[1]]
    Xtrain_noclin <- Xtrain[,-clinindex]
    index <- mrmr_index_inner(Xtrain_noclin, ytrain, nselvars)
    index <- sapply(index, function(x) { x+length(clinindex) })
    return(index)
  }else if (selectseparately == "0" && alwaysincludclin == "1"){    # 2. naive + clinic
    clinindex <- block[[1]]
    Xtrain_noclin <- Xtrain[,-clinindex]
    index <- mrmr_index_inner(Xtrain_noclin, ytrain, nselvars)
    index <- sapply(index, function(x) { x+length(clinindex) })
    index <- c(clinindex,index)
    return(index)
  } else if (selectseparately == "1" && alwaysincludclin == "0"){   # 3. group + non-clinic
    clinindex <- block[[1]]
    sumlength <- ncol(Xtrain[,-clinindex])
    #mirna
    Xtrain_mirna <- Xtrain[,block[[2]]]
    #mrmr
    mrmr_feature<-Xtrain_mirna
    mrmr_feature$y <-ytrain
    target_indices = which(names(mrmr_feature)=='y')
    for (m in which(sapply(mrmr_feature, class)!="numeric")){
      mrmr_feature[,m]=as.numeric(mrmr_feature[,m])
    }
    Data <- mRMR.data(data = data.frame(mrmr_feature))
    mrmr=mRMR.ensemble(data = Data, target_indices = target_indices, 
                       feature_count = ceiling(nselvars/sumlength*length(Xtrain_mirna)), solution_count = 1) #n is here
    index=mrmr@filters[[as.character(mrmr@target_indices)]]#no lable
    mrmr_data <- mrmr_feature[as.numeric(index)]
    index_mirna <- which(colnames(Xtrain_mirna) %in% colnames(mrmr_data))
    mirna <- Xtrain_mirna[index_mirna]
    #mutation
    Xtrain_mutation <- Xtrain[,block[[3]]]
    index_mutation <- mrmr_index_inner(Xtrain_mutation, ytrain, round(nselvars/sumlength*length(Xtrain_mutation)))
    mutation <- Xtrain_mutation[index_mutation]
    #cnv
    Xtrain_cnv <- Xtrain[,block[[4]]]
    index_cnv <- mrmr_index_inner(Xtrain_cnv, ytrain, round(nselvars/sumlength*length(Xtrain_cnv)))
    cnv <- Xtrain_cnv[index_cnv]
    #rna
    Xtrain_rna <- Xtrain[,block[[5]]]
    index_rna <- mrmr_index_inner(Xtrain_rna, ytrain, round(nselvars/sumlength*length(Xtrain_rna)))
    rna <- Xtrain_rna[index_rna]
    # combime
    data_fs <- cbind.data.frame(mirna, mutation, cnv, rna)
    index <- which(colnames(Xtrain) %in% colnames(data_fs))
    #index <- sapply(index, function(x) { x+length(clinindex) })
    return(index)
    
  } else if (selectseparately == "1" && alwaysincludclin == "1"){    # 4. group + clinic
    clinindex <- block[[1]]
    sumlength <- ncol(Xtrain[,-clinindex])
    #mirna
    Xtrain_mirna <- Xtrain[,block[[2]]]
    #mrmr
    mrmr_feature<-Xtrain_mirna
    mrmr_feature$y <-ytrain
    target_indices = which(names(mrmr_feature)=='y')
    for (m in which(sapply(mrmr_feature, class)!="numeric")){
      mrmr_feature[,m]=as.numeric(mrmr_feature[,m])
    }
    Data <- mRMR.data(data = data.frame(mrmr_feature))
    mrmr=mRMR.ensemble(data = Data, target_indices = target_indices, 
                       feature_count = ceiling(nselvars/sumlength*length(Xtrain_mirna)), solution_count = 1) #n is here
    index=mrmr@filters[[as.character(mrmr@target_indices)]]#no lable
    mrmr_data <- mrmr_feature[as.numeric(index)]
    index_mirna <- which(colnames(Xtrain_mirna) %in% colnames(mrmr_data))
    mirna <- Xtrain_mirna[index_mirna]
    #mutation
    Xtrain_mutation <- Xtrain[,block[[3]]]
    index_mutation <- mrmr_index_inner(Xtrain_mutation, ytrain, round(nselvars/sumlength*length(Xtrain_mutation)))
    mutation <- Xtrain_mutation[index_mutation]
    #cnv
    Xtrain_cnv <- Xtrain[,block[[4]]]
    index_cnv <- mrmr_index_inner(Xtrain_cnv, ytrain, round(nselvars/sumlength*length(Xtrain_cnv)))
    cnv <- Xtrain_cnv[index_cnv]
    #rna
    Xtrain_rna <- Xtrain[,block[[5]]]
    index_rna <- mrmr_index_inner(Xtrain_rna, ytrain, round(nselvars/sumlength*length(Xtrain_rna)))
    rna <- Xtrain_rna[index_rna]
    # combime
    data_fs <- cbind.data.frame(mirna, mutation, cnv, rna)
    index <- which(colnames(Xtrain) %in% colnames(data_fs))
    index <- c(clinindex,index)
    return(index)
  }
}
#### 6
#rf

rf_index_inner <- function(Xtrain, ytrain, nselvars){
  
  #t-test was used to get top 10% features
  t_Xtrain <- sapply(Xtrain, function(x) { t.test(x ~  ytrain)$p.value })#no_lable
  t_Xtrain <- Xtrain[ order(t_Xtrain, decreasing = 1)[1:round(ncol(Xtrain )/10)]] #n is here
  
  #RF
  data <- cbind(ytrain, t_Xtrain)
  RF <- ranger(ytrain ~ ., data = data, importance = "impurity")
  ranker <- RF$variable.importance[order(-RF$variable.importance)]
  ranker_s <- data.frame(ranker[1:nselvars])
  index <- which(colnames(Xtrain) %in% rownames(ranker_s))
  
}


rf_index <- function(Xtrain, ytrain, block, nselvars, selectseparately, alwaysincludclin){
  
  if(selectseparately == "0" && alwaysincludclin == "0"){         # 1. naive + non-clinic
    clinindex <- block[[1]]
    Xtrain_noclin <- Xtrain[,-clinindex]
    index <- rf_index_inner(Xtrain_noclin, ytrain, nselvars)
    index <- sapply(index, function(x) { x+length(clinindex) })
    return(index)
  }else if (selectseparately == "0" && alwaysincludclin == "1"){    # 2. naive + clinic
    clinindex <- block[[1]]
    Xtrain_noclin <- Xtrain[,-clinindex]
    index <- rf_index_inner(Xtrain_noclin, ytrain, nselvars)
    index <- sapply(index, function(x) { x+length(clinindex) })
    index <- c(clinindex,index)
    return(index)
  } else if (selectseparately == "1" && alwaysincludclin == "0"){   # 3. group + non-clinic
    clinindex <- block[[1]]
    sumlength <- ncol(Xtrain[,-clinindex])
    #mirna
    Xtrain_mirna <- Xtrain[,block[[2]]]
    #RF
    data <- cbind(ytrain, Xtrain_mirna)
    RF <- ranger(ytrain ~ ., data = data, importance = "impurity")
    ranker <- RF$variable.importance[order(-RF$variable.importance)]
    ranker_s <- data.frame(ranker[1:ceiling(nselvars/sumlength*length(Xtrain_mirna))])
    index_mirna <- which(colnames(Xtrain_mirna) %in% rownames(ranker_s))
    mirna <- Xtrain_mirna[index_mirna]
    #mutation
    Xtrain_mutation <- Xtrain[,block[[3]]]
    index_mutation <- rf_index_inner(Xtrain_mutation, ytrain, round(nselvars/sumlength*length(Xtrain_mutation)))
    mutation <- Xtrain_mutation[index_mutation]
    #cnv
    Xtrain_cnv <- Xtrain[,block[[4]]]
    index_cnv <- rf_index_inner(Xtrain_cnv, ytrain, round(nselvars/sumlength*length(Xtrain_cnv)))
    cnv <- Xtrain_cnv[index_cnv]
    #rna
    Xtrain_rna <- Xtrain[,block[[5]]]
    index_rna <- rf_index_inner(Xtrain_rna, ytrain, round(nselvars/sumlength*length(Xtrain_rna)))
    rna <- Xtrain_rna[index_rna]
    # combime
    data_fs <- cbind.data.frame(mirna, mutation, cnv, rna)
    index <- which(colnames(Xtrain) %in% colnames(data_fs))
    #index <- sapply(index, function(x) { x+length(clinindex) })
    return(index)
    
  } else if (selectseparately == "1" && alwaysincludclin == "1"){    # 4. group + clinic
    clinindex <- block[[1]]
    sumlength <- ncol(Xtrain[,-clinindex])
    #mirna
    Xtrain_mirna <- Xtrain[,block[[2]]]
    #RF
    data <- cbind(ytrain, Xtrain_mirna)
    RF <- ranger(ytrain ~ ., data = data, importance = "impurity")
    ranker <- RF$variable.importance[order(-RF$variable.importance)]
    ranker_s <- data.frame(ranker[1:ceiling(nselvars/sumlength*length(Xtrain_mirna))])
    index_mirna <- which(colnames(Xtrain_mirna) %in% rownames(ranker_s))
    mirna <- Xtrain_mirna[index_mirna]
    #mutation
    Xtrain_mutation <- Xtrain[,block[[3]]]
    index_mutation <- rf_index_inner(Xtrain_mutation, ytrain, round(nselvars/sumlength*length(Xtrain_mutation)))
    mutation <- Xtrain_mutation[,index_mutation]
    #cnv
    Xtrain_cnv <- Xtrain[,block[[4]]]
    index_cnv <- rf_index_inner(Xtrain_cnv, ytrain, round(nselvars/sumlength*length(Xtrain_cnv)))
    cnv <- Xtrain_cnv[index_cnv]
    #rna
    Xtrain_rna <- Xtrain[,block[[5]]]
    index_rna <- rf_index_inner(Xtrain_rna, ytrain, round(nselvars/sumlength*length(Xtrain_rna)))
    rna <- Xtrain_rna[index_rna]
    # combime
    data_fs <- cbind.data.frame(mirna, mutation, cnv,  rna)
    index <- which(colnames(Xtrain) %in% colnames(data_fs))
    index <- c(clinindex,index)
    return(index)
  }
}

