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
  library(glmnet)
  library(caret)
  #library(rJava)
  #library(mRMRe)
  library(e1071)
  library(ranger)
  library(measures)#brier.score
  
#### 1  
  # Initiate lists in which the results
  # will be stored:
  
  runningtime <- list()
  FS_N <- list()
  
  auc_svm <- list()
  brier_svm <- list()
  accuracy_svm <- list()
  
  auc_rf <- list()
  brier_rf <- list()
  accuracy_rf <- list()
  
  
#### 2  
  # Obtain information for the iter-th setting:
  #iter = 1
  dat <- scenariogrid_1$dat[iter]
  seed <- scenariogrid_1$seed[iter]
  featuremethod <- scenariogrid_1$featuremethod[iter]
 # nselvars <- scenariogrid_1$nselvars[iter]
  selectseparately <- scenariogrid_1$selectseparately[iter]
  alwaysincludclin <- scenariogrid_1$alwaysincludclin[iter]
  
  cvind <- scenariogrid_1$cvind[iter]
  cvfoldind <- scenariogrid_1$cvfoldind[iter]
  
#### 3  
  # Load data set:
  
  load(paste("./VarselCompStudy/Additional_file/Data/", dat, sep=""))
  #prepare data, for testing
  #cnvdata <- cnvdata[,sample(ncol(cnvdata),size=1003)]
  #mutationdata <- mutationdata[,sample(ncol(mutationdata),size=1005)]
 # rnadata <- rnadata[,sample(ncol(rnadata),size=1007)]
  
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
  
  if(featuremethod=="CFS") {
    
    index <- CFS_index(Xtrain, ytrain, block, selectseparately, alwaysincludclin)
    
  }else if(featuremethod=="ga"){
    
    index <- ga_index(Xtrain, ytrain, block,  selectseparately, alwaysincludclin)
    
  }else if(featuremethod=="rfe"){
  
    index <- rfe_index(Xtrain, ytrain, block, selectseparately, alwaysincludclin)
    
  }else if(featuremethod=="lasso"){
    
    index <- lasso_index(Xtrain, ytrain, block,  selectseparately, alwaysincludclin)
    
  }
  
  timeend <-Sys.time()
  runningtime[[1]] <-timeend-timestart
  FS_N[[1]] <- length(index)
  
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
  positive = 0
  negative = 1
  auc_rf[[1]] <- AUC(yh[,1], data_test$ytest, negative, positive)
  brier_rf[[1]] <- Brier(yh[,1], data_test$ytest, negative, positive)
  response = as.factor(as.numeric(yh[,2] > 0.5))
  accuracy_rf[[1]] <- ACC(data_test$ytest, response)
    
  #SVM
  fit <- svm(ytrain ~ ., data = data_train,
             type = 'C-classification', kernel = 'radial',
             #gamma = g,cost=c,
             probability=1)
  pred <- predict(fit, newdata = data_test, probability=1)
  yh <- attr(pred, "probabilities")
  positive = 0
  negative = 1
  auc_svm[[1]] <- AUC(yh[,1], data_test$ytest, negative, positive)
  brier_svm[[1]] <- Brier(yh[,1], data_test$ytest, negative, positive)
  response = as.factor(as.numeric(yh[,2] > 0.5))
  accuracy_svm[[1]] <- ACC(data_test$ytest, response)
  
  
  
  
#### 8
  # Combine results in list:
  
  res <- list(runningtime = runningtime,
                FS_N =  FS_N,
              
              auc_svm = auc_svm,
              brier_svm = brier_svm,
              accuracy_svm = accuracy_svm,
              
              auc_rf = auc_rf,
              brier_rf = brier_rf,
              accuracy_rf = accuracy_rf,
             
              settingind=iter)
  
#### 9   
  # Save results in folder:
    save(res, file=paste("./VarselCompStudy/Additional_file/Results_1/res", iter, ".Rda", sep=""))

    
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
#CFS

#infor-gain

CFS_index_inner <- function(Xtrain, ytrain){
  
  #t-test was used to get top 10% features
  t_Xtrain <- sapply(Xtrain, function(x) { t.test(x ~  ytrain)$p.value })#no_lable
  t_Xtrain <- Xtrain[ order(t_Xtrain, decreasing = 1)[1:round(ncol(Xtrain )/10)]] #n is here
  
  
  #CFS
  data <- cbind(ytrain, t_Xtrain)
  subset <- cfs(ytrain~., data)
  index <- which(colnames(Xtrain) %in% subset)
  
}
CFS_index <- function(Xtrain, ytrain, block, selectseparately, alwaysincludclin){
  
  if(selectseparately == "0" && alwaysincludclin == "0"){         # 1. naive + non-clinic
    clinindex <- block[[1]]
    Xtrain_noclin <- Xtrain[,-clinindex]
    index <- CFS_index_inner(Xtrain_noclin, ytrain)
    index <- sapply(index, function(x) { x+length(clinindex) })
    return(index)
  }else if (selectseparately == "0" && alwaysincludclin == "1"){    # 2. naive + clinic
    clinindex <- block[[1]]
    Xtrain_noclin <- Xtrain[,-clinindex]
    index <- CFS_index_inner(Xtrain_noclin, ytrain)
    index <- sapply(index, function(x) { x+length(clinindex) })
    index <- c(clinindex,index)
    return(index)
  } else if (selectseparately == "1" && alwaysincludclin == "0"){   # 3. group + non-clinic
    clinindex <- block[[1]]
    sumlength <- ncol(Xtrain[,-clinindex])
    #mirna
    Xtrain_mirna <- Xtrain[,block[[2]]]
    data <- cbind(ytrain, Xtrain_mirna)
    subset <- cfs(ytrain~., data)
    index_mirna <- which(colnames(Xtrain_mirna) %in% subset)
    mirna <- Xtrain_mirna[index_mirna]
    
    #mutation
    Xtrain_mutation <- Xtrain[,block[[3]]]
    index_mutation <- CFS_index_inner(Xtrain_mutation, ytrain)
    mutation <- Xtrain_mutation[index_mutation]
    #cnv
    Xtrain_cnv <- Xtrain[,block[[4]]]
    index_cnv <- CFS_index_inner(Xtrain_cnv, ytrain)
    cnv <- Xtrain_cnv[index_cnv]
    #rna
    Xtrain_rna <- Xtrain[,block[[5]]]
    index_rna <- CFS_index_inner(Xtrain_rna, ytrain)
    rna <- Xtrain_rna[index_rna]
    # combine
    data_fs <- cbind.data.frame(mirna, mutation, cnv, rna)
    index <- which(colnames(Xtrain) %in% colnames(data_fs))
    #index <- sapply(index, function(x) { x+length(clinindex) })
    return(index)
    
    
  } else if (selectseparately == "1" && alwaysincludclin == "1"){    # 4. group + clinic
    clinindex <- block[[1]]
    sumlength <- ncol(Xtrain[,-clinindex])
    #mirna
    Xtrain_mirna <- Xtrain[,block[[2]]]
    data <- cbind(ytrain, Xtrain_mirna)
    subset <- cfs(ytrain~., data)
    index_mirna <- which(colnames(Xtrain_mirna) %in% subset)
    mirna <- Xtrain_mirna[index_mirna]
    #mutation
    Xtrain_mutation <- Xtrain[,block[[3]]]
    index_mutation <- CFS_index_inner(Xtrain_mutation, ytrain)
    mutation <- Xtrain_mutation[index_mutation]
    #cnv
    Xtrain_cnv <- Xtrain[,block[[4]]]
    index_cnv <- CFS_index_inner(Xtrain_cnv, ytrain)
    cnv <- Xtrain_cnv[index_cnv]
    #rna
    Xtrain_rna <- Xtrain[,block[[5]]]
    index_rna <- CFS_index_inner(Xtrain_rna, ytrain)
    rna <- Xtrain_rna[index_rna]
    # combime
    data_fs <- cbind.data.frame(mirna, mutation, cnv, rna)
    index <- which(colnames(Xtrain) %in% colnames(data_fs))
    index <- c(clinindex,index)
    return(index)
  }
}

#### 2
#GA

ga_index_inner <- function(Xtrain, ytrain){
  
  #t-test was used to get top 10% features
  t_Xtrain <- sapply(Xtrain, function(x) { t.test(x ~  ytrain)$p.value })#no_lable
  t_Xtrain <- Xtrain[ order(t_Xtrain, decreasing = 1)[1:round(ncol(Xtrain )/10)]] #n is here
  
  
  #GA
  ga_ctrl <- gafsControl(functions = rfGA,  # another option is `caretGA`.
                         method = "repeatedcv", # repeated cv
                         repeats = 1, # number of repeats
                         number = 5,
                         allowParallel = FALSE) # number of folds
  ga_obj <- gafs(x=t_Xtrain, 
                 y=ytrain, 
                 iters = 50,   # normally much higher (100+)
                 gafsControl = ga_ctrl)
  index <- which(colnames(Xtrain) %in% ga_obj$optVariables)
  
}

ga_index <- function(Xtrain, ytrain, block, selectseparately, alwaysincludclin){
  
  if(selectseparately == "0" && alwaysincludclin == "0"){         # 1. naive + non-clinic
    clinindex <- block[[1]]
    Xtrain_noclin <- Xtrain[,-clinindex]
    index <- ga_index_inner(Xtrain_noclin, ytrain)
    index <- sapply(index, function(x) { x+length(clinindex) })
    return(index)
  }else if (selectseparately == "0" && alwaysincludclin == "1"){    # 2. naive + clinic
    clinindex <- block[[1]]
    Xtrain_noclin <- Xtrain[,-clinindex]
    index <- ga_index_inner(Xtrain_noclin, ytrain)
    index <- sapply(index, function(x) { x+length(clinindex) })
    index <- c(clinindex,index)
    return(index)
  } else if (selectseparately == "1" && alwaysincludclin == "0"){   # 3. group + non-clinic
    clinindex <- block[[1]]
    sumlength <- ncol(Xtrain[,-clinindex])
    #mirna
    Xtrain_mirna <- Xtrain[,block[[2]]]
    ga_ctrl <- gafsControl(functions = rfGA,  # another option is `caretGA`.
                           method = "repeatedcv", # repeated cv
                           repeats = 1, # number of repeats
                           number = 5,
                           allowParallel = FALSE) # number of folds
    ga_obj <- gafs(x=Xtrain_mirna, 
                   y=ytrain, 
                   iters = 50,   # normally much higher (100+)
                   gafsControl = ga_ctrl)
    index_mirna <- which(colnames(Xtrain_mirna) %in% ga_obj$optVariables)
    mirna <- Xtrain_mirna[index_mirna]
    
    #mutation
    Xtrain_mutation <- Xtrain[,block[[3]]]
    index_mutation <- ga_index_inner(Xtrain_mutation, ytrain)
    mutation <- Xtrain_mutation[index_mutation]
    #cnv
    Xtrain_cnv <- Xtrain[,block[[4]]]
    index_cnv <- ga_index_inner(Xtrain_cnv, ytrain)
    cnv <- Xtrain_cnv[index_cnv]
    #rna
    Xtrain_rna <- Xtrain[,block[[5]]]
    index_rna <- ga_index_inner(Xtrain_rna, ytrain)
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
    ga_ctrl <- gafsControl(functions = rfGA,  # another option is `caretGA`.
                           method = "repeatedcv", # repeated cv
                           repeats = 1, # number of repeats
                           number = 5,
                           allowParallel = FALSE) # number of folds
    ga_obj <- gafs(x=Xtrain_mirna, 
                   y=ytrain, 
                   iters = 50,   # normally much higher (100+)
                   gafsControl = ga_ctrl)
    index_mirna <- which(colnames(Xtrain_mirna) %in% ga_obj$optVariables)
    mirna <- Xtrain_mirna[index_mirna]
    #mutation
    Xtrain_mutation <- Xtrain[,block[[3]]]
    index_mutation <- ga_index_inner(Xtrain_mutation, ytrain)
    mutation <- Xtrain_mutation[index_mutation]
    #cnv
    Xtrain_cnv <- Xtrain[,block[[4]]]
    index_cnv <- ga_index_inner(Xtrain_cnv, ytrain)
    cnv <- Xtrain_cnv[index_cnv]
    #rna
    Xtrain_rna <- Xtrain[,block[[5]]]
    index_rna <- ga_index_inner(Xtrain_rna, ytrain)
    rna <- Xtrain_rna[index_rna]
    # combime
    data_fs <- cbind.data.frame(mirna, mutation, cnv, rna)
    index <- which(colnames(Xtrain) %in% colnames(data_fs))
    index <- c(clinindex,index)
    return(index)
  }
}

#### 3
#rfe

rfe_index_inner <- function(Xtrain, ytrain){
  
  #t-test was used to get top 10% features
  t_Xtrain <- sapply(Xtrain, function(x) { t.test(x ~  ytrain)$p.value })#no_lable
  t_Xtrain <- Xtrain[ order(t_Xtrain, decreasing = 1)[1:round(ncol(Xtrain )/10)]] #n is here
  
  ######rfe
  subsets <- ncol(t_Xtrain)
  # Define the control using a random forest selection function
  control <- rfeControl(functions = rfFuncs, # random forest
                        method = "repeatedcv", # repeated cv
                        repeats = 1, # number of repeats
                        number = 5,
                        allowParallel = FALSE) # number of folds
  result_rfe1 <- rfe(x = t_Xtrain, 
                     y = ytrain, 
                     sizes = seq(from=1,to=subsets,by=100),#retained features
                     rfeControl = control)
  index <- which(colnames(Xtrain) %in% predictors(result_rfe1))
  
}

rfe_index <- function(Xtrain, ytrain, block, selectseparately, alwaysincludclin){
  
  if(selectseparately == "0" && alwaysincludclin == "0"){         # 1. naive + non-clinic
    clinindex <- block[[1]]
    Xtrain_noclin <- Xtrain[,-clinindex]
    index <- rfe_index_inner(Xtrain_noclin, ytrain)
    index <- sapply(index, function(x) { x+length(clinindex) })
    return(index)
  }else if (selectseparately == "0" && alwaysincludclin == "1"){    # 2. naive + clinic
    clinindex <- block[[1]]
    Xtrain_noclin <- Xtrain[,-clinindex]
    index <- rfe_index_inner(Xtrain_noclin, ytrain)
    index <- sapply(index, function(x) { x+length(clinindex) })
    index <- c(clinindex,index)
    return(index)
  } else if (selectseparately == "1" && alwaysincludclin == "0"){   # 3. group + non-clinic
    clinindex <- block[[1]]
    sumlength <- ncol(Xtrain[,-clinindex])
    #mirna
    Xtrain_mirna <- Xtrain[,block[[2]]]
    subsets <- ncol(Xtrain_mirna)
    control <- rfeControl(functions = rfFuncs, # random forest
                          method = "repeatedcv", # repeated cv
                          repeats = 1, # number of repeats
                          number = 5,
                          allowParallel = FALSE) # number of folds
    result_rfe1 <- rfe(x = Xtrain_mirna, 
                       y = ytrain, 
                       sizes = seq(from=1,to=subsets,by=100),#retained features
                       rfeControl = control)
    index_mirna <- which(colnames(Xtrain_mirna) %in% predictors(result_rfe1))
    mirna <- Xtrain_mirna[index_mirna]
    #mutation
    Xtrain_mutation <- Xtrain[,block[[3]]]
    index_mutation <- rfe_index_inner(Xtrain_mutation, ytrain)
    mutation <- Xtrain_mutation[index_mutation]
    #cnv
    Xtrain_cnv <- Xtrain[,block[[4]]]
    index_cnv <- rfe_index_inner(Xtrain_cnv, ytrain)
    cnv <- Xtrain_cnv[index_cnv]
    #rna
    Xtrain_rna <- Xtrain[,block[[5]]]
    index_rna <- rfe_index_inner(Xtrain_rna, ytrain)
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
    subsets <- ncol(Xtrain_mirna)
    control <- rfeControl(functions = rfFuncs, # random forest
                          method = "repeatedcv", # repeated cv
                          repeats = 1, # number of repeats
                          number = 5,
                          allowParallel = FALSE) # number of folds
    result_rfe1 <- rfe(x = Xtrain_mirna, 
                       y = ytrain, 
                       sizes = seq(from=1,to=subsets,by=100),#retained features d<-seq(from=3,to=12,by=3)
                       rfeControl = control)
    index_mirna <- which(colnames(Xtrain_mirna) %in% predictors(result_rfe1))
    mirna <- Xtrain_mirna[index_mirna]
    #mutation
    Xtrain_mutation <- Xtrain[,block[[3]]]
    index_mutation <- rfe_index_inner(Xtrain_mutation, ytrain)
    mutation <- Xtrain_mutation[index_mutation]
    #cnv
    Xtrain_cnv <- Xtrain[,block[[4]]]
    index_cnv <- rfe_index_inner(Xtrain_cnv, ytrain)
    cnv <- Xtrain_cnv[index_cnv]
    #rna
    Xtrain_rna <- Xtrain[,block[[5]]]
    index_rna <- rfe_index_inner(Xtrain_rna, ytrain)
    rna <- Xtrain_rna[index_rna]
    # combime
    data_fs <- cbind.data.frame(mirna, mutation, cnv, rna)
    index <- which(colnames(Xtrain) %in% colnames(data_fs))
    index <- c(clinindex,index)
    return(index)
  }
}

#### 4
#lasso

lasso_index_inner <- function(Xtrain, ytrain){
  
  #t-test was used to get top 10% features
  t_Xtrain <- sapply(Xtrain, function(x) { t.test(x ~  ytrain)$p.value })#no_lable
  t_Xtrain <- Xtrain[ order(t_Xtrain, decreasing = 1)[1:round(ncol(Xtrain )/10)]] #n is here
  
  #lasso
  t_Xtrain <- as.matrix(t_Xtrain)
  cv.lasso <- cv.glmnet(t_Xtrain,ytrain,type.measure = "auc",alpha=1,family="binomial",nfolds = 5)
  c<-coef(cv.lasso,s='lambda.min',exact=TRUE)
  inds<-which(c!=0)
  variables<-row.names(c)[inds]
  index <- which(colnames(Xtrain) %in% variables)

}
lasso_index <- function(Xtrain, ytrain, block, selectseparately, alwaysincludclin){
  
  if(selectseparately == "0" && alwaysincludclin == "0"){         # 1. naive + non-clinic
    clinindex <- block[[1]]
    Xtrain_noclin <- Xtrain[,-clinindex]
    index <- lasso_index_inner(Xtrain_noclin, ytrain)
    index <- sapply(index, function(x) { x+length(clinindex) })
    return(index)
  }else if (selectseparately == "0" && alwaysincludclin == "1"){    # 2. naive + clinic
    clinindex <- block[[1]]
    Xtrain_noclin <- Xtrain[,-clinindex]
    index <- lasso_index_inner(Xtrain_noclin, ytrain)
    index <- sapply(index, function(x) { x+length(clinindex) })
    index <- c(clinindex,index)
    return(index)
  } else if (selectseparately == "1" && alwaysincludclin == "0"){   # 3. group + non-clinic
    clinindex <- block[[1]]
    sumlength <- ncol(Xtrain[,-clinindex])
    #mirna
    Xtrain_mirna <- Xtrain[,block[[2]]]
    Xtrain_mirna_m <- as.matrix(Xtrain_mirna)
    cv.lasso <- cv.glmnet(Xtrain_mirna_m,ytrain,type.measure = "auc",alpha=1,family="binomial",nfolds = 5)
    c<-coef(cv.lasso ,s='lambda.min',exact=TRUE)
    inds<-which(c!=0)
    variables<-row.names(c)[inds]
    index_mirna <- which(colnames(Xtrain_mirna) %in% variables)
    mirna <- Xtrain_mirna[index_mirna]
    #mutation
    Xtrain_mutation <- Xtrain[,block[[3]]]
    index_mutation <- lasso_index_inner(Xtrain_mutation, ytrain)
    mutation <- Xtrain_mutation[index_mutation]
    #cnv
    Xtrain_cnv <- Xtrain[,block[[4]]]
    index_cnv <- lasso_index_inner(Xtrain_cnv, ytrain)
    cnv <- Xtrain_cnv[index_cnv]
    #rna
    Xtrain_rna <- Xtrain[,block[[5]]]
    index_rna <- lasso_index_inner(Xtrain_rna, ytrain)
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
    Xtrain_mirna_m <- as.matrix(Xtrain_mirna)
    cv.lasso <- cv.glmnet(Xtrain_mirna_m,ytrain,type.measure = "auc",alpha=1,family="binomial",nfolds = 5)
    c<-coef(cv.lasso ,s='lambda.min',exact=TRUE)
    inds<-which(c!=0)
    variables<-row.names(c)[inds]
    index_mirna <- which(colnames(Xtrain_mirna) %in% variables)
    mirna <- Xtrain_mirna[index_mirna]
    #mutation
    Xtrain_mutation <- Xtrain[,block[[3]]]
    index_mutation <- lasso_index_inner(Xtrain_mutation, ytrain)
    mutation <- Xtrain_mutation[index_mutation]
    #cnv
    Xtrain_cnv <- Xtrain[,block[[4]]]
    index_cnv <- lasso_index_inner(Xtrain_cnv, ytrain)
    cnv <- Xtrain_cnv[index_cnv]
    #rna
    Xtrain_rna <- Xtrain[,block[[5]]]
    index_rna <- lasso_index_inner(Xtrain_rna, ytrain)
    rna <- Xtrain_rna[index_rna]
    # combime
    data_fs <- cbind.data.frame(mirna, mutation, cnv, rna)
    index <- which(colnames(Xtrain) %in% colnames(data_fs))
    index <- c(clinindex,index)
    return(index)
  }
}

