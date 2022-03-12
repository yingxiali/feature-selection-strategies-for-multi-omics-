#### 1
# Set working directory:

setwd("~/")

#### 2
# Make table of settings:

#dat <- c("COAD.Rda", "LGG.Rda", "OV.Rda", "PAAD.Rda", "SKCM.Rda", "UCEC.Rda")
dat <- c("BLCA.RData","BRCA.RData", "COAD.RData", "ESCA.RData", "HNSC.RData", 
         "LGG.RData", "LIHC.RData", "LUAD.RData", "LUSC.RData", "PAAD.RData", 
         "PRAD.RData", "SARC.RData", "SKCM.RData", "STAD.RData", "UCEC.RData")

featuremethod <- c("t-test","chi-squared","infor","relief", "mrmr","rf") # rank method

cvind <- 1:3

cvfoldind <- 1:5

nselvars <- c(10, 100, 1000, 5000) # 0.01%, 0.1%, 1%, 5%

selectseparately <- c(1,0) # naive or group True is 1, 0 is false.

alwaysincludclin <- c(1,0) # clinic information

scenariogrid <- expand.grid(selectseparately = selectseparately, nselvars=nselvars,
                            alwaysincludclin = alwaysincludclin, cvind=cvind, dat=dat, 
                            cvfoldind=cvfoldind, featuremethod=featuremethod, stringsAsFactors = FALSE)
scenariogrid <- scenariogrid[,ncol(scenariogrid):1]


set.seed(1234)
seeds <- sample(1000:10000000, size=length(dat)*length(cvind)) #repeat 3 times

scenariogrid$seed <- rep(seeds, times=length(featuremethod)*length(nselvars)*
                         length(cvfoldind)*length(selectseparately)*length(alwaysincludclin))


#### 3
# Randomly permute rows of the table containing the settings.
# This is performed to ensure a comparable computational burden for
# the jobs to be performed in parallel:

set.seed(1234)
reorderind <- sample(1:nrow(scenariogrid))
scenariogrid <- scenariogrid[reorderind,]

#### 4
# Save scenariogrid, needed in evaluation of the results:

save(scenariogrid, file="./VarselCompStudy/Additional_file/Results/scenariogrid.Rda")

#### 5
# Source the functions that are used in performing the calculations 
# on the cluster:

source("./VarselCompStudy/Additional_file/Functions/Functions_AnalysisCluster.R")

#### 6
# Start the cluster:

# NOTE: This syntax requires the use of the RMPISNOW script, see the README file
# contained in the root folder "Additional_file_2_HornungWright".

#library(snow)

cl <- makeCluster()


#### 7
# Export the objects in the workspace to the
# parallel jobs:

clusterExport(cl, list=ls())

#### 8
# Perform the calculations:

Results <- parLapply(cl, 1:nrow(scenariogrid), function(z)
  try({evaluatesetting(z)}))

 
save(Results, file="./VarselCompStudy/Additional_file/Results/Results_AnalysisCluster.Rda")

#### 9
# Stop the cluster:

stopCluster(cl)
