# Data prep for SEM
# produces data_for_SEM.rds

library(caret)
library(gbm)
library(parallel)
n_cores=1

dat_no.na <- readRDS("data_for_SEM.RDS")
dat_no.na <- dat_no.na[,-grep("sr_trans", names(dat_no.na))]
names(dat_no.na)

# get number of models to run
args <- commandArgs()
print(args)
nmodels <- as.numeric(args[6])
this_start <- as.numeric(args[7])
var_type <- args[8] 


#### GBM CONTROLS #############################################################

# parameter tuning using caret
gbmControl <- trainControl(method = "repeatedcv", number = 10,
                           repeats = 3, savePredictions = "final",
                           returnResamp ="final")

gbmGrid <-  expand.grid(interaction.depth = c(1, 2, 3), 
                        n.trees = (1:10)*50, 
                        shrinkage = c(0.1, 0.01, 0.001),
                        n.minobsinnode = c(5, 10, 15))


if(var_type=="sr"){
# SR ----------------------------------------------------------------------

# define function
run.gbm <- function(i, seeds, data, trControl=gbmControl,
                    tuneGrid=gbmGrid){
  set.seed(s[i]) 
  if(!i%%1)cat(i,"\r")
  gbm_temp <- train(sr ~ ., data, method = "gbm", 
                  trControl = gbmControl, verbose=FALSE,
                  tuneGrid = gbmGrid)
  temp <- list(summary(gbm_temp)$var, s[i])
  return(temp)
}

s <- seq(540,645,1)
s <- s[this_start:(this_start+nmodels)]
print(s)
system.time(
  sr_list <- mclapply(1:length(s), seeds = s, run.gbm,
                      data = dat_no.na[,-grep("level3", names(dat_no.na))],
                      #response=sr,
                      mc.cores=n_cores)
)
saveRDS(sr_list, paste0("sr_list100_lgm", nmodels, "_", this_start,".rds"))


}


if(var_type=="mrd"){
# MRD ---------------------------------------------------------------------

# redefine function for MRD
run.gbm <- function(i, seeds, data, trControl=gbmControl,
                    tuneGrid=gbmGrid){
  set.seed(s[i]) 
  if(!i%%1)cat(i,"\r")
  gbm_temp <- train(mrd ~ ., data, method = "gbm", 
                    trControl = gbmControl, verbose=FALSE,
                    tuneGrid = gbmGrid)
  temp <- list(summary(gbm_temp)$var, s[i])
  return(temp)
}

s <- seq(670,770,1)
s <- s[this_start:(this_start+nmodels)]
system.time(
  mrd_list <- mclapply(1:length(s), seeds = s, run.gbm,
                       data = dat_no.na[,-grep("sr|level3|mdr",names(dat_no.na))],
                       mc.cores=n_cores)
)
saveRDS(mrd_list, paste0("mrd_list100_lgm", nmodels, "_", this_start,".rds"))

}


if(var_type=="mdr"){
# DivRate -----------------------------------------------------------------

# redefine function for MDR
run.gbm <- function(i, seeds, data, trControl=gbmControl,
                    tuneGrid=gbmGrid){
  set.seed(s[i]) 
  if(!i%%1)cat(i,"\r")
  gbm_temp <- train(mdr ~ ., data, method = "gbm", 
                    trControl = gbmControl, verbose=FALSE,
                    tuneGrid = gbmGrid)
  temp <- list(summary(gbm_temp)$var, s[i])
  return(temp)
}

s <- seq(770,870,1)
s <- s[this_start:(this_start+nmodels)]
system.time(
  mdr_list <- mclapply(1:length(s), seeds = s, run.gbm,
                       data = dat_no.na[,-grep("sr|level3|mrd",names(dat_no.na))],
                       mc.cores=n_cores)
)
saveRDS(mdr_list, paste0("mdr_list100_lgm", nmodels, "_", this_start,".rds"))
}




