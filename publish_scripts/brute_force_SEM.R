# BRUTE FORCE VARIABLE SELECTION ####################################################################################
library(lavaan)
library(stringr)
library(ggplot2)
dat_no.na <- readRDS("data/sem_input_data.rds")

# build regression character strings ###############################################################################
sr_variables <- c("sub_trop_dbf",  "pre_sd", "mrd", "mont_gs" , "pre_m", "sea_sd", "tra_m"
                  , "mat_m" ,"tra_sd", "mat_sd" , "sea_m" , "medit_fws" , "temp_gss" , "tri" , "pet_sd" , "pet_m")
mrd_variables <- c("tra_m", "sea_sd", "soil" , "tri", "tra_sd" , "pre_sd" , "pet_m", "mat_m" , "area"
                   , "mat_sd" , "temp_bmf" , "sub_trop_mbf" , "pet_sd" , "medit_fws" , "sub_trop_cf" , "deserts_x_shrub")
fixed_sr <- c("sr_trans ~ soil + sub_trop_mbf + area +")
fixed_mrd <- c("mrd ~ pre_m + sea_m +")

sr_var_list <- do.call("c", lapply(seq_along(sr_variables[1:9]), function(i) combn(sr_variables[1:9], i, FUN = list)))
sr_combs <- c()
Sys.time() # TAKES 20 secs
for(i in 1:length(sr_var_list)){ #
  temp <- unlist(sr_var_list[i]) 
  sr_combs <- c(sr_combs, str_c(temp, collapse = "+"))
}
Sys.time()
mrd_var_list <- do.call("c", lapply(seq_along(mrd_variables[1:9]), function(i) combn(mrd_variables[1:9], i, FUN = list)))
mrd_combs <- c()
Sys.time() # TAKES 30secs
for(i in 1:length(mrd_var_list)){ #
  temp <- unlist(mrd_var_list[i]) 
  mrd_combs <- c(mrd_combs, str_c(temp, collapse = "+"))
}
Sys.time()

# build all options combining the two regression combos (use)
temp <- expand.grid(sr_combs, mrd_combs)
main_reg = paste(paste(fixed_sr,temp$Var1), paste(fixed_mrd, temp$Var2), sep="\n")
cat(main_reg[54])
# set indirect path dependencies ##################################################################################

dir_path_soil_area <- "soil ~ area"                     # always present
dir_path_sea_sd_area <- "sea_sd ~ area"                 # grepl("sea_sd", mod)
dir_path_tra_sd_area <- "tra_sd ~ area"                 # grepl("tra_sd", mod)
subtrop_path <- "sub_trop_mbf ~ pre_m + tra_m + pet_m"  # always present
desert_path <- "deserts_x_shrub ~ pre_m + sea_m + pet_m"  # if habitat included

# loop through model options  _ TAKES TIME! #####################################################################################

#### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### ~ 5 hours for 2^9 * 2^9 combinations. 
#### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sys.time()
r2list <- list()
mod.list <- list()
for(i in 1:length(main_reg)){
  # build regression models
  mod_main <- paste(c(main_reg[i], dir_path_soil_area, subtrop_path), collapse = "\n")
  mod <- mod_main
  #cat(main_reg)
  
  if(grepl("sea_sd", mod_main)){
    mod <- paste(c(mod_main, dir_path_sea_sd_area), collapse = "\n")
  }
  if(grepl("tra_sd", mod_main)){
    mod <- paste(c(mod_main, dir_path_tra_sd_area), collapse = "\n")
  }
  if(grepl("sea_sd", mod_main) & grepl("tra_sd", mod_main)){
    mod <- paste(c(mod_main, dir_path_sea_sd_area, dir_path_tra_sd_area), collapse = "\n")
  }
  if(grepl("deserts_x_shrub", mod_main)){
    mod <- paste(c(mod_main, dir_path_sea_sd_area, dir_path_tra_sd_area, desert_path), collapse = "\n")
  }
  
  modfit <- sem(mod, data = dat_no.na, estimator="MLM")
  if(i==1){
    res <- fitmeasures(modfit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
  }else{
    res <- rbind(res, fitmeasures(modfit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic")))
  }
  r2list[[i]] <- lavInspect(modfit, "r2")
  mod.list[[i]] <- mod
  
  if(!i%%100)cat(i,"\r")
}
Sys.time()
save(res, r2list, mod.list, file="data/brute_force_mod_selection_results.RData")

# Analyse output ##################################################################################################

head(res)
which.max(res[,4]) # max CFI
which.min(res[,5]) # min RMSEA
cat(mod.list[[50081]])
r2list[50081]

# model fit
out <- as.data.frame(res)
hist(out)
length(which(out$cfi.robust>0.9))/nrow(out) # 753 models have an acceptable robust CFI
length(which(out$rmsea.robust<0.11))/nrow(out) # 22

# explanatory power
hist(unlist(lapply(r2list, length)))
## extract first two
sr.r2 <- sapply(r2list, "[[", 1)
mrd.r2 <- sapply(r2list, "[[", 2)
hist(sr.r2)
hist(mrd.r2)

# get models props
## number of variables connected to CFI
good.mods <- unlist(mod.list)
sr1 <- sub("\n.*", "", good.mods)
gm.sr <- sub(".*~ ", "", sr1)
gm.sr.sp <- strsplit(gm.sr, "+", fixed = TRUE)
gm.sr.n <- unlist(lapply(gm.sr.sp, length))
hist(gm.sr.n, xlab="variables in SR regression") # the number of variables in the good models

mrd1 <- sub(".*\nmrd ~ ", "", good.mods)
gm.mrd <- sub("\nsoil.*", "", mrd1)
gm.mrd.sp <- strsplit(gm.mrd, "+", fixed = TRUE)
gm.mrd.n <- unlist(lapply(gm.mrd.sp, length))
hist(gm.mrd.n, xlab="variables in MRD regression") # the number of variables in the good models
out$gm.sr.n <- gm.sr.n
out$gm.mrd.n <- gm.mrd.n
ggplot(out[out$cfi.robust>0.8, ], aes(factor(gm.sr.n), cfi.robust)) + 
  geom_boxplot()
ggplot(out, aes(factor(gm.mrd.n), cfi.robust)) + 
  geom_boxplot()
out$var.n <- out$gm.mrd.n + out$gm.sr.n
ggplot(out, aes(factor(var.n), cfi.robust)) + 
  geom_boxplot()

ggplot(out[out$cfi.robust>0.8, ], aes(gm.mrd.n, gm.sr.n, col=cfi.robust)) + 
  geom_jitter()

ggplot(out, aes(cfi.robust, rmsea.robust)) + 
  geom_point(alpha=0.1)

# find model that has most variables with best fit
max(out$var.n[which(out$cfi.robust>0.9)])
which(out$cfi.robust>0.9 & out$var.n==15)
fin <- unlist(mod.list[which(out$cfi.robust>0.9 & out$var.n==15)])
cat(fin[1])
cat(fin[2])


