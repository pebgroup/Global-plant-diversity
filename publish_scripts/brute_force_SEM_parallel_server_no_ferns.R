# BRUTE FORCE VARIABLE SELECTION ####################################################################################
rm(list=ls())
library(lavaan)
library(stringr)
library(ggplot2)
library(parallel)
library(foreach)
library(doParallel)
dat_no.na <- readRDS("data/sem_input_data_noferns.rds")

# set max number of variables (to add to the fixed ones) used and register cores
var.max <- 10
registerDoParallel(16)

# build regression character strings ###############################################################################
sr_variables <- c("sub_trop_dbf",  "pre_m", "pre_sd", "sea_sd", "mat_sd", "tra_sd" 
                  , "mat_m", "mont_gs" ,"tra_m", "medit_fws", "tri", "pet_m", "temp_gss" , "pet_sd" , "sea_m" )
mrd_variables <- c("tra_m", "sub_trop_mbf" ,"sea_sd", "tra_sd" , "pre_sd" , "soil" , "pet_m", "mat_m" , "tri",
                   "mat_sd" ,"area", "temp_bmf" ,  "pet_sd" , "deserts_x_shrub", "mont_gs")
fixed_sr <- c("sr_trans ~ soil + sub_trop_mbf + area + mrd +")
fixed_mrd <- c("mrd ~ pre_m + sea_m + tra_m + sub_trop_mbf +")

sr_var_list <- do.call("c", lapply(seq_along(sr_variables[1:var.max]), function(i) combn(sr_variables[1:var.max], i, FUN = list)))
sr_combs <- c()
Sys.time() # TAKES 20 secs
for(i in 1:length(sr_var_list)){ #
  temp <- unlist(sr_var_list[i]) 
  sr_combs <- c(sr_combs, str_c(temp, collapse = "+"))
}
Sys.time()
mrd_var_list <- do.call("c", lapply(seq_along(mrd_variables[1:var.max]), function(i) combn(mrd_variables[1:var.max], i, FUN = list)))
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
subtrop_path <- "sub_trop_mbf ~ pre_m + tra_m + mat_m"  # always present
desert_path <- "deserts_x_shrub ~ pre_m + sea_m + pet_m"  # if habitat included

# loop through model options  _ TAKES TIME! #####################################################################################

#### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### ~ 5 hours for 2^9 * 2^9 combinations for serial analysis
#### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


system.time({
# r2list <- list()
# mod.list <- list()
temp <- foreach(i = 1:length(main_reg)) %dopar% {
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
    #if(i==1){
      res <- fitmeasures(modfit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
    #}else{
    #  res <- rbind(res, fitmeasures(modfit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic")))
    #}
    #r2list[[i]] <- lavInspect(modfit, "r2")
    #mod.list[[i]] <- mod
    list(res, mod, lavInspect(modfit, "r2")[1:2])
  }
})

save(temp, file="brute_force_mod_selection_results_parallel_mat_noferns.RData")

