# BRUTE FORCE VARIABLE SELECTION ####################################################################################
rm(list=ls())
library(lavaan)
library(stringr)
library(ggplot2)
library(parallel)
library(foreach)
library(doParallel)
library(tidyverse)
dat_no.na <- readRDS("data/sem_input_data.rds")

# set max number of variables (to add to the fixed ones) used and register cores
var.max <- 5
registerDoParallel(4)
loop_active <- 1 # set to length(main_reg) to run the entire loop for actual analysis

# DESCRIPTION ######################################################################################################
# Important! Keep in mind that this checks initial model setup fit, not model fit that might happen when manually following the modifications function e.g. Therefore, this is a method to find a good starting point. Manual adjustments might still improve models, so consider also models with less good fit but more variables!
#   
#   Also important: indirect paths make the model fit worse, filtering for goodness of fit is biased against including these paths. Check models with more indirect paths separately
# 
# How it works
# 1. build fixed parts of regression formulas
# ◦ sr_trans ~ soil + sub_trop_mbf + area + mrd
# ◦ mrd ~ pre_m + sea_m 
# 2. build vectors with potential variables to include in each regression separately
# ◦ sr_variables <- c("sub_trop_dbf",  "pre_sd", "mont_gs" , "pre_m", "sea_sd", "tra_m"                  , "mat_m" ,"tra_sd", "mat_sd" , "sea_m" , "medit_fws" , "temp_gss" , "tri" , "pet_sd" , "pet_m")
# ◦ mrd_variables <- c("tra_m", "sea_sd", "soil" , "tri", "tra_sd" , "pre_sd" , "pet_m", "mat_m" , "area" , "mat_sd" , "temp_bmf" , "sub_trop_mbf" , "pet_sd" , "medit_fws" , "sub_trop_cf" , "deserts_x_shrub")
# 3. build all combinations for each regression formulas
# 4. build all combinations of all combinations each formula
# → you end up with 2^(sr variable number)*2^(mrd variable number)
# 5. define indirect path formulas based on correlations in correlation matrix + knowledge
# ◦ dir_path_soil_area <- "soil ~ area"                   	 			# always present
# ◦ dir_path_sea_sd_area <- "sea_sd ~ area"         	 			# grepl("sea_sd", mod)
# ◦ dir_path_tra_sd_area <- "tra_sd ~ area"            				# grepl("tra_sd", mod)
# ◦ subtrop_path <- "sub_trop_mbf ~ pre_m + tra_m + pet_m"  	# always present
# ◦ desert_path <- "deserts_x_shrub ~ pre_m + sea_m + pet_m"   # if biome present
# 6. include secondary regression (indirect paths) according to if variables are present or not using grepl()
# 7. Collect goodness of fit parameters, Rquared for SR and MRD regressions, model formula
# 8. Run!



# build regression character strings ###############################################################################
sr_variables <- c("sub_trop_dbf",  "pre_sd", "mont_gs" , "pre_m", "sea_sd", "tra_m"
                  , "mat_m" ,"tra_sd", "mat_sd" , "sea_m" , "medit_fws" , "temp_gss" , "tri" , "pet_sd" , "pet_m")
mrd_variables <- c("tra_m", "sea_sd", "soil" , "tri", "tra_sd" , "pre_sd" , "pet_m", "mat_m" , "area"
                   , "mat_sd" , "temp_bmf" , "sub_trop_mbf" , "pet_sd" , "medit_fws" , "sub_trop_cf" , "deserts_x_shrub")
fixed_sr <- c("sr_trans ~ soil + sub_trop_mbf + area + mrd +")
fixed_mrd <- c("mrd ~ pre_m + sea_m +")

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
subtrop_path <- "sub_trop_mbf ~ pre_m + tra_m + pet_m"  # always present
desert_path <- "deserts_x_shrub ~ pre_m + sea_m + pet_m"  # if habitat included

# loop through model options  _ TAKES TIME! #####################################################################################

#### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### ~ 5 hours for 2^9 * 2^9 combinations for serial analysis
#### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sys.time()
system.time({
# r2list <- list()
# mod.list <- list()
temp <- foreach(i = 1:loop_active) %dopar% {
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

#save(temp, file="brute_force_mod_selection_results_parallel.RData")

# Analyse output ##################################################################################################
