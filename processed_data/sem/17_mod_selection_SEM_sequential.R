# VARIABLE SELECTION ##########################################################
# 1 h on 32 cores
### ~~~ make sure you set the working directory right and have all files in the
### same directory as this script. Do NOT run this script directly but use the
### bash script <mod_selection_master_sh> ~~~ ####

library(lavaan)
library(doMC)
library(stringr)

# set max number of variables (to add to the fixed ones) used and register cores
# problems with the cluster do not allow parallel processing, happens outside in
# bash scripts for now
var.max <- 10
#n.cpus <- Sys.getenv("SLURM_CPUS_PER_TASK")
#n.cpus <- as.numeric(n.cpus)
n.cpus <- 1
registerDoMC(cores = n.cpus)


# build regression character strings ####
fixed_sr <- c("sr_trans ~ soil + sub_trop_mbf + area + mrd +")
fixed_mrd <- c("mrd ~ pre_m + sub_trop_mbf + prs_m +")

#selected_variables <- readRDS("selected_variables.rds")
#names(selected_variables[[1]][1:15])
sr_variables <- c("pre_sd", "sub_trop_dbf", "elev_range", "pre_m", "mat_lgm_ano_m",
                  "mont_gs", "mat_m", "mat_sd", "prs_sd", "tra_m")

names(selected_variables[[2]][1:15])
mrd_variables <- c("mio_mat_ano_m", "mat_m", "elev_range", "prs_sd", "soil", "tra_m",
                   "tra_sd", "tri", "pre_sd", "temp_bmf")

sr_var_list <- do.call("c", lapply(seq_along(sr_variables[1:var.max]), 
                                   function(i) combn(sr_variables[1:var.max], i, FUN = list)))
sr_combs <- c()
print("Building SR equations...")
for(i in 1:length(sr_var_list)){ #
  temp <- unlist(sr_var_list[i]) 
  sr_combs <- c(sr_combs, str_c(temp, collapse = "+"))
}
print("Building MRD equations...")
mrd_var_list <- do.call("c", lapply(seq_along(mrd_variables[1:var.max]),
                                    function(i) combn(mrd_variables[1:var.max], i, FUN = list)))
mrd_combs <- c()
for(i in 1:length(mrd_var_list)){ #
  temp <- unlist(mrd_var_list[i]) 
  mrd_combs <- c(mrd_combs, str_c(temp, collapse = "+"))
}


# build all options combining the two regression forumulas
temp <- expand.grid(sr_combs, mrd_combs)
main_reg = paste(paste(fixed_sr,temp$Var1), paste(fixed_mrd, temp$Var2), sep="\n")

# set indirect path dependencies ###############################################################################
dir_path_soil_area <- "soil ~ area"                     # always present
dir_path_prs_sd_area <- "prs_sd ~ area"                 # grepl("prs_sd", mod)
dir_path_tra_sd_area <- "tra_sd ~ area"                 # grepl("tra_sd", mod)
dir_path_mat_sd_area <- "mat_sd ~ area"                 # grepl("mat_sd", mod)
dir_path_pre_sd_area <- "pre_sd ~ area"                 # grepl("pre_sd", mod)
subtrop_path <- "sub_trop_mbf ~ pre_m + tra_m + mat_m"  # always present
#desert_path <- "deserts_x_shrub ~ pre_m + prs_m + pet_m"  # if habitat included

# loop through model options ###############################################################################

#### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### ~ 2^9 * 2^9 combinations for serial analysis ~~~~~~~~~~~
#### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# get number of models to run
args <- commandArgs()
print(args)
nmods <- as.numeric(args[6])
nstart <- as.numeric(args[7])
print(paste("Number of models to run",  nmods))
print(paste("Iteration number to start from",  nstart))

# subset model string
main_reg_sub <- main_reg[nstart:(nstart+nmods)]
rm(main_reg)


dat_no.na <- readRDS("sem_input_data.RDS")

print("Starting sequential loop...")
Sys.time()
system.time(
temp <- foreach(i = 1:nmods,
                .packages=c('lavaan', 'stringr')) %do% {
    # build regression models
    mod_main <- paste(c(main_reg_sub[i], dir_path_soil_area, subtrop_path), collapse = "\n")
    mod <- mod_main
    
    if(grepl("prs_sd", mod)){
      mod <- paste(c(mod, dir_path_prs_sd_area), collapse = "\n")
    }
    if(grepl("tra_sd", mod)){
      mod <- paste(c(mod, dir_path_tra_sd_area), collapse = "\n")
    }
    if(grepl("mat_sd", mod)){
      mod <- paste(c(mod, dir_path_mat_sd_area), collapse = "\n")
    }
    if(grepl("pre_sd", mod)){
      mod <- paste(c(mod, dir_path_pre_sd_area), collapse = "\n")
    }
    
    modfit <- sem(mod, data = dat_no.na, estimator="MLM")
    res <- fitmeasures(modfit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
    list(res, mod, lavInspect(modfit, "r2")[1:2])
  }
)

#jobname <- Sys.getenv("SLURM_CPUS_PER_TASK")
#jobid <- Sys.getenv("SLURM_JOB_ID")
print(paste0("Saving files as SEM_results_Wednesday", nstart, "_", nmods, ".RData"))
save(temp, file=paste0("output/SEM_results_Wednesday", nstart, "_", nmods, ".RData"))


# check
#a <- load("output/SEM_results_101382532705.RData")
#b <- load("output/SEM_results_16352132705.RData")
