rm(list=ls())
library(lavaan)
library(stringr)
library(ggplot2)
library(tidyverse)
library(semPlot)
library(gridExtra)
source("scripts/clean_semplot_functions.R") # script for modified plotting functions
theme_set(theme_bw())


var.max <- 10
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

# analysis ######################################################################################################

# EXTRACT STATS ####
load("data/brute_force_mod_selection_results_parallel.RData")
fit <- sapply(temp, "[", 1)
fit2 <- data.frame(matrix(unlist(fit), nrow=length(fit), byrow=T))
rm(fit)
names(fit2) <-  c("chisq", "df","pvalue","cfi.robust","rmsea.robust","aic")

ggplot(pivot_longer(fit2, cols = names(fit2)), aes(value))+
  geom_histogram()+
  facet_wrap(~name, scales = "free")

model.name <- sapply(temp, "[", 2)
fit2$model.name <- unlist(model.name)
rm(model.name)

# model fit
length(which(fit2$cfi.robust>0.9)) # 812 models have an acceptable robust CFI (with 10 variables)
length(which(fit2$rmsea.robust<0.11)) # 14 with 10 variables

# explanatory power for SR and MRD
pow <- sapply(temp, "[", 3)
rm(temp)
sr.r2 <- sapply(pow, "[[", 1)
mrd.r2 <- sapply(pow, "[[", 2)
fit2$sr.r2 <- sr.r2
fit2$mrd.r2 <- mrd.r2
rm(pow, sr.r2, mrd.r2)

# get models props
## number of variables connected to CFI
# fit2$model.name
sr1 <- sub("\n.*", "", fit2$model.name)
gm.sr <- sub(".*~ ", "", sr1)
gm.sr.sp <- strsplit(gm.sr, "+", fixed = TRUE)
# check number of variables
gm.sr.n <- unlist(lapply(gm.sr.sp, length))
hist(gm.sr.n, xlab="variables in SR regression") # the number of variables in the good models
fit2$gm.sr.n <- gm.sr.n
rm(gm.sr.sp)

mrd1 <- sub(".*\nmrd ~ ", "", fit2$model.name)
gm.mrd <- sub("\nsoil.*", "", mrd1)
gm.mrd.sp <- strsplit(gm.mrd, "+", fixed = TRUE)
gm.mrd.n <- unlist(lapply(gm.mrd.sp, length))
hist(gm.mrd.n, xlab="variables in MRD regression") # the number of variables in the good models
fit2$gm.mrd.n <- gm.mrd.n

# check if pet_m and mat_m appear together in MRD
fit2$mat_in_mrd <- grepl("mat_m", gm.mrd.sp)
fit2$pet_in_mrd <- grepl("pet_m", gm.mrd.sp)
fit2$multicolin_in_sr <- (fit2$mat_in_mrd==TRUE & fit2$pet_in_mrd==TRUE)
table(fit2$multicolin_in_sr)

rm(gm.mrd.sp)

rm(gm.mrd, gm.mrd.n, gm.sr, gm.sr.n)


# PLOTS ####
ggplot(fit2, aes(factor(gm.sr.n), cfi.robust)) +
  geom_boxplot()+
  xlab("number of variables in SR regression")
ggplot(fit2, aes(factor(gm.mrd.n), cfi.robust)) +
  geom_boxplot()+
  xlab("number of variables in MRD regression")
fit2$var.n <- fit2$gm.sr.n+fit2$gm.mrd.n
ggplot(fit2, aes(factor(var.n), cfi.robust)) +
  geom_boxplot()

# ggplot(fit2, aes(gm.mrd.n, gm.sr.n, col=cfi.robust)) +
#   geom_jitter()

#ggplot(fit2, aes(cfi.robust, rmsea.robust)) +
#  geom_point(alpha=0.05)

# AIC CFI correlations
grid.arrange(ncol=2,
ggplot(fit2[fit2$aic<2800,], aes(cfi.robust, aic, col=var.n)) +
 geom_point(alpha=0.1),
ggplot(fit2[fit2$aic<2800,], aes(rmsea.robust, aic, col=var.n)) +
  geom_point(alpha=0.1)
)
cor(fit2$cfi.robust[fit2$aic<2800], fit2$aic[fit2$aic<2800])




# GOOD MODEL FITS ####
## BEST CFI, no influence of variable number/AIC ###########################################
dat_no.na <- readRDS("data/sem_input_data.rds")

min(fit2$rmsea.robust[which(fit2$cfi.robust>0.9)]) # min(rmsea with CFI>0.9 is the same as absolute min RMSEA)
# none of our models achieves an excellent RMSEA, regardless of filtering for CFI before or not

best_cfi <- fit2$model.name[which.max(fit2$cfi.robust)]
cat(best_cfi)
best_cfi.fit <- sem(best_cfi, data = dat_no.na, estimator="MLM")
summary(best_cfi.fit, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
modificationindices(best_cfi.fit, sort. = TRUE)[1:15,]

# MODS on the best model acording to CFI #
fitmeasures(best_cfi.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
best_cfi_mod <- paste0(best_cfi, "\nsub_trop_mbf~area")
best_cfi.fit.mod <- sem(best_cfi_mod, data = dat_no.na, estimator="MLM")
fitmeasures(best_cfi.fit.mod, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_cfi.fit.mod, sort. = TRUE)[1:15,]
best_cfi_mod2 <- paste0(best_cfi_mod, "\nmrd  ~ sub_trop_mbf")
best_cfi.fit.mod2 <- sem(best_cfi_mod2, data = dat_no.na, estimator="MLM")
fitmeasures(best_cfi.fit.mod2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_cfi.fit.mod2, sort. = TRUE)[1:15,]
best_cfi_mod3 <- paste0(best_cfi_mod2, "\nsr_trans ~ pet_m")
best_cfi.fit.mod3 <- sem(best_cfi_mod3, data = dat_no.na, estimator="MLM")
fitmeasures(best_cfi.fit.mod3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_cfi.fit.mod3, sort. = TRUE)[1:15,]

cat(best_cfi_mod2)
summary(best_cfi.fit.mod2, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)

png("best_cfi_mod.png", width=900, height=500, units = "px")
par(mfrow=c(1,2))
semPaths(object = best_cfi.fit, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(best_cfi.fit))
title(round(fitmeasures(best_cfi.fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
semPaths(object = best_cfi.fit.mod2, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(best_cfi.fit.mod2))
title(round(fitmeasures(best_cfi.fit.mod2, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
dev.off()


# BEST RMSEA ###############################################
best_rmsea <- fit2$model.name[which.min(fit2$rmsea.robust)]
cat(best_rmsea)
best_rmsea.fit <- sem(best_rmsea, data = dat_no.na, estimator="MLM")
summary(best_rmsea.fit, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
modificationindices(best_rmsea.fit, sort. = TRUE)[1:15,]

# MODS on the best model according to RMSEA #
fitmeasures(best_rmsea.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
best_rmsea_mod <- paste0(best_rmsea, "\nsub_trop_mbf~area")
best_rmsea.fit.mod <- sem(best_rmsea_mod, data = dat_no.na, estimator="MLM")
fitmeasures(best_rmsea.fit.mod, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_rmsea.fit.mod, sort. = TRUE)[1:15,]
best_rmsea_mod2 <- paste0(best_rmsea_mod, "\nmrd  ~ sub_trop_mbf")
best_rmsea.fit.mod2 <- sem(best_rmsea_mod2, data = dat_no.na, estimator="MLM")
fitmeasures(best_rmsea.fit.mod2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_rmsea.fit.mod2, sort. = TRUE)[1:15,]
best_rmsea_mod3 <- paste0(best_rmsea_mod2, "\nsr_trans ~ pet_m")
best_rmsea.fit.mod3 <- sem(best_rmsea_mod3, data = dat_no.na, estimator="MLM")
fitmeasures(best_rmsea.fit.mod3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_rmsea.fit.mod3, sort. = TRUE)[1:15,]

cat(best_rmsea_mod2)

summary(best_rmsea.fit.mod2, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
png("best_rmsea_mod.png", width=900, height=500, units = "px")
par(mfrow=c(1,2))
semPaths(object = best_rmsea.fit, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(best_rmsea.fit))
title(round(fitmeasures(best_rmsea.fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
semPaths(object = best_rmsea.fit.mod2, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(best_rmsea.fit.mod2))
title(round(fitmeasures(best_rmsea.fit.mod2, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
dev.off()




### GLOBAL MIN AIC, regardless of model fit ###############################################################################
# ---> thats the model we discussed once....
cat(fit2$model.name[which.min(fit2$aic)])
fit2[which.min(fit2$aic),-7]
best_global_aic <- fit2$model.name[which.min(fit2$aic)]
best_global_aic.fit <- sem(best_global_aic, data = dat_no.na, estimator="MLM")
summary(best_global_aic.fit, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
fitmeasures(best_global_aic.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_global_aic.fit, sort. = TRUE)[1:15,]

best_global_aic.mod <- paste0(best_global_aic, "\nsub_trop_mbf~area")
best_global_aic.mod.fit <- sem(best_global_aic.mod, data = dat_no.na, estimator="MLM")
fitmeasures(best_global_aic.mod.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_global_aic.mod.fit, sort. = TRUE)[1:15,]

best_global_aic.mod2 <- paste0(best_global_aic.mod, "\nsub_trop_mbf~mat_m")
best_global_aic.mod.fit2 <- sem(best_global_aic.mod2, data = dat_no.na, estimator="MLM")
fitmeasures(best_global_aic.mod.fit2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_global_aic.mod.fit2, sort. = TRUE)[1:15,]

best_global_aic.mod3 <- paste0(best_global_aic.mod2, "\nsub_trop_mbf~mat_sd")
best_global_aic.mod.fit3 <- sem(best_global_aic.mod3, data = dat_no.na, estimator="MLM")
fitmeasures(best_global_aic.mod.fit3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_global_aic.mod.fit3, sort. = TRUE)[1:15,]

best_global_aic.mod4 <- paste0(best_global_aic.mod3, "\nsr_trans~pet_m")
best_global_aic.mod.fit4 <- sem(best_global_aic.mod4, data = dat_no.na, estimator="MLM")
fitmeasures(best_global_aic.mod.fit4, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_global_aic.mod.fit4, sort. = TRUE)[1:20,]
# the RMSEA is rising....

best_global_aic.mod5 <- paste0(best_global_aic.mod4, "\n mrd ~ sub_trop_mbf")
best_global_aic.mod.fit5 <- sem(best_global_aic.mod5, data = dat_no.na, estimator="MLM")
fitmeasures(best_global_aic.mod.fit5, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_global_aic.mod.fit5, sort. = TRUE)[1:25,]

best_global_aic.mod6 <- paste0(best_global_aic.mod5, "\n sr_trans ~ pre_m")
best_global_aic.mod.fit6 <- sem(best_global_aic.mod6, data = dat_no.na, estimator="MLM")
fitmeasures(best_global_aic.mod.fit6, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_global_aic.mod.fit6, sort. = TRUE)[1:25,]

best_global_aic.mod7 <- paste0(best_global_aic.mod6, "\n mrd ~ mat_m")
best_global_aic.mod.fit7 <- sem(best_global_aic.mod7, data = dat_no.na, estimator="MLM")
fitmeasures(best_global_aic.mod.fit7, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# last step too much

png("best_global_aic_mod.png", width=900, height=500, units = "px")
par(mfrow=c(1,2))
semPaths(object = best_global_aic.fit, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(best_global_aic.fit))
title(round(fitmeasures(best_global_aic.fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
semPaths(object = best_global_aic.mod.fit6, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(best_global_aic.mod.fit6))
title(round(fitmeasures(best_global_aic.mod.fit6, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
dev.off()


## Good models  ##################################################################################################
 # More variables while still having acceptable model fit is desirable since 
dev.off()
hist(fit2$rmsea.robust)
plot(fit2$rmsea.robust[fit2$rmsea.robust<=0.15], fit2$var.n[fit2$rmsea.robust<=0.15])
fin_good <- fit2[fit2$cfi.robust>0.9 & fit2$rmsea.robust<0.13,] # leaving some room for improvement RMSEA, and to capture the best cfi
nrow(fin_good)
# 194 models for 0.12, 609 for 0.13


## Good models with best AIC ranking, select the one with most variables ###################################################

min.aic <- min(fin_good$aic)
fin_good_aic <- fin_good[fin_good$aic<min.aic+2,]
nrow(fin_good_aic) # 20 good models with lowest AIC + max 2 more aic points = set of basically identical AICs
fin_good_aic[,-7]
### most variables
hist(fin_good_aic$var.n)

(m.aic <- max(fin_good_aic$var.n))
fin_good_aic_min <- fin_good_aic[which(fin_good_aic$var.n==m.aic),]
cat(fin_good_aic_min$model.name)
fin_good_aic_min[which.min(fin_good_aic_min$aic),-7]

best_aic <- fin_good_aic_min$model.name
best_aic.fit <- sem(best_aic, data = dat_no.na, estimator="MLM")
summary(best_aic.fit, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
fitmeasures(best_aic.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_aic.fit, sort. = TRUE)[1:15,]
best_aic_mod <- paste0(best_aic, "\nsub_trop_mbf~area")
best_aic.fit.mod <- sem(best_aic_mod, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_aic.fit.mod, sort. = TRUE)[1:15,]
best_aic_mod2 <- paste0(best_aic_mod, "\nsub_trop_mbf~mat_m")
best_aic.fit.mod2 <- sem(best_aic_mod2, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_aic.fit.mod2, sort. = TRUE)[1:15,]
best_aic_mod3 <- paste0(best_aic_mod2, "\nmrd ~ sub_trop_mbf")
best_aic.fit.mod3 <- sem(best_aic_mod3, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_aic.fit.mod3, sort. = TRUE)[1:15,]
best_aic_mod4 <- paste0(best_aic_mod3, "\nsr_trans ~ pet_m")
best_aic.fit.mod4 <- sem(best_aic_mod4, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod4, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_aic.fit.mod4, sort. = TRUE)[1:15,]
best_aic_mod5 <- paste0(best_aic_mod4, "\nsub_trop_mbf ~ sea_m")
best_aic.fit.mod5 <- sem(best_aic_mod5, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod5, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# last step is too much

cat(best_aic_mod4)

png("best_aic_mod.png", width=900, height=500, units = "px")
par(mfrow=c(1,2))
semPaths(object = best_aic.fit, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(best_aic.fit))
title(round(fitmeasures(best_aic.fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
semPaths(object = best_aic.fit.mod4, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(best_aic.fit.mod4))
title(round(fitmeasures(best_aic.fit.mod4, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
dev.off()

summary(best_aic.fit.mod4)

## take care of MULTICOLINEARITY
# How man models have multicollinearity issues (from all fits)?
sr_reg <- 
## Number of models with pet_m + mat_m in SR regression
  # ZERO - pet_m is of such minor importance that it is not considered
## Number of models with pet_m + mat_m in MRD regression
  table(fit2$multicolin_in_sr)
  hist(fit2$cfi.robust[fit2$multicolin_in_sr==T])
  hist(fit2$cfi.robust[fit2$multicolin_in_sr==T & fit2$cfi.robust>=0.9])
  hist(fit2$rmsea.robust[fit2$multicolin_in_sr==T & fit2$rmsea.robust<=0.13])
## Number of models with pet_m + mat_m in TRF regression
  # ZERO, that path is pixed to include pet_m. 
  # is that really justfied:
  lm_trf <- lm(sub_trop_mbf ~ pet_m + tra_m + pre_m, data=dat_no.na)
  lm_trf2 <- lm(sub_trop_mbf ~ mat_m + tra_m + pre_m, data=dat_no.na)
  summary(lm_trf2)
  # --> mat_m is the better variable: higher R squared, avoiding highly non-significant variable. pet does not even have correlation with subtropmbf
  
## adjust best model in correct order without adding pet_m


### minimum AIC in good models, regardless of number of variables###############################################################################
cat(fin_good_aic$model.name[which.min(fin_good_aic$aic)])
fin_good_aic[which.min(fin_good_aic$aic),-7]
best_abs_aic <- fin_good_aic$model.name[which.min(fin_good_aic$aic)]
best_abs_aic.fit <- sem(best_abs_aic, data = dat_no.na, estimator="MLM")
summary(best_abs_aic.fit, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
fitmeasures(best_abs_aic.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_abs_aic.fit, sort. = TRUE)[1:15,]

best_abs_aic_mod <- paste0(best_abs_aic, "\nsub_trop_mbf~area")
best_abs_aic.fit.mod <- sem(best_abs_aic_mod, data = dat_no.na, estimator="MLM")
fitmeasures(best_abs_aic.fit.mod, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_abs_aic.fit.mod, sort. = TRUE)[1:15,]

best_abs_aic_mod2 <- paste0(best_abs_aic_mod, "\nsub_trop_mbf~mat_m")
best_abs_aic.fit.mod2 <- sem(best_abs_aic_mod2, data = dat_no.na, estimator="MLM")
fitmeasures(best_abs_aic.fit.mod2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_abs_aic.fit.mod2, sort. = TRUE)[1:15,]

best_abs_aic_mod3 <- paste0(best_abs_aic_mod2, "\nmrd ~ sub_trop_mbf")
best_abs_aic.fit.mod3 <- sem(best_abs_aic_mod3, data = dat_no.na, estimator="MLM")
fitmeasures(best_abs_aic.fit.mod3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_abs_aic.fit.mod3, sort. = TRUE)[1:15,]

best_abs_aic_mod4 <- paste0(best_abs_aic_mod3, "\nmrd ~ mat_m")
best_abs_aic.fit.mod4 <- sem(best_abs_aic_mod4, data = dat_no.na, estimator="MLM")
fitmeasures(best_abs_aic.fit.mod4, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# too much
cat(best_abs_aic_mod3)

png("best_absolute_aic_mod.png", width=900, height=500, units = "px")
par(mfrow=c(1,2))
semPaths(object = best_abs_aic.fit, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(best_abs_aic.fit))
title(round(fitmeasures(best_abs_aic.fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
semPaths(object = best_abs_aic.fit.mod3, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(best_abs_aic.fit.mod3))
title(round(fitmeasures(best_abs_aic.fit.mod3, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
dev.off()





# Max variable number in good models ####################################################################################
(max.var <- max(fin_good$var.n))
fin_good_max.var <- fin_good[fin_good$var.n==max.var,]
nrow(fin_good_max.var) # its just one model

cat(fin_good_max.var$model.name)
max.var.mod <- fin_good_max.var$model.name

max.var.mod.fit <- sem(max.var.mod, data = dat_no.na, estimator="MLM")
summary(max.var.mod.fit, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
fitmeasures(max.var.mod.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit, sort. = TRUE)[1:15,]

max.var.mod_mod <- paste0(max.var.mod, "\nsub_trop_mbf~area")
max.var.mod.fit.mod <- sem(max.var.mod_mod, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit.mod, sort. = TRUE)[1:15,]

max.var.mod_mod2 <- paste0(max.var.mod_mod, "\n mrd  ~ sub_trop_mbf")
max.var.mod.fit.mod2 <- sem(max.var.mod_mod2, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit.mod2, sort. = TRUE)[1:20,]

max.var.mod_mod3 <- paste0(max.var.mod_mod2, "\n soil  ~ tri")
max.var.mod.fit.mod3 <- sem(max.var.mod_mod3, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit.mod3, sort. = TRUE)

max.var.mod_mod4 <- paste0(max.var.mod_mod3, "\nsr_trans~pet_m")
max.var.mod.fit.mod4 <- sem(max.var.mod_mod4, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod4, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit.mod4, sort. = TRUE)
# last step makes AIC go up again

cat(max.var.mod_mod3)

png("max_var.png", width=500, height=500, units = "px")
semPaths(object = max.var.mod.fit.mod3, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(max.var.mod.fit.mod3))
title(round(fitmeasures(max.var.mod.fit.mod3, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
dev.off()







# Highest combined Rsquared in good models ################################################################################
hist(fin_good$sr.r2)
hist(fin_good$mrd.r2)
plot(fin_good$sr.r2, fin_good$mrd.r2)
fin_good$comb.r2 <- fin_good$sr.r2+ fin_good$mrd.r2
best_r2 <- fin_good$model.name[which.max(fin_good$comb.r2)]
cat(best_r2)
fin_good[which.max(fin_good$comb.r2), -7]

max.r2.fit <- sem(best_r2, data = dat_no.na, estimator="MLM")
summary(max.r2.fit, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
fitmeasures(max.r2.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.r2.fit, sort. = TRUE)[1:15,]

best_r2_mod <- paste0(best_r2, "\nsub_trop_mbf~area")
max.r2.fit.mod <- sem(best_r2_mod, data = dat_no.na, estimator="MLM")
fitmeasures(max.r2.fit.mod, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.r2.fit.mod, sort. = TRUE)[1:15,]

best_r2_mod2 <- paste0(best_r2_mod, "\nsub_trop_mbf~mat_m")
max.r2.fit.mod2 <- sem(best_r2_mod2, data = dat_no.na, estimator="MLM")
fitmeasures(max.r2.fit.mod2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.r2.fit.mod2, sort. = TRUE)[1:15,]

best_r2_mod3 <- paste0(best_r2_mod2, "\nmrd~sub_trop_mbf")
max.r2.fit.mod3 <- sem(best_r2_mod3, data = dat_no.na, estimator="MLM")
fitmeasures(max.r2.fit.mod3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.r2.fit.mod3, sort. = TRUE)[1:15,]

best_r2_mod4 <- paste0(best_r2_mod3, "\nsr_trans~pet_m")
max.r2.fit.mod4 <- sem(best_r2_mod4, data = dat_no.na, estimator="MLM")
fitmeasures(max.r2.fit.mod4, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.r2.fit.mod4, sort. = TRUE)[1:15,]

best_r2_mod5 <- paste0(best_r2_mod4, "\nsub_trop_mbf~sea_m")
max.r2.fit.mod5 <- sem(best_r2_mod5, data = dat_no.na, estimator="MLM")
fitmeasures(max.r2.fit.mod5, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.r2.fit.mod4, sort. = TRUE)[1:15,]
# last step is too much
cat(best_r2_mod4)

png("best_r2.png", width=900, height=500, units = "px")
par(mfrow=c(1,2))
semPaths(object = max.r2.fit, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(max.r2.fit))
title(round(fitmeasures(max.r2.fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
semPaths(object = max.r2.fit.mod4, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(max.r2.fit.mod4))
title(round(fitmeasures(max.r2.fit.mod4, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
dev.off()
















# MORE INDIRECT PATHS ############################################################################################################
b <- gregexpr("\n", fit2$model.name)
d <- lapply(b, length)
rm(b)
fit2$reg.n <- unlist(d)+1

ggplot(fit2, aes(factor(reg.n), cfi.robust)) +
  geom_boxplot()+
  xlab("Number of regressions")

mcfi <- max(fit2$cfi.robust[which(fit2$reg.n==6)])
which(fit2$cfi.robust==mcfi & fit2$reg.n==6)
max_reg_cif_mod <- fit2$model.name[which(fit2$cfi.robust==mcfi & fit2$reg.n==6)]
cat(max_reg_cif_mod)
# remove sea_sd
max_reg_cif_mod<- gsub("\nsea_sd ~ area", "", max_reg_cif_mod)
max_reg_cif_mod<-gsub("\nmrd ~ pre_m + sea_m + tra_m+sea_sd+soil+tra_sd+pet_m+area", "\nmrd ~ pre_m + sea_m + tra_m+soil+tra_sd+pet_m+area", max_reg_cif_mod)

max.reg.fit <- sem(max_reg_cif_mod, data = dat_no.na, estimator="MLM")
fitmeasures(max.reg.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.reg.fit, sort. = TRUE)[1:15,]

# modifications are somewhat questionable.... these variables are correlated, but need to be non causal if included as correlations. think thats hard to justify, that the mean and sd of one variable are non-causally correlated. however why would the mean predict the sd?
max_reg_cif_mod2 <- paste0(max_reg_cif_mod, "\nsea_sd~sea_m")
max.reg.fit2 <- sem(max_reg_cif_mod2, data = dat_no.na, estimator="MLM")
fitmeasures(max.reg.fit2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.reg.fit2, sort. = TRUE)[1:15,]

max_reg_cif_mod3 <- paste0(max_reg_cif_mod2, "\ntra_sd~tra_m")
max.reg.fit3 <- sem(max_reg_cif_mod3, data = dat_no.na, estimator="MLM")
fitmeasures(max.reg.fit3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.reg.fit3, sort. = TRUE)[1:15,]

semPaths(object = max.reg.fit3, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(max.reg.fit3))
title(round(fitmeasures(max.reg.fit3, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)




# sea_sd only influences mrd, and the effect is almost 0. # ignore sea_sd
mcfi2<- max(fit2$cfi.robust[which(fit2$reg.n==5)])
which(fit2$cfi.robust==mcfi2 & fit2$reg.n==5)
max_reg_cif_mod2 <- fit2$model.name[which(fit2$cfi.robust==mcfi2 & fit2$reg.n==5)]
cat(max_reg_cif_mod2)


# FINAL MODEL, closerr inspection ##################################################################################
cat(best_aic_mod4)
fin_model <- "
  sr_trans ~ soil + sub_trop_mbf + area + mrd + mont_gs + tra_m + mat_m + sea_m + pet_m
  mrd ~ pre_m + sea_m + tra_m + soil + pet_m + area + sub_trop_mbf
  sub_trop_mbf ~ pre_m + tra_m + pet_m + area + mat_m
  soil ~ area
"
fin_model_fit <- sem(fin_model, data = dat_no.na, estimator="MLM")
fitmeasures(fin_model_fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))


semPaths(object = fin_model_fit, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(fin_model_fit))
title(round(fitmeasures(fin_model_fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)

# pet and mat have contrasting influence on TRf and SR, thats a bit suspicious
# multicollinearty in the model:
lm_sr <- lm(sr_trans ~ soil + sub_trop_mbf + area + mrd + mont_gs + tra_m + mat_m + sea_m, data=dat_no.na) # removing pet_m solves problem
sort(car::vif(lm_sr))

lm_mrd <- lm(sub_trop_mbf ~ pre_m + tra_m + pet_m + area, data=dat_no.na) # removing mat_m solves issue
sort(car::vif(lm_mrd))

fin_model_nomulti <- "
  sr_trans ~ soil + sub_trop_mbf + area + mrd + mont_gs + tra_m + mat_m + sea_m
  mrd ~ pre_m + sea_m + tra_m + soil + pet_m + area + sub_trop_mbf
  sub_trop_mbf ~ pre_m + tra_m + pet_m + area
  soil ~ area
"
fin_model_nomulti_fit <- sem(fin_model_nomulti, data = dat_no.na, estimator="MLM")

fin_model_nomulti2 <- "
  sr_trans ~ soil + sub_trop_mbf + area + mrd + mont_gs + tra_m + mat_m + sea_m
  mrd ~ pre_m + sea_m + tra_m + soil + pet_m + area + sub_trop_mbf
  sub_trop_mbf ~ pre_m + tra_m + mat_m + area
  soil ~ area
"
fin_model_nomulti2_fit <- sem(fin_model_nomulti2, data = dat_no.na, estimator="MLM")
fitmeasures(fin_model_nomulti2_fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))

# compare 
par(mfrow=c(1,2))
semPaths(object = fin_model_fit, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(fin_model_fit))
title(round(fitmeasures(fin_model_fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
semPaths(object = fin_model_nomulti_fit, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(fin_model_nomulti_fit))
title(round(fitmeasures(fin_model_nomulti_fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
semPaths(object = fin_model_nomulti2_fit, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(fin_model_nomulti2_fit))
title(round(fitmeasures(fin_model_nomulti2_fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)



# decision:
# removing pet_m from SR regression does not cause changes in regression coefficients: keep pet_m in SR for better model fit and actually showing the (lack of)influence
# removing pet_or mat_m does not affect MRD regression
# removing both pet_m or mat_m from TRF regression makes turn tra_m significant! does not affect SR or MRD
## removing mat_m makes pet_m lose importance, but stays NONsign
## removing pet_m makes mat_m lose importance, but stays sign.
## PET_m has to leave TRF regression since its non-significant anyway
## can we really keep it in the SR regression? it has a significant negative influence, while mat_m is positive - ?











# Rsquared main regressions versus number of parameters included in regressions
fit2$sr.r2 <- sr.r2
fit2$mrd.r2 <- mrd.r2
rm(mrd.r2, sr.r2)
ggplot(fit2, aes(sr.r2, mrd.r2))+
  geom_point(alpha=0.01)

ggplot(fit2, aes(sr.r2, gm.sr.n))+
  geom_point(alpha=0.01)

ggplot(fit2, aes(mrd.r2, gm.mrd.n))+
  geom_point(alpha=0.01)


cat(aic_mods$model.name[max(aic_mods$var.n)])
# plot model 
# models from CFI>0.9, RMSEA<0.12, with AIC filtered to only inlcude the minimum+2 points, model with most variables

# model with best Max variable number, filtered for high CFI (>0.9), low RMSEA (<0.12) 
# i=461697
mod_mods <- paste(fin$model.name, "\nsub_trop_mbf~area")
sem_fit <- sem(fin$model.name, data = dat_no.na, estimator="MLM")
sem_fit_mod <- sem(mod_mods, data = dat_no.na, estimator="MLM")
summary(sem_fit)
library(semPlot)
source("scripts/clean_semplot_functions.R") # script for modified plotting functions
fitmeasures(sem_fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(sem_fit, sort. = TRUE)[1:15,]
semPaths(object = sem_fit, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(sem_fit))
semPaths(object = sem_fit_mod, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(sem_fit_mod))





