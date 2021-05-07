# analyse output from server brute force runs
# produces brute_force_models_new_order.RData
# next script: sem_results.rmd
rm(list=ls())
library(lavaan)
library(stringr)
library(ggplot2)
library(tidyverse)
library(semPlot)
library(gridExtra)
source("scripts/clean_semplot_functions.R") # script for modified plotting functions


var.max <- 10
# build regression character strings ###############################################################################
sr_variables <- c("sub_trop_dbf",  "pre_sd", "mont_gs", "mat_sd", "sea_sd", "tra_m", "pre_m", "tri", "pet_sd", "tra_sd",
                  "mat_m", "pet_m", "sea_m")
mrd_variables <- c("mat_m", "tra_m", "tri", "soil", "sea_sd", "tra_sd", "pre_sd", "temp_bmf", "pet_m",  "pet_sd",   
                   "mat_sd", "area")
fixed_sr <- c("sr_trans ~ soil + sub_trop_mbf + area + mrd +")
fixed_mrd <- c("mrd ~ pre_m + sub_trop_mbf + sea_m +")

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

# possible variables unique display
all_vars <- c(sr_variables[1:10], mrd_variables[1:10],
              unlist(strsplit(fixed_sr, "+", fixed = T)), 
              unlist(strsplit(fixed_mrd, "+", fixed = T)))
all_vars <- (unlist(strsplit(all_vars,"~",fixed=T)))
unique(sort(gsub(" ", "", all_vars)))

# build all options combining the two regression combos (use)
temp <- expand.grid(sr_combs, mrd_combs)
main_reg = paste(paste(fixed_sr,temp$Var1), paste(fixed_mrd, temp$Var2), sep="\n")

# set indirect path dependencies ##################################################################################

dir_path_soil_area <- "soil ~ area"                     # always present
dir_path_sea_sd_area <- "sea_sd ~ area"                 # grepl("sea_sd", mod)
dir_path_tra_sd_area <- "tra_sd ~ area"                 # grepl("tra_sd", mod)
subtrop_path <- "sub_trop_mbf ~ pre_m + tra_m + pet_m"  # always present
desert_path <- "deserts_x_shrub ~ pre_m + sea_m + pet_m"  # if habitat included




# EXTRACT STATS ####
load("processed_data/brute_force_mod_selection_Apr2021.RData")
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
length(which(fit2$cfi.robust>0.9)) # 257 with acceptable cfi
length(which(fit2$rmsea.robust<0.12)) # 4 with RMSEA < 0.12

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
fit2$multicolin_in_mrd <- (fit2$mat_in_mrd==TRUE & fit2$pet_in_mrd==TRUE)
table(fit2$multicolin_in_mrd)

rm(gm.mrd.sp, gm.mrd, gm.mrd.n, gm.sr, gm.sr.n)


ggplot(fit2, aes(factor(gm.sr.n), cfi.robust)) +
  geom_boxplot()+
  xlab("number of variables in SR regression")
ggplot(fit2, aes(factor(gm.mrd.n), cfi.robust)) +
  geom_boxplot()+
  xlab("number of variables in MRD regression")
fit2$var.n <- fit2$gm.sr.n+fit2$gm.mrd.n
ggplot(fit2, aes(factor(var.n), cfi.robust)) +
  geom_boxplot()


# AIC CFI correlations
grid.arrange(ncol=2,
             ggplot(fit2[fit2$aic<2800,], aes(cfi.robust, aic, col=var.n)) +
               geom_point(alpha=0.1),
             ggplot(fit2[fit2$aic<2800,], aes(rmsea.robust, aic, col=var.n)) +
               geom_point(alpha=0.1)
)






# MODEL MANUAL FITTING ####


## 1) BEST CFI, no influence of variable number/AIC ###########################################
dat_no.na <- readRDS("processed_data/sem_input_data_noferns.rds")

min(fit2$rmsea.robust)
min(fit2$rmsea.robust[which(fit2$cfi.robust>0.9)]) # min(rmsea with CFI>0.9 is the same as absolute min RMSEA)
# none of our models achieves an excellent RMSEA, regardless of filtering for CFI before or not

best_cfi <- fit2$model.name[which.max(fit2$cfi.robust)]
cat(best_cfi)
best_cfi.fit <- sem(best_cfi, data = dat_no.na, estimator="MLM")
summary(best_cfi.fit, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
fitmeasures(best_cfi.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_cfi.fit, sort. = TRUE)[1:15,]

# MODS on the best model according to CFI #
fitmeasures(best_cfi.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
best_cfi_mod <- paste0(best_cfi, "\nsub_trop_mbf~area")
best_cfi.fit.mod <- sem(best_cfi_mod, data = dat_no.na, estimator="MLM")
fitmeasures(best_cfi.fit.mod, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_cfi.fit.mod, sort. = TRUE)[1:15,]
best_cfi_mod2 <- paste0(best_cfi_mod, "\nmrd  ~ area")
best_cfi.fit.mod2 <- sem(best_cfi_mod2, data = dat_no.na, estimator="MLM")
fitmeasures(best_cfi.fit.mod2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_cfi.fit.mod2, sort. = TRUE)[1:15,]
best_cfi_mod3 <- paste0(best_cfi_mod2, "\nsr_trans ~ mat_m")
best_cfi.fit.mod3 <- sem(best_cfi_mod3, data = dat_no.na, estimator="MLM")
fitmeasures(best_cfi.fit.mod3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_cfi.fit.mod3, sort. = TRUE)[1:15,]
best_cfi_mod4 <- paste0(best_cfi_mod3, "\nsr_trans ~ sea_m")
best_cfi.fit.mod4 <- sem(best_cfi_mod4, data = dat_no.na, estimator="MLM")
fitmeasures(best_cfi.fit.mod4, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_cfi.fit.mod4, sort. = TRUE)[1:15,]
best_cfi_mod5 <- paste0(best_cfi_mod4, "\nmrd ~ tra_m")
best_cfi.fit.mod5 <- sem(best_cfi_mod5, data = dat_no.na, estimator="MLM")
fitmeasures(best_cfi.fit.mod5, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# rmsea goes up
modificationindices(best_cfi.fit.mod5, sort. = TRUE)[1:20,]
best_cfi_mod6 <- paste0(best_cfi_mod5, "\nsub_trop_mbf ~ sea_m")
best_cfi.fit.mod6 <- sem(best_cfi_mod6, data = dat_no.na, estimator="MLM")
fitmeasures(best_cfi.fit.mod6, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# AIC goes up

cat(best_cfi_mod4)


png("publish_figures/best_cfi_model.png", width=900, height=500, units = "px")
par(mfrow=c(1,2))
semPaths(object = best_cfi.fit, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(best_cfi.fit))
title(round(fitmeasures(best_cfi.fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
semPaths(object = best_cfi.fit.mod4, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(best_cfi.fit.mod4))
title(round(fitmeasures(best_cfi.fit.mod4, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
dev.off()


# 2) BEST RMSEA ###############################################
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
best_rmsea_mod2 <- paste0(best_rmsea_mod, "\nsoil  ~ mont_gs")
best_rmsea.fit.mod2 <- sem(best_rmsea_mod2, data = dat_no.na, estimator="MLM")
fitmeasures(best_rmsea.fit.mod2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_rmsea.fit.mod2, sort. = TRUE)[1:15,]
best_rmsea_mod3 <- paste0(best_rmsea_mod2, "\nmrd ~ area")
best_rmsea.fit.mod3 <- sem(best_rmsea_mod3, data = dat_no.na, estimator="MLM")
fitmeasures(best_rmsea.fit.mod3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_rmsea.fit.mod3, sort. = TRUE)[1:15,]
best_rmsea_mod4 <- paste0(best_rmsea_mod3, "\nsub_trop_mbf~tri")
best_rmsea.fit.mod4 <- sem(best_rmsea_mod4, data = dat_no.na, estimator="MLM")
fitmeasures(best_rmsea.fit.mod4, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_rmsea.fit.mod4, sort. = TRUE)[1:15,]
best_rmsea_mod5 <- paste0(best_rmsea_mod4, "\nsr_trans ~ mat_m")
best_rmsea.fit.mod5 <- sem(best_rmsea_mod5, data = dat_no.na, estimator="MLM")
fitmeasures(best_rmsea.fit.mod5, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_rmsea.fit.mod5, sort. = TRUE)[1:15,]
best_rmsea_mod6 <- paste0(best_rmsea_mod5, "\nsr_trans ~ sea_m")
best_rmsea.fit.mod6 <- sem(best_rmsea_mod6, data = dat_no.na, estimator="MLM")
fitmeasures(best_rmsea.fit.mod6, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_rmsea.fit.mod6, sort. = TRUE)[1:15,]
best_rmsea_mod7 <- paste0(best_rmsea_mod6, "\nsr_trans ~ tra_m")
best_rmsea.fit.mod7 <- sem(best_rmsea_mod7, data = dat_no.na, estimator="MLM")
fitmeasures(best_rmsea.fit.mod7, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_rmsea.fit.mod7, sort. = TRUE)[1:15,]
best_rmsea_mod8 <- paste0(best_rmsea_mod7, "\nmrd ~ tra_m")
best_rmsea.fit.mod8 <- sem(best_rmsea_mod8, data = dat_no.na, estimator="MLM")
fitmeasures(best_rmsea.fit.mod8, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_rmsea.fit.mod8, sort. = TRUE)[1:15,]
best_rmsea_mod9 <- paste0(best_rmsea_mod8, "\nsub_trop_mbf~sea_m")
best_rmsea.fit.mod9 <- sem(best_rmsea_mod9, data = dat_no.na, estimator="MLM")
fitmeasures(best_rmsea.fit.mod9, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# too much

cat(best_rmsea_mod8)

png("publish_figures/best_rmsea_model.png", width=900, height=500, units = "px")
par(mfrow=c(1,2))
semPaths(object = best_rmsea.fit, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(best_rmsea.fit))
title(round(fitmeasures(best_rmsea.fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
semPaths(object = best_rmsea.fit.mod6, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(best_rmsea.fit.mod6))
title(round(fitmeasures(best_rmsea.fit.mod6, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
dev.off()




### 3) GLOBAL MIN AIC, regardless of model fit ###################################
cat(fit2$model.name[which.min(fit2$aic)])
fit2[which.min(fit2$aic),-7]
best_global_aic <- fit2$model.name[which.min(fit2$aic)]
best_global_aic.fit <- sem(best_global_aic, data = dat_no.na, estimator="MLM")
#summary(best_global_aic.fit, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)

fitmeasures(best_global_aic.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_global_aic.fit, sort. = TRUE)[1:15,]
best_global_aic.mod <- paste0(best_global_aic, "\nsr_trans ~ temp_bmf")
best_global_aic.mod.fit <- sem(best_global_aic.mod, data = dat_no.na, estimator="MLM")
fitmeasures(best_global_aic.mod.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_global_aic.mod.fit, sort. = TRUE)[1:15,]
best_global_aic.mod2 <- paste0(best_global_aic.mod, "\nsub_trop_mbf~area")
best_global_aic.mod.fit2 <- sem(best_global_aic.mod2, data = dat_no.na, estimator="MLM")
fitmeasures(best_global_aic.mod.fit2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_global_aic.mod.fit2, sort. = TRUE)[1:20,]
best_global_aic.mod3 <- paste0(best_global_aic.mod2, "\nsr_trans~mat_m")
best_global_aic.mod.fit3 <- sem(best_global_aic.mod3, data = dat_no.na, estimator="MLM")
fitmeasures(best_global_aic.mod.fit3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_global_aic.mod.fit3, sort. = TRUE)[1:25,]
best_global_aic.mod4 <- paste0(best_global_aic.mod3, "\nmrd~area")
best_global_aic.mod.fit4 <- sem(best_global_aic.mod4, data = dat_no.na, estimator="MLM")
fitmeasures(best_global_aic.mod.fit4, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# RMSEA GOES UP
modificationindices(best_global_aic.mod.fit4, sort. = TRUE)[1:25,]
best_global_aic.mod5 <- paste0(best_global_aic.mod4, "\nsub_trop_mbf~tri")
best_global_aic.mod.fit5 <- sem(best_global_aic.mod5, data = dat_no.na, estimator="MLM")
fitmeasures(best_global_aic.mod.fit5, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_global_aic.mod.fit5, sort. = TRUE)[1:25,]
best_global_aic.mod6 <- paste0(best_global_aic.mod5, "\nsr_trans~sea_m")
best_global_aic.mod.fit6 <- sem(best_global_aic.mod6, data = dat_no.na, estimator="MLM")
fitmeasures(best_global_aic.mod.fit6, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_global_aic.mod.fit6, sort. = TRUE)[1:25,]
best_global_aic.mod7 <- paste0(best_global_aic.mod6, "\nmrd~tra_m")
best_global_aic.mod.fit7 <- sem(best_global_aic.mod7, data = dat_no.na, estimator="MLM")
fitmeasures(best_global_aic.mod.fit7, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_global_aic.mod.fit7, sort. = TRUE)[1:25,]
# RMSEA just keeps going up, no more meaningful changes

png("publish_figures/best_global_aic_model.png", width=900, height=500, units = "px")
par(mfrow=c(1,2))
semPaths(object = best_global_aic.fit, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(best_global_aic.fit))
title(round(fitmeasures(best_global_aic.fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
semPaths(object = best_global_aic.mod.fit7, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(best_global_aic.mod.fit7))
title(round(fitmeasures(best_global_aic.mod.fit7, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
dev.off()


## 4) pre-select good CFI and RMSEA#####################################################################
# More variables while still having acceptable model fit is desirable since 
dev.off()
plot(fit2$rmsea.robust[fit2$rmsea.robust<=0.15], fit2$var.n[fit2$rmsea.robust<=0.15])
fin_good <- fit2[fit2$cfi.robust>0.9 & fit2$rmsea.robust<0.13,] # leaving some room for improvement RMSEA, and to capture the best cfi

nrow(fin_good) # models with potential


## pre-select, best AIC ranking, most variables ######################
min.aic <- min(fin_good$aic)
fin_good_aic <- fin_good[fin_good$aic<min.aic+2,]
nrow(fin_good_aic) # 27# 21 good models with lowest AIC + max 2 more aic points = set of basically identical AICs
fin_good_aic[,-7]
### most variables
hist(fin_good_aic$var.n)

(m.aic <- max(fin_good_aic$var.n)) # 12 variables
fin_good_aic_min <- fin_good_aic[which(fin_good_aic$var.n==m.aic),] 
nrow(fin_good_aic_min) # 12 good models
cat(fin_good_aic_min$model.name)
fin_good_aic_min[which.min(fin_good_aic_min$aic),-7]

# select the one with lowest AIC
best_aic <- fin_good_aic_min$model.name[which.min(fin_good_aic_min$aic)]
best_aic.fit <- sem(best_aic, data = dat_no.na, estimator="MLM")
summary(best_aic.fit, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
fitmeasures(best_aic.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_aic.fit, sort. = TRUE)[1:15,]
best_aic_mod <- paste0(best_aic, "\nsub_trop_mbf~area")
best_aic.fit.mod <- sem(best_aic_mod, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_aic.fit.mod, sort. = TRUE)[1:15,]
best_aic_mod2 <- paste0(best_aic_mod, "\nmrd ~ area")
best_aic.fit.mod2 <- sem(best_aic_mod2, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_aic.fit.mod2, sort. = TRUE)[1:15,]
best_aic_mod3 <- paste0(best_aic_mod2, "\nsub_trop_mbf~tri")
best_aic.fit.mod3 <- sem(best_aic_mod3, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_aic.fit.mod3, sort. = TRUE)[1:15,]
best_aic_mod4 <- paste0(best_aic_mod3, "\nsr_trans~mat_m")
best_aic.fit.mod4 <- sem(best_aic_mod4, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod4, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_aic.fit.mod4, sort. = TRUE)[1:15,]
# RMSEA goes up
best_aic_mod5 <- paste0(best_aic_mod4, "\nsr_trans~sea_m")
best_aic.fit.mod5 <- sem(best_aic_mod5, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod5, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_aic.fit.mod5, sort. = TRUE)[1:15,]
best_aic_mod6 <- paste0(best_aic_mod5, "\nsr_trans~tri")
best_aic.fit.mod6 <- sem(best_aic_mod6, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod6, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_aic.fit.mod6, sort. = TRUE)[1:15,]
best_aic_mod7 <- paste0(best_aic_mod6, "\nsr_trans~tra_m")
best_aic.fit.mod7 <- sem(best_aic_mod7, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod7, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_aic.fit.mod7, sort. = TRUE)[1:15,]
best_aic_mod8 <- paste0(best_aic_mod7, "\nsub_trop_mbf~sea_m")
best_aic.fit.mod8 <- sem(best_aic_mod8, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod8, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# CFI down, RMSEA and AIC up

cat(best_aic_mod7)

png("publish_figures/best_aic_model.png", width=900, height=500, units = "px")
par(mfrow=c(1,2))
semPaths(object = best_aic.fit, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(best_aic.fit))
title(round(fitmeasures(best_aic.fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
semPaths(object = best_aic.fit.mod7, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(best_aic.fit.mod7))
title(round(fitmeasures(best_aic.fit.mod7, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
dev.off()




### minimum AIC in good models, regardless of number of variables###############################################################################
#THIS MODEL PERFOMRS WORSE THAN THE ONE WITH MORE VIARBLES; OUT
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

best_abs_aic_mod2 <- paste0(best_abs_aic_mod, "\nmrd~area")
best_abs_aic.fit.mod2 <- sem(best_abs_aic_mod2, data = dat_no.na, estimator="MLM")
fitmeasures(best_abs_aic.fit.mod2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_abs_aic.fit.mod2, sort. = TRUE)[1:15,]

best_abs_aic_mod3 <- paste0(best_abs_aic_mod2, "\nsr_trans ~ mat_m")
best_abs_aic.fit.mod3 <- sem(best_abs_aic_mod3, data = dat_no.na, estimator="MLM")
fitmeasures(best_abs_aic.fit.mod3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_abs_aic.fit.mod3, sort. = TRUE)[1:15,]

best_abs_aic_mod4 <- paste0(best_abs_aic_mod3, "\nsr_trans ~ sea_m")
best_abs_aic.fit.mod4 <- sem(best_abs_aic_mod4, data = dat_no.na, estimator="MLM")
fitmeasures(best_abs_aic.fit.mod4, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_abs_aic.fit.mod4, sort. = TRUE)[1:15,]

best_abs_aic_mod5 <- paste0(best_abs_aic_mod4, "\nmrd ~ tra_m")
best_abs_aic.fit.mod5 <- sem(best_abs_aic_mod5, data = dat_no.na, estimator="MLM")
fitmeasures(best_abs_aic.fit.mod5, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_abs_aic.fit.mod5, sort. = TRUE)[1:15,]
# RMSEA goes up
best_abs_aic_mod6 <- paste0(best_abs_aic_mod5, "\nsr_trans ~ tra_m")
best_abs_aic.fit.mod6 <- sem(best_abs_aic_mod6, data = dat_no.na, estimator="MLM")
fitmeasures(best_abs_aic.fit.mod6, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# aic goes up

cat(best_abs_aic_mod5)

png("publish_figures/best_absolute_aic_model.png", width=900, height=500, units = "px")
par(mfrow=c(1,2))
semPaths(object = best_abs_aic.fit, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(best_abs_aic.fit))
title(round(fitmeasures(best_abs_aic.fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
semPaths(object = best_abs_aic.fit.mod5, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(best_abs_aic.fit.mod5))
title(round(fitmeasures(best_abs_aic.fit.mod5, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
dev.off()





# Max variable number in good models ####################################################################################
(max.var <- max(fin_good$var.n))
fin_good_max.var <- fin_good[fin_good$var.n==max.var,]
nrow(fin_good_max.var) # good models with max var

cat(fin_good_max.var$model.name)
# multicolineariy? no. take the one with best fit
fin_good_max.var[,-7]
max.var.mod <- fin_good_max.var$model.name[which.min(fin_good_max.var$rmsea.robust)]

max.var.mod.fit <- sem(max.var.mod, data = dat_no.na, estimator="MLM")
summary(max.var.mod.fit, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
fitmeasures(max.var.mod.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit, sort. = TRUE)[1:15,]

max.var.mod_mod <- paste0(max.var.mod, "\nsoil~mont_gs")
max.var.mod.fit.mod <- sem(max.var.mod_mod, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit.mod, sort. = TRUE)[1:15,]

max.var.mod_mod2 <- paste0(max.var.mod_mod, "\n sub_trop_mbf~area")
max.var.mod.fit.mod2 <- sem(max.var.mod_mod2, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit.mod2, sort. = TRUE)[1:15,]

max.var.mod_mod3 <- paste0(max.var.mod_mod2, "\n mrd~area")
max.var.mod.fit.mod3 <- sem(max.var.mod_mod3, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit.mod3, sort. = TRUE)[1:15,]

max.var.mod_mod4 <- paste0(max.var.mod_mod3, "\nsub_trop_mbf  ~ tri")
max.var.mod.fit.mod4 <- sem(max.var.mod_mod4, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod4, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit.mod4, sort. = TRUE)[1:15,]

max.var.mod_mod5 <- paste0(max.var.mod_mod4, "\nsr_trans~mat_m")
max.var.mod.fit.mod5 <- sem(max.var.mod_mod5, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod5, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit.mod5, sort. = TRUE)[1:15,]

max.var.mod_mod6 <- paste0(max.var.mod_mod5, "\nsr_trans~sea_m")
max.var.mod.fit.mod6 <- sem(max.var.mod_mod6, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod6, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit.mod6, sort. = TRUE)[1:15,]

max.var.mod_mod7 <- paste0(max.var.mod_mod6, "\nsr_trans~tri")
max.var.mod.fit.mod7 <- sem(max.var.mod_mod7, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod7, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit.mod7, sort. = TRUE)[1:15,]

max.var.mod_mod8 <- paste0(max.var.mod_mod7, "\nsub_trop_mbf~sea_m")
max.var.mod.fit.mod8 <- sem(max.var.mod_mod8, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod8, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# CFI down, RMSEA + AIC up

cat(max.var.mod_mod7)

png("publish_figures/max_var_model.png", width=900, height=500, units = "px")
par(mfrow=c(1,2))
semPaths(object = max.var.mod.fit, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(max.var.mod.fit))
title(round(fitmeasures(max.var.mod.fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
semPaths(object = max.var.mod.fit.mod7, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(max.var.mod.fit.mod7))
title(round(fitmeasures(max.var.mod.fit.mod7, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
dev.off()



# 
# 
# # Highest combined Rsquared in good models ################################################################################
# hist(fin_good$sr.r2)
# hist(fin_good$mrd.r2)
# plot(fin_good$sr.r2, fin_good$mrd.r2)
# fin_good$comb.r2 <- fin_good$sr.r2+ fin_good$mrd.r2
# best_r2 <- fin_good$model.name[which.max(fin_good$comb.r2)]
# cat(best_r2)
# fin_good[which.max(fin_good$comb.r2), -7]
# 
# max.r2.fit <- sem(best_r2, data = dat_no.na, estimator="MLM")
# summary(max.r2.fit, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
# fitmeasures(max.r2.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# modificationindices(max.r2.fit, sort. = TRUE)[1:15,]
# 
# best_r2_mod <- paste0(best_r2, "\nsub_trop_mbf~area")
# max.r2.fit.mod <- sem(best_r2_mod, data = dat_no.na, estimator="MLM")
# fitmeasures(max.r2.fit.mod, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# modificationindices(max.r2.fit.mod, sort. = TRUE)[1:15,]
# 
# best_r2_mod2 <- paste0(best_r2_mod, "\nmrd ~ sub_trop_mbf")
# max.r2.fit.mod2 <- sem(best_r2_mod2, data = dat_no.na, estimator="MLM")
# fitmeasures(max.r2.fit.mod2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# modificationindices(max.r2.fit.mod2, sort. = TRUE)[1:15,]
# 
# best_r2_mod3 <- paste0(best_r2_mod2, "\nsr_trans~mat_m")
# max.r2.fit.mod3 <- sem(best_r2_mod3, data = dat_no.na, estimator="MLM")
# fitmeasures(max.r2.fit.mod3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# modificationindices(max.r2.fit.mod3, sort. = TRUE)[1:20,]
# 
# best_r2_mod4 <- paste0(best_r2_mod3, "\nsr_trans~sea_m")
# max.r2.fit.mod4 <- sem(best_r2_mod4, data = dat_no.na, estimator="MLM")
# fitmeasures(max.r2.fit.mod4, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# modificationindices(max.r2.fit.mod4, sort. = TRUE)[1:20,]
# 
# best_r2_mod5 <- paste0(best_r2_mod4, "\nsoil~sub_trop_mbf")
# max.r2.fit.mod5 <- sem(best_r2_mod5, data = dat_no.na, estimator="MLM")
# fitmeasures(max.r2.fit.mod5, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# modificationindices(max.r2.fit.mod5, sort. = TRUE)[1:20,]
# # RMSEA goes up, no use
# cat(best_r2_mod4)
# 
# png("best_r2_new_order.png", width=900, height=500, units = "px")
# par(mfrow=c(1,2))
# semPaths(object = max.r2.fit, layout = "circle2", rotation = 1,
#          whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE,
#          edgeLabels = sem_sig_labels(max.r2.fit))
# title(round(fitmeasures(max.r2.fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
# semPaths(object = max.r2.fit.mod4, layout = "circle2", rotation = 1,
#          whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE,
#          edgeLabels = sem_sig_labels(max.r2.fit.mod4))
# title(round(fitmeasures(max.r2.fit.mod4, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
# dev.off()
# 




save.image("processed_data/models.RData")






# 
# # MORE INDIRECT PATHS ############################################################################################################
# b <- gregexpr("\n", fit2$model.name)
# d <- lapply(b, length)
# rm(b)
# fit2$reg.n <- unlist(d)+1
# 
# ggplot(fit2, aes(factor(reg.n), cfi.robust)) +
#   geom_boxplot()+
#   xlab("Number of regressions")
# 
# mcfi <- max(fit2$cfi.robust[which(fit2$reg.n==6)])
# which(fit2$cfi.robust==mcfi & fit2$reg.n==6)
# max_reg_cif_mod <- fit2$model.name[which(fit2$cfi.robust==mcfi & fit2$reg.n==6)]
# cat(max_reg_cif_mod)
# # remove sea_sd
# max_reg_cif_mod<- gsub("\nsea_sd ~ area", "", max_reg_cif_mod)
# max_reg_cif_mod<-gsub("\nmrd ~ pre_m + sea_m + tra_m+sea_sd+soil+tra_sd+pet_m+area", "\nmrd ~ pre_m + sea_m + tra_m+soil+tra_sd+pet_m+area", max_reg_cif_mod)
# 
# max.reg.fit <- sem(max_reg_cif_mod, data = dat_no.na, estimator="MLM")
# fitmeasures(max.reg.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# modificationindices(max.reg.fit, sort. = TRUE)[1:15,]
# 
# # modifications are somewhat questionable.... these variables are correlated, but need to be non causal if included as correlations. think thats hard to justify, that the mean and sd of one variable are non-causally correlated. however why would the mean predict the sd?
# max_reg_cif_mod2 <- paste0(max_reg_cif_mod, "\nsea_sd~sea_m")
# max.reg.fit2 <- sem(max_reg_cif_mod2, data = dat_no.na, estimator="MLM")
# fitmeasures(max.reg.fit2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# modificationindices(max.reg.fit2, sort. = TRUE)[1:15,]
# 
# max_reg_cif_mod3 <- paste0(max_reg_cif_mod2, "\ntra_sd~tra_m")
# max.reg.fit3 <- sem(max_reg_cif_mod3, data = dat_no.na, estimator="MLM")
# fitmeasures(max.reg.fit3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# modificationindices(max.reg.fit3, sort. = TRUE)[1:15,]
# 
# semPaths(object = max.reg.fit3, layout = "circle2", rotation = 1,
#          whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE,
#          edgeLabels = sem_sig_labels(max.reg.fit3))
# title(round(fitmeasures(max.reg.fit3, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
# 
# 
# 
# 
# # sea_sd only influences mrd, and the effect is almost 0. # ignore sea_sd
# mcfi2<- max(fit2$cfi.robust[which(fit2$reg.n==5)])
# which(fit2$cfi.robust==mcfi2 & fit2$reg.n==5)
# max_reg_cif_mod2 <- fit2$model.name[which(fit2$cfi.robust==mcfi2 & fit2$reg.n==5)]
# cat(max_reg_cif_mod2)
# 
# 
# # BEST MODELS, closer inspection ##################################################################################
# cat(max.var.mod_mod5)
# cat(best_aic_mod3)
# 
# 
# fin_model <- "
#   sr_trans ~ soil + sub_trop_mbf + area + mrd + mont_gs + tra_m + mat_m + sea_m + pet_m
#   mrd ~ pre_m + sea_m + tra_m + soil + pet_m + area + sub_trop_mbf
#   sub_trop_mbf ~ pre_m + tra_m + pet_m + area + mat_m
#   soil ~ area
# "
# fin_model_fit <- sem(fin_model, data = dat_no.na, estimator="MLM")
# fitmeasures(fin_model_fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# 
# 
# semPaths(object = fin_model_fit, layout = "circle2", rotation = 1,
#          whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE,
#          edgeLabels = sem_sig_labels(fin_model_fit))
# title(round(fitmeasures(fin_model_fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
# 
# # pet and mat have contrasting influence on TRf and SR, thats a bit suspicious
# # multicollinearty in the model:
# lm_sr <- lm(sr_trans ~ soil + sub_trop_mbf + area + mrd + mont_gs + tra_m + mat_m + sea_m, data=dat_no.na) # removing pet_m solves problem
# sort(car::vif(lm_sr))
# 
# lm_mrd <- lm(sub_trop_mbf ~ pre_m + tra_m + pet_m + area, data=dat_no.na) # removing mat_m solves issue
# sort(car::vif(lm_mrd))
# 
# fin_model_nomulti <- "
#   sr_trans ~ soil + sub_trop_mbf + area + mrd + mont_gs + tra_m + mat_m + sea_m
#   mrd ~ pre_m + sea_m + tra_m + soil + pet_m + area + sub_trop_mbf
#   sub_trop_mbf ~ pre_m + tra_m + pet_m + area
#   soil ~ area
# "
# fin_model_nomulti_fit <- sem(fin_model_nomulti, data = dat_no.na, estimator="MLM")
# 
# fin_model_nomulti2 <- "
#   sr_trans ~ soil + sub_trop_mbf + area + mrd + mont_gs + tra_m + mat_m + sea_m
#   mrd ~ pre_m + sea_m + tra_m + soil + pet_m + area + sub_trop_mbf
#   sub_trop_mbf ~ pre_m + tra_m + mat_m + area
#   soil ~ area
# "
# fin_model_nomulti2_fit <- sem(fin_model_nomulti2, data = dat_no.na, estimator="MLM")
# fitmeasures(fin_model_nomulti2_fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# 
# # compare 
# par(mfrow=c(1,2))
# semPaths(object = fin_model_fit, layout = "circle2", rotation = 1,
#          whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE,
#          edgeLabels = sem_sig_labels(fin_model_fit))
# title(round(fitmeasures(fin_model_fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
# semPaths(object = fin_model_nomulti_fit, layout = "circle2", rotation = 1,
#          whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE,
#          edgeLabels = sem_sig_labels(fin_model_nomulti_fit))
# title(round(fitmeasures(fin_model_nomulti_fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
# semPaths(object = fin_model_nomulti2_fit, layout = "circle2", rotation = 1,
#          whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE,
#          edgeLabels = sem_sig_labels(fin_model_nomulti2_fit))
# title(round(fitmeasures(fin_model_nomulti2_fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
# 
# 
# 
# 
# 
# 
# 
# 




