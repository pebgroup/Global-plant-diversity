# get stats from model selection runs, manually adjust models
# produces models.RData
# next script: sem_results.rmd

wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
library(lavaan)
library(stringr)
library(ggplot2)
library(tidyr)
library(semPlot)
source("0_functions.R")
theme_set(theme_bw())
library(cowplot)


# Load data DivRate -----------------------------------------------------------

rm(list = setdiff(ls(), lsf.str())) 
filnames <- dir("../processed_data/sem_output", full.names = T)
filnames <- filnames[grep("results_DivRate", filnames)]
length(filnames)
myfun <- function(x){
  a <- get(load(x))
  return(a)
}
models <- sapply(filnames, myfun, simplify = F) 
# combine files
models <- unlist(models, recursive = F)


# EXTRACT STATS ####
fit <- sapply(models, "[", 1)
fit2 <- data.frame(matrix(unlist(fit), nrow=length(fit), byrow=T))
rm(fit)
names(fit2) <-  c("chisq", "df","pvalue","cfi.robust","rmsea.robust","aic")

ggplot(pivot_longer(fit2, cols = names(fit2)), aes(value))+
  geom_histogram()+
  facet_wrap(~name, scales = "free")

model.name <- sapply(models, "[", 2)
models[grep("NA", model.name)] # is fine now....
# first 2 models of batch 32705 are ok, third one is missing the main regression

fit2$model.name <- unlist(model.name)
rm(model.name)

# model fit
length(which(fit2$cfi.robust>0.9)) # 772 with acceptable CFI
length(which(fit2$rmsea.robust<0.11)) # 112 with RMSEA < 0.11

# explanatory power for SR and mdr
pow <- sapply(models, "[", 3)
rm(models)
sr.r2 <- sapply(pow, "[[", 1)
mdr.r2 <- sapply(pow, "[[", 2)
fit2$sr.r2 <- sr.r2
fit2$mdr.r2 <- mdr.r2
rm(pow, sr.r2, mdr.r2)

# get models props
## number of variables connected to CFI
# fit2$model.name
sr1 <- sub("\n.*", "", fit2$model.name)
gm.sr <- sub(".*~ ", "", sr1)
gm.sr.sp <- strsplit(gm.sr, "+", fixed = TRUE)
# check number of variables
gm.sr.n <- unlist(lapply(gm.sr.sp, length))
#hist(gm.sr.n, xlab="variables in SR regression") # the number of variables in the good models
fit2$gm.sr.n <- gm.sr.n
rm(gm.sr.sp)

mdr1 <- sub(".*\nmdr ~ ", "", fit2$model.name)
gm.mdr <- sub("\nsoil.*", "", mdr1)
gm.mdr.sp <- strsplit(gm.mdr, "+", fixed = TRUE)
gm.mdr.n <- unlist(lapply(gm.mdr.sp, length))
#hist(gm.mdr.n, xlab="variables in mdr regression") # the number of variables in the good models
fit2$gm.mdr.n <- gm.mdr.n

# check if pet_m and mat_m appear together in mdr
fit2$mat_in_mdr <- grepl("mat_m", gm.mdr.sp)
fit2$pet_in_mdr <- grepl("pet_m", gm.mdr.sp)
fit2$multicolin_in_mdr <- (fit2$mat_in_mdr==TRUE & fit2$pet_in_mdr==TRUE)
table(fit2$multicolin_in_mdr)

rm(gm.mdr.sp, gm.mdr, gm.mdr.n, gm.sr, gm.sr.n)


ggplot(fit2, aes(factor(gm.sr.n), cfi.robust)) +
  geom_boxplot()+
  xlab("number of variables in SR regression")
ggplot(fit2, aes(factor(gm.mdr.n), cfi.robust)) +
  geom_boxplot()+
  xlab("number of variables in mdr regression")
fit2$var.n <- fit2$gm.sr.n+fit2$gm.mdr.n
ggplot(fit2, aes(factor(var.n), cfi.robust)) +
  geom_boxplot()






# MODEL MANUAL FITTING ########################################################################


## 1) BEST CFI, no influence of variable number/AIC ###########################################
dat_no.na <- readRDS("../processed_data/sem_input_data.rds")

min(fit2$rmsea.robust)
min(fit2$rmsea.robust[which(fit2$cfi.robust>0.9)]) # min(rmsea with CFI>0.9 is the same as absolute min RMSEA)
# none of our models achieves an excellent RMSEA, regardless of filtering for CFI before or not

best_cfi <- fit2$model.name[which.max(fit2$cfi.robust)]
cat(best_cfi)
best_cfi.fit <- sem(best_cfi, data = dat_no.na, estimator="MLM")
summary(best_cfi.fit, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
fitmeasures(best_cfi.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_cfi.fit, sort. = TRUE)[1:15,]

fitmeasures(best_cfi.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
best_cfi.mod <-  paste0(best_cfi, "\nsub_trop_mbf ~ area")
best_cfi.fit.mod <- sem(best_cfi.mod, data = dat_no.na, estimator="MLM")
fitmeasures(best_cfi.fit.mod, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_cfi.fit.mod, sort. = TRUE)[1:15,]

best_cfi.mod2 <-  paste0(best_cfi.mod, "\nmdr~sub_trop_mbf")
best_cfi.fit.mod2 <- sem(best_cfi.mod2, data = dat_no.na, estimator="MLM")
fitmeasures(best_cfi.fit.mod2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_cfi.fit.mod2, sort. = TRUE)[1:15,]

best_cfi.mod3 <-  paste0(best_cfi.mod2, "\nmdr ~ tra_m")
best_cfi.fit.mod3 <- sem(best_cfi.mod3, data = dat_no.na, estimator="MLM")
fitmeasures(best_cfi.fit.mod3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_cfi.fit.mod3, sort. = TRUE)[1:15,]

best_cfi.mod4 <-  paste0(best_cfi.mod3, "\nsr_trans  ~ boreal_f_taiga")
best_cfi.fit.mod4 <- sem(best_cfi.mod4, data = dat_no.na, estimator="MLM")
fitmeasures(best_cfi.fit.mod4, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_cfi.fit.mod4, sort. = TRUE)[1:15,]

best_cfi.fit.mod5 <- update(best_cfi.fit.mod4, add="sr_trans ~ tra_m")
fitmeasures(best_cfi.fit.mod5, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# rmsea goes up
# modificationindices(best_cfi.fit.mod5, sort. = TRUE)[1:20,]
# 
# best_cfi.fit.mod6 <- update(best_cfi.fit.mod4, add="sr_trans ~ tra_m")
# fitmeasures(best_cfi.fit.mod6, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# # all stable
# modificationindices(best_cfi.fit.mod6, sort. = TRUE)[1:20,]
# 
# best_cfi.fit.mod7 <- update(best_cfi.fit.mod6, add="mdr ~ tra_m")
# fitmeasures(best_cfi.fit.mod7, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# # RMSEA goes up

model <- best_cfi.fit.mod4
png("../figures/DivRate_best_cfi_model_2022.png", width=500, height=500, units = "px")
semPaths(object = model, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(model))
title(round(fitmeasures(model, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
dev.off()


# 2) BEST RMSEA ###############################################
best_rmsea <- fit2$model.name[which.min(fit2$rmsea.robust)]
cat(best_rmsea)
best_rmsea.fit <- sem(best_rmsea, data = dat_no.na, estimator="MLM")
summary(best_rmsea.fit, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
modificationindices(best_rmsea.fit, sort. = TRUE)[1:15,]

# MODS on the best model according to RMSEA #
fitmeasures(best_rmsea.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
best_rmsea_mod <- paste0(best_rmsea, "\nsub_trop_mbf  ~ area") #mdr~area
best_rmsea.fit.mod <- sem(best_rmsea_mod, data = dat_no.na, estimator="MLM")
fitmeasures(best_rmsea.fit.mod, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_rmsea.fit.mod, sort. = TRUE)[1:15,]
best_rmsea_mod2 <- paste0(best_rmsea_mod, "\nsub_trop_mbf  ~  mio_pre_ano_m")
best_rmsea.fit.mod2 <- sem(best_rmsea_mod2, data = dat_no.na, estimator="MLM")
fitmeasures(best_rmsea.fit.mod2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_rmsea.fit.mod2, sort. = TRUE)[1:15,]
best_rmsea_mod3 <- paste0(best_rmsea_mod2, "\nmdr ~ tra_m")
best_rmsea.fit.mod3 <- sem(best_rmsea_mod3, data = dat_no.na, estimator="MLM")
fitmeasures(best_rmsea.fit.mod3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_rmsea.fit.mod3, sort. = TRUE)[1:15,]
best_rmsea_mod4 <- paste0(best_rmsea_mod3, "\nmdr ~ sub_trop_mbf")
best_rmsea.fit.mod4 <- sem(best_rmsea_mod4, data = dat_no.na, estimator="MLM")
fitmeasures(best_rmsea.fit.mod4, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_rmsea.fit.mod4, sort. = TRUE)[1:15,]
best_rmsea_mod5 <- paste0(best_rmsea_mod4, "\nsr_trans  ~ boreal_f_taiga")
best_rmsea.fit.mod5 <- sem(best_rmsea_mod5, data = dat_no.na, estimator="MLM")
fitmeasures(best_rmsea.fit.mod5, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_rmsea.fit.mod5, sort. = TRUE)[1:20,]
best_rmsea_mod6 <- paste0(best_rmsea_mod5, "\nsr_trans  ~      medit_fws")
best_rmsea.fit.mod6 <- sem(best_rmsea_mod6, data = dat_no.na, estimator="MLM")
fitmeasures(best_rmsea.fit.mod6, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_rmsea.fit.mod6, sort. = TRUE)[1:25,]
best_rmsea_mod7 <- paste0(best_rmsea_mod6, "\nsr_trans  ~          pre_m")
best_rmsea.fit.mod7 <- sem(best_rmsea_mod7, data = dat_no.na, estimator="MLM")
fitmeasures(best_rmsea.fit.mod7, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
#AIC+RMSEA going up
#modificationindices(best_rmsea.fit.mod7, sort. = TRUE)[1:20,]
# stop, nothing meaninful to add


cat(best_rmsea_mod6)

png("../figures/DivRatebest_rmsea_model2022.png", width=500, height=500, units = "px")
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

fitmeasures(best_global_aic.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# modificationindices(best_global_aic.fit, sort. = TRUE)[1:15,]
# best_global_aic.mod <- paste0(best_global_aic, "\nsub_trop_mbf  ~ elev_range")
# best_global_aic.mod.fit <- sem(best_global_aic.mod, data = dat_no.na, estimator="MLM")
# fitmeasures(best_global_aic.mod.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# modificationindices(best_global_aic.mod.fit, sort. = TRUE)[1:15,]
# best_global_aic.mod2 <- paste0(best_global_aic.mod, "\nmdr~area")
# best_global_aic.mod.fit2 <- sem(best_global_aic.mod2, data = dat_no.na, estimator="MLM")
# fitmeasures(best_global_aic.mod.fit2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# modificationindices(best_global_aic.mod.fit2, sort. = TRUE)[1:20,]
# best_global_aic.mod3 <- paste0(best_global_aic.mod2, "\nsub_trop_mbf~area")
# best_global_aic.mod.fit3 <- sem(best_global_aic.mod3, data = dat_no.na, estimator="MLM")
# fitmeasures(best_global_aic.mod.fit3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# modificationindices(best_global_aic.mod.fit3, sort. = TRUE)[1:25,]
# best_global_aic.mod4 <- paste0(best_global_aic.mod3, "\nsr_trans~mat_m")
# best_global_aic.mod.fit4 <- sem(best_global_aic.mod4, data = dat_no.na, estimator="MLM")
# fitmeasures(best_global_aic.mod.fit4, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# modificationindices(best_global_aic.mod.fit4, sort. = TRUE)[1:25,]
# best_global_aic.mod5 <- paste0(best_global_aic.mod4, "\nsr_trans  ~      tra_m")
# best_global_aic.mod.fit5 <- sem(best_global_aic.mod5, data = dat_no.na, estimator="MLM")
# fitmeasures(best_global_aic.mod.fit5, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# modificationindices(best_global_aic.mod.fit5, sort. = TRUE)[1:25,]
# best_global_aic.mod6 <- paste0(best_global_aic.mod5, "\nsr_trans~prs_m")
# best_global_aic.mod.fit6 <- sem(best_global_aic.mod6, data = dat_no.na, estimator="MLM")
# fitmeasures(best_global_aic.mod.fit6, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# modificationindices(best_global_aic.mod.fit6, sort. = TRUE)[1:25,]
# best_global_aic.mod7 <- paste0(best_global_aic.mod6, "\nsub_trop_mbf  ~        prs_m")
# best_global_aic.mod.fit7 <- sem(best_global_aic.mod7, data = dat_no.na, estimator="MLM")
# fitmeasures(best_global_aic.mod.fit7, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# modificationindices(best_global_aic.mod.fit7, sort. = TRUE)[1:25,]
# # RMSEA just keeps going up, CFI goes down. Stop. 
# 
# png("figures/best_global_aic_model_2022.png", width=900, height=500, units = "px")
# par(mfrow=c(1,2))
# semPaths(object = best_global_aic.fit, layout = "circle2", rotation = 1,
#          whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE,
#          edgeLabels = sem_sig_labels(best_global_aic.fit))
# title(round(fitmeasures(best_global_aic.fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
# semPaths(object = best_global_aic.mod.fit6, layout = "circle2", rotation = 1,
#          whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE,
#          edgeLabels = sem_sig_labels(best_global_aic.mod.fit6))
# title(round(fitmeasures(best_global_aic.mod.fit6, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
# dev.off()
# 

## 4) pre-select good CFI and RMSEA#####################################################################
# More variables while still having acceptable model fit is desirable since 
dev.off()
plot(fit2$rmsea.robust[fit2$rmsea.robust<=0.12], fit2$var.n[fit2$rmsea.robust<=0.12])
fin_good <- fit2[fit2$cfi.robust>0.9 & fit2$rmsea.robust<0.105,] # leaving some room for improvement RMSEA, and to capture the best cfi

nrow(fin_good) # models with potential


## pre-select, best AIC ranking, most variables ######################
min.aic <- min(fin_good$aic)
fin_good_aic <- fin_good[fin_good$aic<min.aic+2,]
nrow(fin_good_aic) # good models with lowest AIC + max 2 more aic points = set of basically identical AICs
fin_good_aic[,-7]
### most variables
hist(fin_good_aic$var.n)

(m.aic <- max(fin_good_aic$var.n)) 
fin_good_aic_min <- fin_good_aic[which(fin_good_aic$var.n==m.aic),] 
nrow(fin_good_aic_min) # 6 good models
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
best_aic_mod2 <- paste0(best_aic_mod, "\nmdr  ~   sub_trop_mbf")
best_aic.fit.mod2 <- sem(best_aic_mod2, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_aic.fit.mod2, sort. = TRUE)[1:15,]
best_aic_mod3 <- paste0(best_aic_mod2, "\nsr_trans  ~          prs_m")
best_aic.fit.mod3 <- sem(best_aic_mod3, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_aic.fit.mod3, sort. = TRUE)[1:15,]
best_aic_mod4 <- paste0(best_aic_mod3, "\nsub_trop_mbf  ~  mio_pre_ano_m")
best_aic.fit.mod4 <- sem(best_aic_mod4, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod4, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# RMSEA goes up
modificationindices(best_aic.fit.mod4, sort. = TRUE)[1:15,]
best_aic_mod5 <- paste0(best_aic_mod4, "\nsr_trans  ~      medit_fws") 
best_aic.fit.mod5 <- sem(best_aic_mod5, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod5, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_aic.fit.mod5, sort. = TRUE)[1:15,]
best_aic_mod6 <- paste0(best_aic_mod5, "\nsr_trans  ~        temp_cf")
best_aic.fit.mod6 <- sem(best_aic_mod6, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod6, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_aic.fit.mod6, sort. = TRUE)[1:15,]
best_aic_mod7 <- paste0(best_aic_mod6, "\nsub_trop_mbf  ~        temp_cf")
best_aic.fit.mod7 <- sem(best_aic_mod7, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod7, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_aic.fit.mod7, sort. = TRUE)[1:20,]
best_aic_mod8 <- paste0(best_aic_mod7, "\nsr_trans  ~  mio_pre_ano_m")
best_aic.fit.mod8 <- sem(best_aic_mod8, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod8, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# # RMSEA up, AIC down
modificationindices(best_aic.fit.mod8, sort. = TRUE)[1:20,]
best_aic_mod9 <- paste0(best_aic_mod8, "\nsr_trans  ~          tra_m")
best_aic.fit.mod9 <- sem(best_aic_mod9, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod9, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
 # STOP
# best_aic_mod5 <- paste0(best_aic_mod4, "\nsr_trans~prs_m")
# best_aic.fit.mod5 <- sem(best_aic_mod5, data = dat_no.na, estimator="MLM")
# fitmeasures(best_aic.fit.mod5, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# modificationindices(best_aic.fit.mod5, sort. = TRUE)[1:15,]
# best_aic_mod6 <- paste0(best_aic_mod5, "\nsr_trans~tri")
# best_aic.fit.mod6 <- sem(best_aic_mod6, data = dat_no.na, estimator="MLM")
# fitmeasures(best_aic.fit.mod6, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# modificationindices(best_aic.fit.mod6, sort. = TRUE)[1:15,]
# best_aic_mod7 <- paste0(best_aic_mod6, "\nsr_trans~tra_m")
# best_aic.fit.mod7 <- sem(best_aic_mod7, data = dat_no.na, estimator="MLM")
# fitmeasures(best_aic.fit.mod7, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# modificationindices(best_aic.fit.mod7, sort. = TRUE)[1:15,]
# best_aic_mod8 <- paste0(best_aic_mod7, "\nsub_trop_mbf~prs_m")
# best_aic.fit.mod8 <- sem(best_aic_mod8, data = dat_no.na, estimator="MLM")
# fitmeasures(best_aic.fit.mod8, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# # CFI down, RMSEA and AIC up

cat(best_aic_mod7)

png("../figures/DivRate_best_aic_model_2022.png", width=1000, height=1000, units = "px", res = 600)
title(round(fitmeasures(best_aic.fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
semPaths(object = best_aic.fit.mod7, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(best_aic.fit.mod7))
title(round(fitmeasures(best_aic.fit.mod7, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
dev.off()


### minimum AIC in good models, regardless of number of variables########################################
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

best_abs_aic_mod2 <- paste0(best_abs_aic_mod, "\nmdr~sub_trop_mbf")
best_abs_aic.fit.mod2 <- sem(best_abs_aic_mod2, data = dat_no.na, estimator="MLM")
fitmeasures(best_abs_aic.fit.mod2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_abs_aic.fit.mod2, sort. = TRUE)[1:15,]

best_abs_aic_mod3 <- paste0(best_abs_aic_mod2, "\nsr_trans  ~          prs_m")
best_abs_aic.fit.mod3 <- sem(best_abs_aic_mod3, data = dat_no.na, estimator="MLM")
fitmeasures(best_abs_aic.fit.mod3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))

modificationindices(best_abs_aic.fit.mod3, sort. = TRUE)[1:15,]
best_abs_aic_mod4 <- paste0(best_abs_aic_mod3, "\nsub_trop_mbf  ~  mio_pre_ano_m")
best_abs_aic.fit.mod4 <- sem(best_abs_aic_mod4, data = dat_no.na, estimator="MLM")
fitmeasures(best_abs_aic.fit.mod4, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))

modificationindices(best_abs_aic.fit.mod4, sort. = TRUE)[1:15,]
best_abs_aic_mod5 <- paste0(best_abs_aic_mod4, "\nmdr  ~          tra_m")
best_abs_aic.fit.mod5 <- sem(best_abs_aic_mod5, data = dat_no.na, estimator="MLM")
fitmeasures(best_abs_aic.fit.mod5, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_abs_aic.fit.mod5, sort. = TRUE)[1:15,]

best_abs_aic_mod6 <- paste0(best_abs_aic_mod5, "\nsr_trans  ~      medit_fws")
best_abs_aic.fit.mod6 <- sem(best_abs_aic_mod6, data = dat_no.na, estimator="MLM")
fitmeasures(best_abs_aic.fit.mod6, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_abs_aic.fit.mod6, sort. = TRUE)[1:15,]

best_abs_aic_mod7 <- paste0(best_abs_aic_mod6, "\nsr_trans ~ temp_cf")
best_abs_aic.fit.mod7 <- sem(best_abs_aic_mod7, data = dat_no.na, estimator="MLM")
fitmeasures(best_abs_aic.fit.mod7, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_abs_aic.fit.mod7, sort. = TRUE)[1:20,]

best_abs_aic_mod8 <- paste0(best_abs_aic_mod7, "\nsr_trans  ~  mio_pre_ano_m")
best_abs_aic.fit.mod8 <- sem(best_abs_aic_mod8, data = dat_no.na, estimator="MLM")
fitmeasures(best_abs_aic.fit.mod8, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_abs_aic.fit.mod8, sort. = TRUE)[1:20,]

best_abs_aic_mod9 <- paste0(best_abs_aic_mod8, "\nsr_trans  ~          tra_m")
best_abs_aic.fit.mod9 <- sem(best_abs_aic_mod9, data = dat_no.na, estimator="MLM")
fitmeasures(best_abs_aic.fit.mod9, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_abs_aic.fit.mod9, sort. = TRUE)[1:20,]
# STOP

cat(best_abs_aic_mod8)

model <- best_abs_aic.fit.mod8
png("../figures/DivRate_best_absolute_aic_model_2022.png", width=5000, height=5000, 
    res=600, units = "px")
semPaths(object = model, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(model))
title(round(fitmeasures(model, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
dev.off()





# Max variable number in good models ####################################################################################
(max.var <- max(fin_good$var.n))
fin_good_max.var <- fin_good[fin_good$var.n==max.var,]
nrow(fin_good_max.var) # good models with max var

cat(fin_good_max.var$model.name)
fin_good_max.var[,-7]
max.var.mod <- fin_good_max.var$model.name[which(!fin_good_max.var$multicolin_in_mdr)]

max.var.mod.fit <- sem(max.var.mod, data = dat_no.na, estimator="MLM")
summary(max.var.mod.fit, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
fitmeasures(max.var.mod.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit, sort. = TRUE)[1:15,]

max.var.mod_mod <- paste0(max.var.mod, "\nsub_trop_mbf~area")
max.var.mod.fit.mod <- sem(max.var.mod_mod, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit.mod, sort. = TRUE)[1:15,]

max.var.mod_mod2 <- paste0(max.var.mod_mod, "\n  mdr  ~   sub_trop_mbf")
max.var.mod.fit.mod2 <- sem(max.var.mod_mod2, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit.mod2, sort. = TRUE)[1:15,]

max.var.mod_mod3 <- paste0(max.var.mod_mod2, "\n sr_trans  ~prs_m")
max.var.mod.fit.mod3 <- sem(max.var.mod_mod3, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit.mod3, sort. = TRUE)[1:15,]

max.var.mod_mod4 <- paste0(max.var.mod_mod3, "\n sr_trans  ~      medit_fws")
max.var.mod.fit.mod4 <- sem(max.var.mod_mod4, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod4, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit.mod4, sort. = TRUE)[1:15,]

max.var.mod_mod5 <- paste0(max.var.mod_mod4, "\nsr_trans  ~        temp_cf")
max.var.mod.fit.mod5 <- sem(max.var.mod_mod5, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod5, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit.mod5, sort. = TRUE)[1:15,]

max.var.mod_mod6 <- paste0(max.var.mod_mod5, "\nsub_trop_mbf  ~  mio_pre_ano_m")
max.var.mod.fit.mod6 <- sem(max.var.mod_mod6, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod6, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# RMSE goes up, AIC down
modificationindices(max.var.mod.fit.mod6, sort. = TRUE)[1:15,]

max.var.mod_mod7 <- paste0(max.var.mod_mod6, "\nsub_trop_mbf  ~        temp_cf")
max.var.mod.fit.mod7 <- sem(max.var.mod_mod7, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod7, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit.mod7, sort. = TRUE)[1:15,]

max.var.mod_mod8 <- paste0(max.var.mod_mod7, "\nsr_trans  ~  mio_pre_ano_m")
max.var.mod.fit.mod8 <- sem(max.var.mod_mod8, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod8, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# RMSE goes up, AIC down
modificationindices(max.var.mod.fit.mod8, sort. = TRUE)[1:15,]

max.var.mod_mod9 <- paste0(max.var.mod_mod8, "\nsr_trans  ~ boreal_f_taiga")
max.var.mod.fit.mod9 <- sem(max.var.mod_mod9, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod9, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# stop
# modificationindices(max.var.mod.fit.mod9, sort. = TRUE)[1:15,]
# 
# max.var.mod_mod10 <- paste0(max.var.mod_mod9, "\nsub_trop_mbf~prs_m")
# max.var.mod.fit.mod10 <- sem(max.var.mod_mod10, data = dat_no.na, estimator="MLM")
# fitmeasures(max.var.mod.fit.mod10, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# # RMSEA jumps up

cat(max.var.mod_mod8)
model <- max.var.mod.fit.mod8
png("../figures/DivRate_max_var_model_2022.png", width=5000, height=5000, 
    res=600, units = "px")
semPaths(object = model, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(model))
title(round(fitmeasures(model, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
dev.off()


save(best_abs_aic_mod8, best_abs_aic.fit.mod8, file="../processed_data/best_model_DivRate.RData")
save.image("../processed_data/all_models_DivRate.RData")



# Save DivRate fit measures -----------------------------------------------

pres2 <- standardizedSolution(best_abs_aic.fit.mod8)
# filter for direct and indirect paths
pres2 <- pres2[pres2$op=="~" |pres2$op==":=" ,-grep("z|ci", names(pres2))]
# round estimates and p-values separately
pres2[,grep("est|se", names(pres2))] <- round(pres2[,grep("est|se", names(pres2))],2)
pres2[,grep("pvalue", names(pres2))] <- round(pres2[,grep("pvalue", names(pres2))],3)
temp <- pres2[order(pres2$op, pres2$lhs, abs(pres2$est.std), decreasing = T),]

write.csv(temp,"../processed_data/sem_model_output_DivRate.csv") 



# Interaction effects (piecewiseSEM) ----------------------
rm(list=ls())
library(piecewiseSEM)
load("../processed_data/best_model_DivRate.RData")
dat_no.na <- readRDS("../processed_data/sem_input_data.rds")

cat(best_abs_aic_mod8)

lm1 <- lm(sr_trans ~ soil + sub_trop_mbf + area + mdr + mat_m + medit_fws + 
            temp_cf + mio_pre_ano_m + prs_m, dat_no.na)
lm2 <- lm(mdr ~ mat_m + mio_mat_ano_m + boreal_f_taiga + prs_m+soil+temp_cf+
            area+medit_fws+mio_pre_ano_m+sub_trop_mbf+tra_m, dat_no.na)
lm3 <- lm(soil ~ area + mont_gs + sub_trop_mbf, dat_no.na)
lm4 <- lm(sub_trop_mbf ~ pre_m + tra_m + mat_m + area + mio_pre_ano_m, dat_no.na)

pmodel <- psem(lm1, lm2, lm3, lm4)
summary(pmodel, .progressBar = F) 

lm1.int <- lm(sr_trans ~ soil*area + sub_trop_mbf*area + area + mdr*area + mat_m*area + medit_fws*area + 
                temp_cf*area + mio_pre_ano_m*area + prs_m*area, dat_no.na)
lm2.int <- lm(mdr ~ mat_m*area + mio_mat_ano_m*area + boreal_f_taiga*area + prs_m*area + 
                soil*area + temp_cf*area + area + medit_fws*area + mio_pre_ano_m*area +
                sub_trop_mbf*area + tra_m*area, dat_no.na)
lm3.int <- lm(soil ~ area + mont_gs*area + sub_trop_mbf*area, dat_no.na)
lm4.int <- lm(sub_trop_mbf ~ pre_m*area + tra_m*area + mat_m*area + area + mio_pre_ano_m*area, dat_no.na)

pmodel.int <- psem(lm1.int, lm2.int, lm3.int, lm4.int)
summary(pmodel.int, .progressBar = F) 

rsquared(pmodel); rsquared(pmodel.int)
infCrit(pmodel); infCrit(pmodel.int)
plot(pmodel, digits=2)
plot(pmodel.int, digits=2)

fisherC(pmodel)


temp <- summary(pmodel.int)
temp <- temp$coefficients
names(temp)[ncol(temp)] <- "sig"
temp <- temp[order(temp$Response, abs(temp$Std.Estimate), decreasing = T),]
temp$Std.Estimate <- round(temp$Std.Estimate,2)
temp$Std.Error <- round(temp$Std.Error,2)
temp$P.Value <- round(temp$P.Value,3)

temp <- temp[,c("Response", "Predictor", "Std.Estimate", "Std.Error", "P.Value", "sig")]

write.csv(temp, "../processed_data/sem_model_coefficients_area_interactions_Div.csv")

# Interaction plots -------------------------------------------------------

#range(dat_no.na$area)
#seq(range(dat_no.na$area)[1],range(dat_no.na$area)[2])
#dat_no.na$area_g <- cut_interval(dat_no.na$area, n = 8)
#table(dat_no.na$area_g)
plot_grid(ncol=2,
ggplot(dat_no.na, aes(x=mat_m, y=sr_trans, group=area<0, col=area<0))+
  geom_point()+ 
  geom_smooth(method="lm")
,
ggplot(dat_no.na, aes(x=prs_m, y=sr_trans, group=area<0, col=area<0))+
  geom_point()+ 
  geom_smooth(method="lm")
)


# Figure  -----------------------------------------------------------------

tmp = 
grViz("
digraph SEM {

  # Node statements
//  subgraph clustercenter{
    node [shape = box, fontname = Helvetica, fontsize=10, style=filled, 
    color=white, peripheries=2, fillcolor=grey90]
      SR[label='SR \n R&#xb2;=0.65']
      RD [label='RD \n R&#xb2;=0.56']
//  }
//  subgraph clusterclimate {
    node [shape = ellipse, style=filled, fillcolor = '#71A8F4', 
    color=white, peripheries=2, margin=0.01, fontsize=8]
    mat[label='mean annual\n temperature']
    pre[label='total annual\n precipitation']
//  }
//  subgraph clusterbio {
    node [style=filled, fillcolor = '#9DF471', 
    color=white, peripheries=2]
    sub_trop_mbf[label='(sub)tropical moist\nbroadleaf forest']
    medit_fws[label='medit. forests, \n woodlands, scrub']
    temp_cf[label='temperate\nconiferous forest']
    boreal[label='boreal \n forests/taiga']
//  }
//  subgraph clusterseason {
    node [style=filled, fillcolor = '#F4F471', 
    color=white, peripheries=2]
      prs[label='precipitation \n seasonality']
      tra[label='annual \n temperature range']
//  }
//  subgraph clusterEH {
    node [style=filled, fillcolor = '#F471B2', 
    color=white, peripheries=2]
    soil[label='n soil types']
    area
//  }
    node [style=filled, fillcolor = '#F48771', 
    color=white, peripheries=2]
    mat_ano[label='Miocene \n temperature anomaly']
    pre_ano[label='Miocene \n precipitation anomaly']

  # Edge statements
  
  edge [color = red, penwidth=1, len=1]
    prs -> SR
    {prs; sub_trop_mbf; pre_ano} -> RD
  
  edge [color = red, penwidth=1, style=dashed, len=1]
    {tra; area} -> RD
    tra -> sub_trop_mbf
  
  edge [color = red, penwidth=3, style=normal, len=1, arrowsize=0.7]
    mat -> RD
    
  edge [color = black, penwidth=1, style=normal, len=1, arrowsize=1]
    {mat; area; medit_fws; temp_cf} -> SR
    {mat_ano; boreal; soil}  -> RD
    {area mat}  -> sub_trop_mbf
  
  edge [color = black, penwidth=1, style=dashed, minlen=1]
    {pre_ano; RD} -> SR
    {temp_cf; medit_fws} -> RD
  
  edge [color = black, penwidth=3, style=normal, len=1, arrowsize=0.7]
    sub_trop_mbf -> SR
    area -> soil
  
  edge [color = black, penwidth=5, style=normal, arrowsize=0.4]
    soil -> SR
    pre -> sub_trop_mbf[weight=0]
    
  edge [style=invis]
//    prs -> {RD}
//    prs -> SR
    pre -> SR
    tra -> prs

//  graph[nodesep=0.1, layout=dot, ranksep=.3, rankdir=BT, splines=true]
//  graph[nodesep=0.01, layout=circo, mindist=.7]  
  graph[layout=twopi, ranksep=1, root=SR, splines=true, overlap=false]
    
}
")
## 2. Convert to SVG, then save as png ####
tmp = DiagrammeRsvg::export_svg(tmp)
tmp = charToRaw(tmp) # flatten
rsvg::rsvg_svg(tmp, "../figures/sem_DivRates.svg")
