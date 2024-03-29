# get stats from model selection runs, manually adjust models
# produces models.RData



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


# Load data for MRD -----------------------------------------------------------

rm(list = setdiff(ls(), lsf.str())) 
filnames <- dir("../processed_data/sem_output/Wednesday", full.names = T)
filnames <- filnames[grep("results_", filnames)]
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
models[grep("NA", model.name)] # should be zero

fit2$model.name <- unlist(model.name)
rm(model.name)

# model fit
length(which(fit2$cfi.robust>0.9)) # acceptable CFI
length(which(fit2$rmsea.robust<0.11)) # RMSEA < 0.11

# explanatory power for SR and MRD
pow <- sapply(models, "[", 3)
rm(models)
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
#hist(gm.sr.n, xlab="variables in SR regression") # the number of variables in the good models
fit2$gm.sr.n <- gm.sr.n
rm(gm.sr.sp)

mrd1 <- sub(".*\nmrd ~ ", "", fit2$model.name)
gm.mrd <- sub("\nsoil.*", "", mrd1)
gm.mrd.sp <- strsplit(gm.mrd, "+", fixed = TRUE)
gm.mrd.n <- unlist(lapply(gm.mrd.sp, length))
#hist(gm.mrd.n, xlab="variables in MRD regression") # the number of variables in the good models
fit2$gm.mrd.n <- gm.mrd.n

# check if pet_m and mat_m appear together in MRD
fit2$mat_in_mrd <- grepl("mat_m", gm.mrd.sp)
fit2$pet_in_mrd <- grepl("pet_m", gm.mrd.sp)
fit2$multicolin_in_mrd <- (fit2$mat_in_mrd==TRUE & fit2$pet_in_mrd==TRUE)
table(fit2$multicolin_in_mrd)
# is all FALSE since removed pet due do colinerarity

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
# fitmeasures(best_cfi.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_cfi.fit, sort. = TRUE)[1:15,]

# MODS on the best model according to CFI #
fitmeasures(best_cfi.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
best_cfi.mod <-  paste0(best_cfi, "\nsub_trop_mbf ~ area")
best_cfi.mod.fit <- sem(best_cfi.mod, data = dat_no.na, estimator="MLM")
fitmeasures(best_cfi.mod.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_cfi.mod.fit, sort. = TRUE)[1:15,]

best_cfi.mod2 <-  paste0(best_cfi.mod, "\nmrd ~ area")
best_cfi.fit.mod2 <- sem(best_cfi.mod2, data = dat_no.na, estimator="MLM")
fitmeasures(best_cfi.fit.mod2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_cfi.fit.mod2, sort. = TRUE)[1:15,]

best_cfi.mod3 <-  paste0(best_cfi.mod2, "\nsr_trans ~ prs_m")
best_cfi.fit.mod3 <- sem(best_cfi.mod3, data = dat_no.na, estimator="MLM")
fitmeasures(best_cfi.fit.mod3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_cfi.fit.mod3, sort. = TRUE)[1:15,]

best_cfi.mod4 <-  paste0(best_cfi.mod3, "\nsr_trans ~ tra_m")
best_cfi.fit.mod4 <- sem(best_cfi.mod4, data = dat_no.na, estimator="MLM")
fitmeasures(best_cfi.fit.mod4, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_cfi.fit.mod4, sort. = TRUE)[1:20,]
# RMSEA + AIC slightly up

best_cfi.mod5 <-  paste0(best_cfi.mod4, "\nsub_trop_mbf ~soil")
best_cfi.fit.mod5 <- sem(best_cfi.mod5, data = dat_no.na, estimator="MLM")
fitmeasures(best_cfi.fit.mod5, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# rmsea goes up
modificationindices(best_cfi.fit.mod5, sort. = TRUE)[1:20,]

best_cfi.mod6 <-  paste0(best_cfi.mod5, "\nsub_trop_mbf ~prs_m")
best_cfi.fit.mod6 <- sem(best_cfi.mod6, data = dat_no.na, estimator="MLM")
fitmeasures(best_cfi.fit.mod6, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# stop
# modificationindices(best_cfi.fit.mod6, sort. = TRUE)[1:20,]
# 
# best_cfi.fit.mod7 <- update(best_cfi.fit.mod6, add="mrd ~ tra_m")
# fitmeasures(best_cfi.fit.mod7, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# # RMSEA goes up

png("../figures/best_cfi_model_2022.png", width=500, height=500, units = "px")
#par(mfrow=c(1,2))
# semPaths(object = best_cfi.fit, layout = "circle2", rotation = 1,
#          whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE,
#          edgeLabels = sem_sig_labels(best_cfi.fit))
# title(round(fitmeasures(best_cfi.fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
model <- best_cfi.fit.mod4
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
best_rmsea_mod <- paste0(best_rmsea, "\nsub_trop_mbf  ~ area") #mrd~area
best_rmsea.fit.mod <- sem(best_rmsea_mod, data = dat_no.na, estimator="MLM")
fitmeasures(best_rmsea.fit.mod, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_rmsea.fit.mod, sort. = TRUE)[1:15,]
best_rmsea_mod2 <- paste0(best_rmsea_mod, "\nmrd  ~ area")
best_rmsea.fit.mod2 <- sem(best_rmsea_mod2, data = dat_no.na, estimator="MLM")
fitmeasures(best_rmsea.fit.mod2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_rmsea.fit.mod2, sort. = TRUE)[1:15,]
best_rmsea_mod3 <- paste0(best_rmsea_mod2, "\nsr_trans ~ prs_m")
best_rmsea.fit.mod3 <- sem(best_rmsea_mod3, data = dat_no.na, estimator="MLM")
fitmeasures(best_rmsea.fit.mod3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_rmsea.fit.mod3, sort. = TRUE)[1:15,]
best_rmsea_mod4 <- paste0(best_rmsea_mod3, "\nsr_trans~tra_m")
best_rmsea.fit.mod4 <- sem(best_rmsea_mod4, data = dat_no.na, estimator="MLM")
fitmeasures(best_rmsea.fit.mod4, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_rmsea.fit.mod4, sort. = TRUE)[1:20,]
best_rmsea_mod5 <- paste0(best_rmsea_mod4, "\nsub_trop_mbf ~ soil")
best_rmsea.fit.mod5 <- sem(best_rmsea_mod5, data = dat_no.na, estimator="MLM")
fitmeasures(best_rmsea.fit.mod5, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_rmsea.fit.mod5, sort. = TRUE)[1:20,]
best_rmsea_mod6 <- paste0(best_rmsea_mod5, "\nsub_trop_mbf ~ prs_m")
best_rmsea.fit.mod6 <- sem(best_rmsea_mod6, data = dat_no.na, estimator="MLM")
fitmeasures(best_rmsea.fit.mod6, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# stop

# modificationindices(best_rmsea.fit.mod6, sort. = TRUE)[1:25,]
# best_rmsea_mod7 <- paste0(best_rmsea_mod6, "\n")
# best_rmsea.fit.mod7 <- sem(best_rmsea_mod7, data = dat_no.na, estimator="MLM")
# fitmeasures(best_rmsea.fit.mod7, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# #AIC going up
# modificationindices(best_rmsea.fit.mod7, sort. = TRUE)[1:20,]
# # stop, nothing meaninful to add


cat(best_rmsea_mod3)

png("../figures/best_rmsea_model2022.png", width=500, height=500, units = "px")
# par(mfrow=c(1,2))
# semPaths(object = best_rmsea.fit, layout = "circle2", rotation = 1,
#          whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE,
#          edgeLabels = sem_sig_labels(best_rmsea.fit))
# title(round(fitmeasures(best_rmsea.fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
model <- best_rmsea.fit.mod3
semPaths(object = model, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(model))
title(round(fitmeasures(model, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
dev.off()




### 3) GLOBAL MIN AIC, regardless of model fit ###################################
cat(fit2$model.name[which.min(fit2$aic)])
fit2[which.min(fit2$aic),-7]
best_global_aic <- fit2$model.name[which.min(fit2$aic)]
best_global_aic.fit <- sem(best_global_aic, data = dat_no.na, estimator="MLM")
#summary(best_global_aic.fit, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)

fitmeasures(best_global_aic.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_global_aic.fit, sort. = TRUE)[1:15,]
best_global_aic.mod <- paste0(best_global_aic, "\nsub_trop_mbf  ~ temp_bmf")
best_global_aic.mod.fit <- sem(best_global_aic.mod, data = dat_no.na, estimator="MLM")
fitmeasures(best_global_aic.mod.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_global_aic.mod.fit, sort. = TRUE)[1:15,]
best_global_aic.mod2 <- paste0(best_global_aic.mod, "\nsr_trans  ~ temp_bmf")
best_global_aic.mod.fit2 <- sem(best_global_aic.mod2, data = dat_no.na, estimator="MLM")
fitmeasures(best_global_aic.mod.fit2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_global_aic.mod.fit2, sort. = TRUE)[1:20,]
best_global_aic.mod3 <- paste0(best_global_aic.mod2, "\nmrd~area")
best_global_aic.mod.fit3 <- sem(best_global_aic.mod3, data = dat_no.na, estimator="MLM")
fitmeasures(best_global_aic.mod.fit3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_global_aic.mod.fit3, sort. = TRUE)[1:25,]
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
# # par(mfrow=c(1,2))
# # semPaths(object = best_global_aic.fit, layout = "circle2", rotation = 1,
# #          whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
# #          exoCov = FALSE, exoVar = FALSE,
# #          nCharNodes = 0, fade=FALSE,
# #          edgeLabels = sem_sig_labels(best_global_aic.fit))
# # title(round(fitmeasures(best_global_aic.fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
# model <- 
# semPaths(object = model, layout = "circle2", rotation = 1,
#          whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE,
#          edgeLabels = sem_sig_labels(model))
# title(round(fitmeasures(model, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
# dev.off()


## 4) pre-select good CFI and RMSEA#####################################################################
# More variables while still having acceptable model fit is desirable since 
dev.off()
plot(fit2$rmsea.robust[fit2$rmsea.robust<=0.15], fit2$var.n[fit2$rmsea.robust<=0.15])
fin_good <- fit2[fit2$cfi.robust>0.9 & fit2$rmsea.robust<0.11,] # leaving some room for improvement RMSEA, and to capture the best cfi

nrow(fin_good) # models with potential


## pre-select, best AIC ranking, most variables ######################
min.aic <- min(fin_good$aic)
fin_good_aic <- fin_good[fin_good$aic<min.aic+2,]
nrow(fin_good_aic) # good models with lowest AIC + max 2 more aic points = set of basically identical AICs
fin_good_aic[,-7]
### most variables
hist(fin_good_aic$var.n)

(m.aic <- max(fin_good_aic$var.n)) # 12 variables
fin_good_aic_min <- fin_good_aic[which(fin_good_aic$var.n==m.aic),] 
nrow(fin_good_aic_min) # 6 good models
cat(fin_good_aic_min$model.name)
fin_good_aic_min[which.min(fin_good_aic_min$aic),-7]

# select the one with lowest AIC
best_aic <- fin_good_aic_min$model.name[which.min(fin_good_aic_min$aic)]
best_aic.fit <- sem(best_aic, data = dat_no.na, estimator="MLM")
#summary(best_aic.fit, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
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
best_aic_mod3 <- paste0(best_aic_mod2, "\nsr_trans~prs_m")
best_aic.fit.mod3 <- sem(best_aic_mod3, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_aic.fit.mod3, sort. = TRUE)[1:15,]
best_aic_mod4 <- paste0(best_aic_mod3, "\nsub_trop_mbf~tri")
best_aic.fit.mod4 <- sem(best_aic_mod4, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod4, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_aic.fit.mod4, sort. = TRUE)[1:15,]
best_aic_mod5 <- paste0(best_aic_mod4, "\nsr_trans~tri") 
best_aic.fit.mod5 <- sem(best_aic_mod5, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod5, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_aic.fit.mod5, sort. = TRUE)[1:15,]
best_aic_mod6 <- paste0(best_aic_mod5, "\nsr_trans~tra_m")
best_aic.fit.mod6 <- sem(best_aic_mod6, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod6, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_aic.fit.mod6, sort. = TRUE)[1:15,]
best_aic_mod7 <- paste0(best_aic_mod6, "\nsub_trop_mbf~prs_m")
best_aic.fit.mod7 <- sem(best_aic_mod7, data = dat_no.na, estimator="MLM")
fitmeasures(best_aic.fit.mod7, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# stop
# modificationindices(best_aic.fit.mod7, sort. = TRUE)[1:15,]
# best_aic_mod8 <- paste0(best_aic_mod7, "\nsoil~sub_trop_mbf")
# best_aic.fit.mod8 <- sem(best_aic_mod8, data = dat_no.na, estimator="MLM")
# fitmeasures(best_aic.fit.mod8, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# # # CFI up, RMSEA up, AIC down
# modificationindices(best_aic.fit.mod8, sort. = TRUE)[1:15,]
# best_aic_mod9 <- paste0(best_aic_mod8, "\nsub_trop_mbf~prs_m")
# best_aic.fit.mod9 <- sem(best_aic_mod9, data = dat_no.na, estimator="MLM")
# fitmeasures(best_aic.fit.mod9, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))


cat(best_aic_mod6)

png("../figures/best_aic_model_2022.png", width=500, height=500, units = "px")
# par(mfrow=c(1,2))
# semPaths(object = best_aic.fit, layout = "circle2", rotation = 1,
#          whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE,
#          edgeLabels = sem_sig_labels(best_aic.fit))
# title(round(fitmeasures(best_aic.fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
model = best_aic.fit.mod6
semPaths(object = model, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(model))
title(round(fitmeasures(model, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
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

best_abs_aic_mod2 <- paste0(best_abs_aic_mod, "\nmrd~area")
best_abs_aic.fit.mod2 <- sem(best_abs_aic_mod2, data = dat_no.na, estimator="MLM")
fitmeasures(best_abs_aic.fit.mod2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_abs_aic.fit.mod2, sort. = TRUE)[1:15,]

best_abs_aic_mod3 <- paste0(best_abs_aic_mod2, "\nsr_trans ~ prs_m")
best_abs_aic.fit.mod3 <- sem(best_abs_aic_mod3, data = dat_no.na, estimator="MLM")
fitmeasures(best_abs_aic.fit.mod3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))

modificationindices(best_abs_aic.fit.mod3, sort. = TRUE)[1:15,]
best_abs_aic_mod4 <- paste0(best_abs_aic_mod3, "\nsub_trop_mbf ~ tri")
best_abs_aic.fit.mod4 <- sem(best_abs_aic_mod4, data = dat_no.na, estimator="MLM")
fitmeasures(best_abs_aic.fit.mod4, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))

modificationindices(best_abs_aic.fit.mod4, sort. = TRUE)[1:15,]
best_abs_aic_mod5 <- paste0(best_abs_aic_mod4, "\nsr_trans ~ tri")
best_abs_aic.fit.mod5 <- sem(best_abs_aic_mod5, data = dat_no.na, estimator="MLM")
fitmeasures(best_abs_aic.fit.mod5, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_abs_aic.fit.mod5, sort. = TRUE)[1:15,]

best_abs_aic_mod6 <- paste0(best_abs_aic_mod5, "\nsr_trans ~ tra_m")
best_abs_aic.fit.mod6 <- sem(best_abs_aic_mod6, data = dat_no.na, estimator="MLM")
fitmeasures(best_abs_aic.fit.mod6, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(best_abs_aic.fit.mod6, sort. = TRUE)[1:15,]

best_abs_aic_mod7 <- paste0(best_abs_aic_mod6, "\nsub_trop_mbf  ~prs_m")
best_abs_aic.fit.mod7 <- sem(best_abs_aic_mod7, data = dat_no.na, estimator="MLM")
fitmeasures(best_abs_aic.fit.mod7, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# STOP
# modificationindices(best_abs_aic.fit.mod7, sort. = TRUE)[1:15,]
# 
# best_abs_aic_mod8 <- paste0(best_abs_aic_mod7, "\nsub_trop_mbf ~ prs_m")
# best_abs_aic.fit.mod8 <- sem(best_abs_aic_mod8, data = dat_no.na, estimator="MLM")
# fitmeasures(best_abs_aic.fit.mod8, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# # rmsea jumps up, stop
# 

cat(best_abs_aic_mod6)

png("../figures/best_absolute_aic_model_2022.png", width=500, height=500, units = "px")
# par(mfrow=c(1,2))
# semPaths(object = best_abs_aic.fit, layout = "circle2", rotation = 1,
#          whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE,
#          edgeLabels = sem_sig_labels(best_abs_aic.fit))
# title(round(fitmeasures(best_abs_aic.fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
model = best_abs_aic.fit.mod6
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
# multicolineariy? no. take the one with best fit
fin_good_max.var[,-7]
max.var.mod <- fin_good_max.var$model.name[which.min(fin_good_max.var$rmsea.robust)]

max.var.mod.fit <- sem(max.var.mod, data = dat_no.na, estimator="MLM")
summary(max.var.mod.fit, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
fitmeasures(max.var.mod.fit, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit, sort. = TRUE)[1:15,]

max.var.mod_mod <- paste0(max.var.mod, "\nsub_trop_mbf~area")
max.var.mod.fit.mod <- sem(max.var.mod_mod, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit.mod, sort. = TRUE)[1:15,]

max.var.mod_mod2 <- paste0(max.var.mod_mod, "\n soil  ~  mont_gs")
max.var.mod.fit.mod2 <- sem(max.var.mod_mod2, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod2, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit.mod2, sort. = TRUE)[1:15,]

max.var.mod_mod3 <- paste0(max.var.mod_mod2, "\n mrd~area")
max.var.mod.fit.mod3 <- sem(max.var.mod_mod3, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod3, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit.mod3, sort. = TRUE)[1:15,]

max.var.mod_mod4 <- paste0(max.var.mod_mod3, "\nsr_trans  ~ prs_m")
max.var.mod.fit.mod4 <- sem(max.var.mod_mod4, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod4, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit.mod4, sort. = TRUE)[1:15,]

max.var.mod_mod5 <- paste0(max.var.mod_mod4, "\nsub_trop_mbf~tri")
max.var.mod.fit.mod5 <- sem(max.var.mod_mod5, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod5, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit.mod5, sort. = TRUE)[1:15,]

max.var.mod_mod6 <- paste0(max.var.mod_mod5, "\nsr_trans~tri")
max.var.mod.fit.mod6 <- sem(max.var.mod_mod6, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod6, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit.mod6, sort. = TRUE)[1:15,]

max.var.mod_mod7 <- paste0(max.var.mod_mod6, "\nsr_trans~tra_m") 
max.var.mod.fit.mod7 <- sem(max.var.mod_mod7, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod7, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(max.var.mod.fit.mod7, sort. = TRUE)[1:15,]

max.var.mod_mod8 <- paste0(max.var.mod_mod7, "\nsub_trop_mbf~prs_m")
max.var.mod.fit.mod8 <- sem(max.var.mod_mod8, data = dat_no.na, estimator="MLM")
fitmeasures(max.var.mod.fit.mod8, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# STOP
# modificationindices(max.var.mod.fit.mod8, sort. = TRUE)[1:15,]
# 
# max.var.mod_mod9 <- paste0(max.var.mod_mod8, "\nsr_trans~tri")
# max.var.mod.fit.mod9 <- sem(max.var.mod_mod9, data = dat_no.na, estimator="MLM")
# fitmeasures(max.var.mod.fit.mod9, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# 
# modificationindices(max.var.mod.fit.mod9, sort. = TRUE)[1:15,]
# 
# max.var.mod_mod10 <- paste0(max.var.mod_mod9, "\nsub_trop_mbf~prs_m")
# max.var.mod.fit.mod10 <- sem(max.var.mod_mod10, data = dat_no.na, estimator="MLM")
# fitmeasures(max.var.mod.fit.mod10, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# # RMSEA jumps up

cat(max.var.mod_mod7)

png("../figures/max_var_model_2022.png", width=5000, height=5000, res=600, units = "px")
#par(mfrow=c(1,2))
# semPaths(object = max.var.mod.fit, layout = "circle2", rotation = 1,
#          whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE,
#          edgeLabels = sem_sig_labels(max.var.mod.fit))
# title(round(fitmeasures(max.var.mod.fit, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
semPaths(object = max.var.mod.fit.mod7, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(max.var.mod.fit.mod7))
title(round(fitmeasures(max.var.mod.fit.mod7, c("cfi.robust", "rmsea.robust", "aic")),3), cex.main=1)
dev.off()

cat(max.var.mod_mod7)
save(max.var.mod_mod7, max.var.mod.fit.mod7, file="../processed_data/best_model.RData")
save.image("../processed_data/all_models.RData")



# Interaction effects (piecewiseSEM) ----------------------
rm(list=ls())
library(piecewiseSEM)
library(lavaan)
load("../processed_data/best_model.RData")
cat(max.var.mod_mod7)
dat_no.na <- readRDS("../processed_data/sem_input_data.rds")
lm1 <- lm(sr_trans ~ soil + sub_trop_mbf + area + mrd + mont_gs + mat_m + prs_m + 
            tri + tra_m, dat_no.na)
lm2 <- lm(mrd ~ pre_m + sub_trop_mbf + prs_m + mio_mat_ano_m+soil+tri+area, dat_no.na)
lm3 <- lm(soil ~ area + mont_gs, dat_no.na)
lm4 <- lm(sub_trop_mbf ~ pre_m + tra_m + mat_m + area + tri, dat_no.na)

pmodel <- psem(lm1, lm2, lm3, lm4)
summary(pmodel, .progressBar = F) 

lm1.int <- lm(sr_trans ~ soil*area + sub_trop_mbf*area + area + mrd*area + mont_gs*area+
                mat_m*area + prs_m*area + tra_m*area + tri*area, dat_no.na)
lm2.int <- lm(mrd ~ pre_m*area + sub_trop_mbf*area + prs_m*area + mio_mat_ano_m*area+
                soil*area+tri*area+area, dat_no.na)
lm3.int <- lm(soil ~ area + mont_gs*area, dat_no.na)
lm4.int <- lm(sub_trop_mbf ~ pre_m*area + tra_m*area + mat_m*area + area + tri*area, dat_no.na)

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

write.csv(temp, "../processed_data/sem/sem_model_coefficients_ares_interactions.csv")


# Interaction plots -------------------------------------------------------

#range(dat_no.na$area)
#seq(range(dat_no.na$area)[1],range(dat_no.na$area)[2])
#dat_no.na$area_g <- cut_interval(dat_no.na$area, n = 8)
#table(dat_no.na$area_g)
plot_grid(ncol=3,
          ggplot(dat_no.na, aes(x=prs_m, y=sr_trans, group=area<1, col=area<1))+
            geom_point()+ 
            geom_smooth(method="lm")+
            theme(legend.position = c(0.85,0.85)),
          ggplot(dat_no.na, aes(x=mat_m, y=sr_trans, group=area<0, col=area<0))+
            geom_point()+ 
            geom_smooth(method="lm")+
            theme(legend.position = c(0.15,0.85)),
          ggplot(dat_no.na, aes(x=tri, y=sr_trans, group=area<0, col=area<0))+
            geom_point()+ 
            geom_smooth(method="lm")+
            theme(legend.position = c(0.85,0.85))
)
plot_grid(ncol=2,
          ggplot(dat_no.na, aes(x=pre_m, y=mrd, group=area<0, col=area<0))+
            geom_point()+ 
            geom_smooth(method="lm")+
            theme(legend.position = c(0.85,0.85)),
          ggplot(dat_no.na, aes(x=prs_m, y=mrd, group=area<1, col=area<1))+
            geom_point()+ 
            geom_smooth(method="lm")+
            theme(legend.position = c(0.85,0.85))
)



# Islands  --------------------------------------------------------------
## cite http://dx.doi.org/10.1111/j.1466-8238.2011.00728.x for data 
islands2 <- read.csv("../data/Islands_TDWG_AllData.txt")
dim(islands2)
table(islands2$GeologicalOrigin)
# chose atoll + volcanic
islands2$origins[which(islands2$GeologicalOrigin %in% c("atoll", "volcanic"))] <- "island"
islands2$origins[which(!islands2$GeologicalOrigin %in% c("atoll", "volcanic"))] <- "other"
dat <- dat_no.na
dat$LEVEL_3_CO <- row.names(dat)
dat <- merge(dat, islands2[,c("LEVEL_3_CO", "origins")], all.x=TRUE)
dat$origins[is.na(dat$origins)] <- "other"

regular <- sem(mod, data = dat, estimator="MLM")
multigroup1.volcanic <- sem(mod, data = dat, estimator="MLM", group = "origins")
# error: no variance in mont_gs
# solution: add tiny variation
dat$r.mont_gs <- dat$mont_gs + 10
dat$r.mont_gs <- rnorm(length(dat$r.mont_gs), mean = dat$r.mont_gs, sd = 0.0001*dat$r.mont_gs)
dat$r.mont_gs <- dat$r.mont_gs-10
range(dat$r.mont_gs-dat$mont_gs)
range(dat$mont_gs)
dat$mont_gs <- dat$r.mont_gs

multigroup1 <- sem(mod, data = dat, estimator="MLM", group = "origins")
#summary(multigroup1, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
fitmeasures(multigroup1, c("cfi.robust", "rmsea.robust", "aic", "pvalue.scaled", "chisq.scaled", "df", "chisq.scaling.factor"))

multigroup1.constrained <- sem(mod, dat, estimator="MLM", group = "origins", group.equal = c("intercepts", "regressions"))
#summary(multigroup1.constrained, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
fitmeasures(multigroup1, c("cfi.robust", "rmsea.robust", "aic", "pvalue.scaled", "chisq.scaled", "df", "chisq.scaling.factor"))
fitmeasures(regular, c("cfi.robust", "rmsea.robust", "aic", "pvalue.scaled", "chisq.scaled", "df", "chisq.scaling.factor"))

write.csv(standardizedsolution(multigroup1), "../processed_data/sem/island_model_SEM_output.csv")

anova(multigroup1, multigroup1.constrained)
# the free and contrained model are significantly different: some paths vary between groups while others do not
# Example: let prs_m vary between groups

multigroup2 <- sem(mod2, dat, estimator="MLM", group = "origins")
summary(multigroup2, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
anova(multigroup1, multigroup2)
# sig not diff; global contraints are valid, this path does not vary between groups

pmodel <- psem(lm1, lm2, lm3, lm4)
(pmultigroup_volcanic <- multigroup(pmodel, group = "origins"))
# notice the constrained to the global model paths
# diffs between groups:
## sr_trans~area
## sr_trans~prs_m
## mrd~pre_m
## soil~area
## sub_trop_mbf~pre_m
## sub_trop_mbf~mat_m
## sub_trop_mbf~areA

## all islands: 
# notice the constrained to the global model paths
# diffs between groups:
## mrd~pre_m
## mrd_prs_m
## sub_trop_mbf~pre_m
## sub_trop_mbf~mat_m
## sub_trop_mbf~area
## sub_trop_mbf~tri

fitmeasures(multigroup1, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
fitmeasures(multigroup1.volcanic, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))

m1 <- standardizedsolution(multigroup1)
m1[,grep("est|se", names(m1))] <- round(m1[,grep("est|se", names(m1))],2)
m1[,grep("pvalue", names(m1))] <- round(m1[,grep("pvalue", names(m1))],3)
temp1 <- m1[order(m1$op, m1$group, m1$lhs, abs(m1$est.std), decreasing = T),]
write.csv(temp1, "../processed_data/sem/island_model_SEM_output.csv")
m1 <- standardizedsolution(multigroup1.volcanic)
m1[,grep("est|se", names(m1))] <- round(m1[,grep("est|se", names(m1))],2)
m1[,grep("pvalue", names(m1))] <- round(m1[,grep("pvalue", names(m1))],3)
temp2 <- m1[order(m1$op, m1$group, m1$lhs, abs(m1$est.std), decreasing = T),]
write.csv(temp2, "../processed_data/sem/island_volcanic_model_SEM_output.csv")

temp1$all_islands <- "all islands"
temp2$all_islands <- "volcanic islands"
temp <- merge(temp1, temp2, all=TRUE)
temp$group[temp$group==1] <- "mainland"
temp$group[temp$group==2] <- "island"
head(temp)
temp$path <- paste(temp$lhs, temp$op, temp$rhs)
ggplot(temp[temp$op=="~",], aes(x=path, y=est.std, col=all_islands, 
                                size=2, shape=factor(group)))+
  geom_point()+
  theme(axis.text.x = element_text(angle=45, hjust=1))+
  geom_hline(yintercept=0)

## all island models have worse fit, volcanic islands only has less worse fit
## check csvs for order

