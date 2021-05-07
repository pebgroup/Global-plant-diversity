rm(list=ls())
library(sf)
library(lavaan)
library(semPlot)
library(modEvA)
library(fitdistrplus)
source("scripts/clean_semplot_functions.R") # script for modified plotting functions

# load data ###################################################################################################
dat_no.na <-  readRDS("data/data_for_SEM.rds")
names(dat_no.na) <- sub("_mean$", "_m", names(dat_no.na)) # shorten mean to _m

# sqrt transform species richness for better normal distribution
dat_no.na$sr_trans <- sqrt(dat_no.na$sr)
#hist(dat_no.na$sr_trans)


# test distribution of potential response variables ############################################################
## distributions SR and MRD 
hist(dat_no.na[,c("sr", "mrd")])
hist(dat_no.na[,c("sr_trans", "mrd")])
hist(dat_no.na$sr, breaks=100)
# mrd is kinda normally distributed, SR not. Lets test:
descdist(dat_no.na$mrd, discrete = FALSE, boot=500)
shapiro.test(dat_no.na$mrd) # not far away from normal distribution
mrd_norm_fit <- fitdist(dat_no.na$mrd, "norm",method="mle")
summary(mrd_norm_fit)
plot(mrd_norm_fit)

descdist(dat_no.na$sr, discrete = TRUE, boot=500) # Poisson regression or negative binomial regression
descdist(log(dat_no.na$sr), discrete = FALSE, boot=500) # doesnt make it better
descdist(sqrt(dat_no.na$sr), discrete = FALSE, boot=500) # lognormal?
fitp <- fitdist(dat_no.na$sr,"pois") # untransformed poisson
fnbinom <- fitdist(dat_no.na$sr,"nbinom", discrete = TRUE) # untransformed neg.binomial
fnorm_nt <- fitdist(dat_no.na$sr,"norm") # untransformed normal
flnorm <- fitdist(sqrt(dat_no.na$sr),"lnorm")
fnorm <- fitdist(sqrt(dat_no.na$sr),"norm")

summary(fitp)
summary(fnbinom)
summary(fnorm_nt)
#_____________
summary(flnorm)
summary(fnorm)

plot(fitp) # numbers too large
plot(fnbinom)
plot(fnorm_nt)
#_____________
plot(flnorm)
plot(fnorm)

# test with soil + area
lmp <- lm(sr~soil + area, dat_no.na)
lms <- lm(sqrt(sr)~soil + area, dat_no.na)
m1 <- glm.nb(sr ~ soil + area, dat_no.na)
hist(lmp$residuals, breaks=20)
hist(lms$residuals, breaks=20)
hist(m1$residuals, breaks=20)
modEvA::RsqGLM(m1)

# distributions secondary response variables ##
## sub_trop_mbf is a zero inflated variable, however will be treated as the others since the linear regression captures the basic connection
descdist(dat_no.na$sub_trop_mbf, discrete = TRUE, boot=500)
hist(dat_no.na$sub_trop_mbf, breaks=30)
hist(log(dat_no.na$sub_trop_mbf), breaks=30)
plot(dat_no.na$sub_trop_mbf~dat_no.na$pre_m)
plot(dat_no.na$sub_trop_mbf~dat_no.na$tra_m)

# ## mont_gs: zero inflated too, but same reasoning as in sub_trop_mbf
# descdist(dat_no.na$mont_gs, discrete = TRUE, boot=500)
# hist(dat_no.na$mont_gs, breaks=30)
# plot(dat_no.na$mont_gs~dat_no.na$elev_m)
# abline(lm(dat_no.na$mont_gs~dat_no.na$elev_m))

## deserts_x_shrub: zero inflated
descdist(dat_no.na$deserts_x_shrub, discrete = TRUE, boot=500)
plot(dat_no.na$deserts_x_shrub~dat_no.na$pet_m)
abline(lm(dat_no.na$deserts_x_shrub~dat_no.na$pet_m))

## tra_sd
descdist(dat_no.na$tra_sd, discrete = TRUE, boot=500)
hist(dat_no.na$tra_sd, breaks=50)
hist(sqrt(dat_no.na$tra_sd), breaks=50)
shapiro.test(sqrt(dat_no.na$tra_sd))

## sea_sd
descdist(dat_no.na$sea_sd, discrete = TRUE, boot=500)
hist(log10(dat_no.na$sea_sd), breaks=50)
shapiro.test(log(dat_no.na$sea_sd))

## soil
descdist(dat_no.na$soil, discrete = TRUE, boot=500)
hist(dat_no.na$soil)
shapiro.test(dat_no.na$soil)

# transformations:
dat_no.na$tra_sd <- sqrt(dat_no.na$tra_sd)  
dat_no.na$sea_sd <- sqrt(dat_no.na$sea_sd)  

# scale data #################################################################################################
ztrans <- function(x){(x - mean(x)) / sd(x)} #z-transformation to avoid variance issues
dat_no.na <- apply(dat_no.na, 2, ztrans)
dat_no.na <- as.data.frame(dat_no.na)

saveRDS(dat_no.na, "data/sem_input_data_noferns.rds")










# interaction between precipitation and temperature #######################################################################
lm_simp <- lm(sub_trop_mbf ~ mat_m + pre_m, dat_no.na)
lm_int <- lm(sub_trop_mbf ~ mat_m + pre_m + mat_m*pre_m, dat_no.na)
summary(lm_simp)
summary(lm_int)
# model building description ###############################################################################################
# procedure: 
# 1) check variables for influence using GBM to avoid multicollinerarity as a first step of analysis
# 2) Check most likely regressions for potential vif() issues
# 3) Build SEMs following brute force selection, including the most important variables
# 4) check final model again

# estimator and bootstrap comparison ###########################################################################
tsem_full <- "
          # regressions
          sr_trans ~ soil + sub_trop_mbf + sub_trop_dbf+ area+ pre_sd+ mrd+ mont_gs + pre_m+ sea_sd+ tra_m
          +  mat_m +tra_sd+ mat_sd + sea_m + medit_fws + temp_gss + tri + pet_sd + pet_m + temp_bmf
            # all variables with releative influence>2
          mrd ~ pre_m+ sea_m+ tra_m+ sea_sd+ soil + tri+ tra_sd + pre_sd+  pet_m+ mat_m + area
          + mat_sd + temp_bmf + sub_trop_mbf + pet_sd + medit_fws + sub_trop_cf + deserts_x_shrub
          
          # secondary  regressions
          sub_trop_mbf ~ mat_m + pre_m + tra_m
          deserts_x_shrub ~ pet_m
          tra_sd ~ area
          sea_sd ~ area
          soil ~ area
"
sem_full <- sem(tsem_full, data = dat_no.na,  estimator="ML") 
sem_full_mlm <- sem(tsem_full, data = dat_no.na,  estimator="MLM") # use the robust output which deals well with non-normality and  corrects chi square rsults and standard errors
sem_full_se <- sem(tsem_full, data = dat_no.na,  se = "bootstrap", bootstrap = 500)
sem_full_test <- sem(tsem_full, data = dat_no.na, test="bootstrap", bootstrap = 500)
sem_full_both <- sem(tsem_full, data = dat_no.na,  se = "bootstrap", test="bootstrap", bootstrap = 500)

#summary(sem_full, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
#compare p value
mlm <- parameterestimates(sem_full_mlm, standardized = TRUE)[1:56,c("lhs", "op", "rhs", "pvalue")]
both <- parameterestimates(sem_full_both, boot.ci.type = "bca.simple", standardized = TRUE)[1:56,c("lhs", "op", "rhs", "pvalue")]
se <- parameterestimates(sem_full_se, boot.ci.type = "bca.simple", standardized = TRUE)[1:56,c("lhs", "op", "rhs", "pvalue")]
test <- parameterestimates(sem_full_test, boot.ci.type = "bca.simple", standardized = TRUE)[1:56,c("lhs", "op", "rhs", "pvalue")]
se$pvalue_se <- se$pvalue
test$pvalue_test <- test$pvalue
mlm$pvalue_mlm <- mlm$pvalue
se <- se[,-grep("pvalue$",names(se))]
test <- test[,-grep("pvalue$",names(test))]
mlm <- mlm[,-grep("pvalue$",names(mlm))]
res <- merge(both, se, all.x=TRUE)
res <- merge(res, test, all = TRUE)
res <- merge(res, mlm, all = TRUE)
res$variable <- paste0(res$lhs,res$op, res$rhs)
library(tidyverse)
library(ggplot2)
res2 <- pivot_longer(res, starts_with("pvalue"))
ggplot(res2, aes(x=name, y=value, col=variable))+
  geom_line(aes(group=variable, col=variable))+
  scale_y_log10()+
  geom_hline(yintercept=0.05, col="red")

# FULL MODEL ##########################################################################################################

## STEP 2
tsem_full <- "
          # regressions
          sr_trans ~ soil + sub_trop_mbf + sub_trop_dbf+ area+ pre_sd+ mrd+ mont_gs + pre_m+ sea_sd+ tra_m
          + mat_m +tra_sd+ mat_sd + sea_m + medit_fws + temp_gss + tri + pet_sd + pet_m + temp_bmf
          mrd ~ pre_m+ sea_m+ tra_m+ sea_sd+ soil + tri+ tra_sd + pre_sd + pet_m+ mat_m + area
          + mat_sd + temp_bmf + sub_trop_mbf + pet_sd + medit_fws + sub_trop_cf + deserts_x_shrub
          
          # secondary  regressions
          sub_trop_mbf ~ mat_m + pre_m + tra_m
          deserts_x_shrub ~ pet_m
          tra_sd ~ area
          sea_sd ~ area
          soil ~ area
"
tsem_full_mods <- "
      # regressions
          sr_trans ~ soil + sub_trop_mbf + sub_trop_dbf+ area+ pre_sd+ mrd+ mont_gs + pre_m+ sea_sd+ tra_m
          + mat_m +tra_sd+ mat_sd + sea_m + medit_fws + temp_gss + tri + pet_sd + pet_m + temp_bmf
          mrd ~ pre_m+ sea_m+ tra_m+ sea_sd+ soil + tri+ tra_sd + pre_sd + pet_m+ mat_m + area
          + mat_sd + temp_bmf + sub_trop_mbf + pet_sd + medit_fws + sub_trop_cf + deserts_x_shrub
      # secondary  regressions
          sub_trop_mbf ~ mat_m + pre_m + tra_m
          deserts_x_shrub ~ pet_m
          tra_sd ~ area
          sea_sd ~ area
          soil ~ area
      # mods
          sea_sd ~ sea_m 
          tra_sd ~ tra_m
          deserts_x_shrub ~ mat_m
          tra_sd ~ sea_sd
          sea_sd ~ pre_sd
          sea_sd ~ pet_sd
          tra_sd~~sea_sd
          deserts_x_shrub ~ pre_m
          deserts_x_shrub ~ mrd
          sea_sd ~ mat_sd
"
sem_full_mods <- sem(tsem_full_mods, data = dat_no.na,  estimator="MLM") 
fitmeasures(sem_full_mods, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(sem_full_mods, sort. = TRUE)[1:15,]
summary(sem_full_mods, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)

## STEP 2 #####
# check potential full regressions for multicollinearity
temp <- dat_no.na[,!grepl("sr$|elev_sd|number_biomes", names(dat_no.na))]
lm_sr <- lm(sr_trans ~ ., data=temp)
sort(car::vif(lm_sr))
# check out mat_m and pet_m
lm_sr_red <- lm(sr_trans ~ soil + sub_trop_mbf + area + mrd + sub_trop_dbf+ pre_sd+, 
            data=dat_no.na)
sort(car::vif(lm_sr_red))
lm_sr_red2 <- lm(sr_trans ~ soil + sub_trop_mbf + sub_trop_dbf+ area+ pre_sd+ mrd+ mont_gs + pre_m+ sea_sd+ tra_m
            + mat_m +tra_sd+ mat_sd + sea_m + medit_fws + temp_gss + tri + pet_sd + temp_bmf, 
            data=dat_no.na)
sort(car::vif(lm_sr_red2))

summary(lm_sr)
summary(lm_sr_red)
summary(lm_sr_red2)
# removing mat_m is more efficient in reducing VIFs, removing pet_m leaves more explanatory power.
# --> Remove pet_m


lm_mrd <- lm(mrd ~ pre_m+ sea_m+ tra_m+ sea_sd+ soil + tri+ tra_sd + pre_sd+ pet_m+ mat_m + area
             + mat_sd + temp_bmf + sub_trop_mbf + pet_sd + medit_fws + sub_trop_cf + deserts_x_shrub,
             data=dat_no.na)
sort(car::vif(lm_mrd))
# check out mat_m and pet_m
lm_mrd_red <- lm(mrd ~ pre_m+ sea_m+ tra_m+ sea_sd+ soil + tri+ tra_sd + pre_sd+ pet_m + area
                 + mat_sd + temp_bmf + sub_trop_mbf + pet_sd + medit_fws + sub_trop_cf + deserts_x_shrub,
                 data=dat_no.na)
sort(car::vif(lm_mrd_red))
lm_mrd_red2 <- lm(mrd ~ pre_m+ sea_m+ tra_m+ sea_sd+ soil + tri+ tra_sd + pre_sd+ mat_m + area
                  + mat_sd + temp_bmf + sub_trop_mbf + pet_sd + medit_fws + sub_trop_cf + deserts_x_shrub,
                  data=dat_no.na)
sort(car::vif(lm_mrd_red2))
summary(lm_mrd)$r.squared
summary(lm_mrd_red)$r.squared
summary(lm_mrd_red2)$r.squared
# removing mat_m is more efficient in reducing VIFs, and leaves slightly more explanatory power. 
# --> Remove mat_m

lm_subtropmbf <- lm( sub_trop_mbf ~ mat_m + pre_m + tra_m, data=dat_no.na)
sort(car::vif(lm_subtropmbf))
# no action required

## STEP 4
step.sr <- step(lm_sr_red2)
#variables to remove: pet_sd, tra_m, temp_gss, tri, mrd, sea_sd
# we keep mrd though as its part of our theoretical model structure

step.mrd <- step(lm_mrd_red)
#variables to remove: pre_sd, pet_sd, sea_m, pet_m, tri, mat_sd

## STEP 5
# mod summary: removed pet_m from SR, mat_m from MRD based on vif() results
# removed all suggested variables from  step()
tsem_mod <- "
          # regressions
          sr_trans ~ soil + sub_trop_mbf + sub_trop_dbf+ area+ pre_sd+ mrd+ mont_gs + pre_m
          + mat_m +tra_sd+ mat_sd + medit_fws + pet_m + temp_bmf
          mrd ~ pre_m+  tra_m+ sea_sd+ soil +  tra_sd + mat_m + area
           + temp_bmf + sub_trop_mbf + medit_fws + sub_trop_cf + deserts_x_shrub
          
          # secondary regressions
          sub_trop_mbf ~ mat_m + pre_m + tra_m
          deserts_x_shrub ~ pet_m
          tra_sd ~ area
          sea_sd ~ area
          soil ~ area
          
          # mod indicies mods
          deserts_x_shrub ~ mat_m + pre_m
          tra_sd ~ tra_m
          tra_sd ~ sea_sd
          sea_sd ~ mat_m
          sea_sd ~ pet_m
"
sem_mod <- sem(tsem_mod, data = dat_no.na,  estimator="MLM") 
fitmeasures(sem_mod, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(sem_mod, sort. = TRUE)[1:15,]
summary(sem_mod, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)


tsem_full_no_sd <- "
      # regressions
          sr_trans ~ soil + sub_trop_mbf + sub_trop_dbf+ area+ mrd+ mont_gs + pre_m+ tra_m
          + mat_m + sea_m + medit_fws + temp_gss + tri + pet_m + temp_bmf
          mrd ~ pre_m+ sea_m+ tra_m+ soil + tri+ pet_m+ mat_m + area
           + temp_bmf + sub_trop_mbf + medit_fws + sub_trop_cf
      # secondary  regressions
          sub_trop_mbf ~ mat_m + pre_m + tra_m
          deserts_x_shrub ~ pet_m
          soil ~ area
      # mods
          deserts_x_shrub  ~ mat_m+ pre_m
          sub_trop_mbf ~ temp_bmf
          #sub_trop_mbf ~ area
          sub_trop_mbf ~ pet_m
          #sub_trop_mbf ~ tri
          sub_trop_mbf  ~    sub_trop_dbf
          
"
sem_full_no_sd <- sem(tsem_full_no_sd, data = dat_no.na,  estimator="MLM") 
fitmeasures(sem_full_no_sd, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
summary(sem_full_no_sd, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
modificationindices(sem_full_no_sd, sort. = TRUE)[1:15,]


tsem_full_no_spatial_temp_hg <- " # (sea_sd + tra_sd)
          # regressions
          sr_trans ~ soil + sub_trop_mbf + sub_trop_dbf+ area+ pre_sd+ mrd+ mont_gs + pre_m+ tra_m
          + mat_m + mat_sd + sea_m + medit_fws + temp_gss + tri + pet_sd + pet_m + temp_bmf
          mrd ~ pre_m+ sea_m+ tra_m+  soil + tri+  pre_sd + pet_m+ mat_m + area
          + mat_sd + temp_bmf + sub_trop_mbf + pet_sd + medit_fws + sub_trop_cf + deserts_x_shrub
          
          # secondary  regressions
          sub_trop_mbf ~ mat_m + pre_m + tra_m
          deserts_x_shrub ~ pet_m
          soil ~ area
          
          # mods 
          #deserts_x_shrub  ~ mat_m + pre_m
"
sem_full_no_spatial_temp_hg <- sem(tsem_full_no_spatial_temp_hg, data = dat_no.na,  estimator="MLM") 
fitmeasures(sem_full_no_spatial_temp_hg, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(sem_full_no_spatial_temp_hg, sort. = TRUE)[1:15,]


tsem_mod_no_spatial_temp_hg <- " # (sea_sd + tra_sd)
          # regressions
          sr_trans ~ soil + sub_trop_mbf + sub_trop_dbf+ area+ pre_sd+ mrd+ mont_gs + pre_m+ tra_m
          + mat_m + mat_sd + sea_m + medit_fws + temp_gss + tri + pet_sd + pet_m + temp_bmf
          mrd ~ pre_m+ sea_m+ tra_m+  soil + tri+  pre_sd + pet_m+ mat_m + area
          + mat_sd + temp_bmf + sub_trop_mbf + pet_sd + medit_fws + sub_trop_cf + deserts_x_shrub
          
          # secondary  regressions
          sub_trop_mbf ~ mat_m + pre_m + tra_m
          deserts_x_shrub ~ pet_m
          #soil ~ area
          
          # mods 
          deserts_x_shrub  ~ mat_m + pre_m
"
sem_mod_no_spatial_temp_hg <- sem(tsem_mod_no_spatial_temp_hg, data = dat_no.na,  estimator="MLM") 
fitmeasures(sem_mod_no_spatial_temp_hg, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(sem_mod_no_spatial_temp_hg, sort. = TRUE)[1:10,]










tsem_full_hg <- "
      # latent variabes
          hg_spatial =~ pre_sd + pet_sd + mat_sd + soil
      # regressions
          sr_trans ~ hg_spatial + sub_trop_mbf + sub_trop_dbf+ area+ mrd+ mont_gs + pre_m+ sea_sd+ tra_m
          + mat_m +tra_sd+ sea_m + medit_fws + temp_gss + tri + pet_m + temp_bmf
          mrd ~ pre_m+ sea_m+ tra_m+ sea_sd + tri+ tra_sd + pet_m+ mat_m + area
           + temp_bmf + sub_trop_mbf + medit_fws + sub_trop_cf + deserts_x_shrub + hg_spatial
          
      # secondary  regressions
          sub_trop_mbf ~ mat_m + pre_m + tra_m
          deserts_x_shrub ~ pet_m
          hg_spatial ~ area
          #sea_sd ~ area
          #soil ~ area
"
sem_full_hg <- sem(tsem_full_hg, data = dat_no.na,  estimator="MLM") 
fitmeasures(sem_full_hg, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
cor(dat_no.na[,c("pre_sd", "pet_sd", "mat_sd", "tri", "soil", "sr_trans")])
summary(sem_full_hg, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)



tsem_full_mods <- " # removing sea_sd + tra_sd
      # regressions
          sr_trans ~ soil + sub_trop_mbf + sub_trop_dbf+ area+ pre_sd+ mrd+ mont_gs + pre_m+ tra_m
          + mat_m +mat_sd + sea_m + medit_fws + temp_gss + tri + pet_sd + pet_m + temp_bmf
          mrd ~ pre_m+ sea_m+ tra_m+ soil + tri+  pre_sd + pet_m+ mat_m + area
          + mat_sd + temp_bmf + sub_trop_mbf + pet_sd + medit_fws + sub_trop_cf + deserts_x_shrub
      # secondary  regressions
          sub_trop_mbf ~ mat_m + pre_m + tra_m + area
          deserts_x_shrub ~ pet_m
          soil ~ area
      # mods
          deserts_x_shrub ~ mat_m
          deserts_x_shrub ~ pre_m
          sub_trop_mbf ~ pet_m
          sub_trop_mbf ~ tri
"
sem_full_mods <- sem(tsem_full_mods, data = dat_no.na,  estimator="MLM") 
fitmeasures(sem_full_mods, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(sem_full_mods, sort. = TRUE)[1:15,]



# mod_base5_ind <- "
#           # regressions
#           sr_trans ~ soil + sub_trop_mbf + sub_trop_dbf+ area+ pre_sd + mrd 
#           mrd ~ pre_m+ sea_m+ tra_m + sea_sd + soil + tri
#           
#           # area influence on spatial heterogeneity
#           soil ~ area
# #          pre_sd ~ area
# #          tra_sd ~ area
#           sea_sd ~ area
# #          mat_sd ~ area
# 
#           # climate influence on habitat
#           sub_trop_mbf ~ mat_m + pre_m + tra_m 
# "
# fsem_base5_ind <- sem(mod_base5_ind, data = dat_no.na, estimator="MLM")
# fitmeasures(fsem_base5_ind, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
# modificationindices(fsem_base5_ind, sort = TRUE)[1:10,]

# WORK HERE ON FORWARD SELECTION
# theory: base 5 + mrd + modificatins to get acceptable model fit
# remove non-significant variables that are not connected with others (subtrop_dbf from SR, tra_m from MRD)
# + one pet_m for temperature varibale since MAT is causing lots of VI, 
# add more variables as long as they do not make model fit worse, preferrably ones that are present already 
mod_base5_ind <- "
          # regressions
          sr_trans ~ soil + sub_trop_mbf + area+ pre_sd + mrd 
          + pet_m
          + sea_sd + sea_m + tri
           
          mrd ~ pre_m+ sea_m + sea_sd + soil
          + pet_m
          #+ pre_sd
          
          # area influence on spatial heterogeneity
          soil ~ area
          sea_sd ~ area + sea_m

          # climate influence on habitat
          sub_trop_mbf ~ pre_m + tra_m +  pet_m + area
          
          # mods
          # sea_sd  ~ sea_m
          # sea_sd ~ pre_sd
          # sub_trop_mbf  ~ area
          # mrd  ~ sub_trop_mbf
          # sea_sd  ~        pre_m
          # mrd ~ area
          # soil ~ tri
"

fsem_base5_ind <- sem(mod_base5_ind, data = dat_no.na, estimator="MLM")
fitmeasures(fsem_base5_ind, c("chisq", "df", "pvalue", "cfi.robust", "rmsea.robust", "aic"))
modificationindices(fsem_base5_ind, sort = TRUE)[1:15,]
summary(fsem_base5_ind, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
semPaths(object = fsem_base5_ind, layout = "circle2", rotation = 1,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
         exoCov = FALSE, exoVar = FALSE,
         nCharNodes = 0, fade=FALSE,
         edgeLabels = sem_sig_labels(fsem_base5_ind))



