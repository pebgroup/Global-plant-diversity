# Other model tests
# produces data_for_SEM.rds


org <- st_drop_geometry(shp)

# remove not needed columns 
## mangroves and flooded grasslands get excluded because they are hierarchical very fine constructs compared to other biome types and therefore pretty rare

rownames(org) <- org$LEVEL_3_CO
dat <- org[,-grep("_n|lng|lat|mangroves|flooded|MRD.sd|LEVEL_3_CO|ID|LEVEL_NAM|REGION_NAM|CONTINENT",
                  names(org))]


dat_no.na <- na.omit(dat)
dat_no.na$level3 <- row.names(dat_no.na)
#saveRDS(dat_no.na, "processed_data/data_for_SEM.rds")



#### CORRELATION #############################################################################

cor.dat_no.na<- cor(dat_no.na[,-grep("level3",names(dat_no.na))])
p.dat_no.na <- cor_pmat(dat_no.na[,-grep("level3",names(dat_no.na))])

my_corrplot(cor.dat_no.na, lab=TRUE, p.mat = p.dat_no.na, insig = "blank",
            tl.cex = 8, digits = 1, type = "lower", hc.order = FALSE,
            lab_size = 3, highlight = TRUE)+
  theme(axis.text.x = element_text(margin=margin(0,0,0,0)),  # Order: top, right, bottom, left
        axis.text.y = element_text(margin=margin(0,0,0,0)),
        plot.margin = margin(0, 0, 0, 0, "cm"))
ggsave("publish_figures/correlation_april2021.png", width=7, height=6, units = "in", dpi = 600)


# monitor pet_m + mat_m, mat_m + tra_m for multicollinearity



#### GBM CONTROLS ############################################################################

# parameter tuning using caret
gbmControl <- trainControl(method = "repeatedcv", number = 10,
                           repeats = 3, savePredictions = "final",
                           returnResamp ="final")

gbmGrid <-  expand.grid(interaction.depth = c(1, 2, 3), 
                        n.trees = (1:10)*50, 
                        shrinkage = c(0.1, 0.01, 0.001),
                        n.minobsinnode = c(5, 10, 15))


# GBMs  #####################################################################################

set.seed(123)
gbm_SR <- train(sr ~ ., data = dat_no.na[,-grep("level3", names(dat_no.na))], method = "gbm", 
                     trControl = gbmControl, verbose=FALSE,
                     tuneGrid=gbmGrid)
gbm_SR$results[row.names(gbm_SR$bestTune),]
gbm_sr <- ggplot(summary(gbm_SR), aes(x=reorder(var, rel.inf), y=rel.inf))+
  geom_bar(stat="identity")+
  coord_flip()+
  xlab("")+ylab("Relative influence")

set.seed(100) 
gbm_MRD <- train(mrd ~ ., data = dat_no.na[,-grep("sr|level3",names(dat_no.na))], method = "gbm", 
                     trControl = gbmControl, verbose=FALSE,
                     tuneGrid=gbmGrid)
gbm_MRD$results[row.names(gbm_MRD$bestTune),]
gbm_mrd <- ggplot(summary(gbm_MRD), aes(x=reorder(var, rel.inf), y=rel.inf))+
  geom_bar(stat="identity")+
  coord_flip()+
  xlab("")+ylab("Relative influence")

plot_grid(gbm_sr, gbm_mrd, labels = c("a)", "b)"), ncol = 2)
ggsave("publish_figures/varImp_gbm_April2021.png", width=7, height=4, units = "in", dpi = 600)







# rm(list=ls())
# library(sf)
# library(lavaan)
# library(semPlot)
# library(modEvA)
# library(fitdistrplus)
# source("scripts/clean_semplot_functions.R") # script for modified plotting functions

# load data ###################################################################################################
#dat_no.na <-  readRDS("processed_data/data_for_SEM.rds")
names(dat_no.na) <- sub("_mean$", "_m", names(dat_no.na)) # shorten mean to _m

# sqrt transform species richness for better normal distribution
dat_no.na$sr_trans <- sqrt(dat_no.na$sr)
#hist(dat_no.na$sr_trans)


# test distribution of potential response variables ############################################################
## distributions SR and MRD 
hist(dat_no.na[,c("sr", "mrd")])
hist(dat_no.na[,c("sr_trans", "mrd")])
hist(dat_no.na$sr, breaks=100)
# mrd seems normally distributed, SR not. Lets test:
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

## prs_sd
descdist(dat_no.na$prs_sd, discrete = TRUE, boot=500)
hist(log10(dat_no.na$prs_sd), breaks=50)
shapiro.test(log(dat_no.na$prs_sd))

## soil
descdist(dat_no.na$soil, discrete = TRUE, boot=500)
hist(dat_no.na$soil)
shapiro.test(dat_no.na$soil)

# transformations:
dat_no.na$tra_sd <- sqrt(dat_no.na$tra_sd)  
dat_no.na$prs_sd <- sqrt(dat_no.na$prs_sd)  

#  scale data #############################################################################################
#z-transformation to avoid variance issues
dat_no.na <- dat_no.na[,-grep("level3",names(dat_no.na))]
dat_no.na <- apply(dat_no.na, 2, ztrans)
dat_no.na <- as.data.frame(dat_no.na)

saveRDS(dat_no.na, "processed_data/sem_input_data.rds")
# feed this into model selection




## multicollinearity ######################################################################################
# check potential full regressions for multicollinearity
temp <- dat_no.na[,!grepl("sr$|deserts_x_shrub|medit_fws|temp_cf|
                          sub_trop_cf|sub_trop_gss|temp_gss|boreal|tundra", names(dat_no.na))]
lm_sr <- lm(sr_trans ~ . , data=temp)
sort(car::vif(lm_sr))

# check out mat_m and pet_m
lm_sr1 <- lm(sr_trans ~ . -mat_m, data=temp)
sort(car::vif(lm_sr1))
lm_sr2 <- lm(sr_trans ~ . -pet_m, data=temp)
sort(car::vif(lm_sr2))

summary(lm_sr)
summary(lm_sr1)
summary(lm_sr2)
# removing mat_m is more efficient in reducing VIFs, removing pet_m leaves more explanatory power.
# --> If problems occurr, remove pet_m and keep mat_m




lm_mrd <- lm(mrd ~ . -sr_trans, data=temp)
sort(car::vif(lm_mrd))

# check out mat_m and pet_m
lm_mrd1 <- lm(mrd ~ . -sr_trans -mat_m, data=temp)
sort(car::vif(lm_mrd1))
lm_mrd2 <- lm(mrd ~ . -sr_trans -pet_m, data=temp)
sort(car::vif(lm_mrd2))
summary(lm_mrd)$r.squared
summary(lm_mrd1)$r.squared
summary(lm_mrd2)$r.squared
# removing mat_m is more efficient in reducing VIFs, removing pet_m leaves slightly more explanatory power.
# Same decision as above if necessary


# indirect effects:
lm_subtropmbf <- lm(sub_trop_mbf ~ mat_m + pre_m + tra_m, data=dat_no.na)
summary(lm_subtropmbf)
sort(car::vif(lm_subtropmbf))
# no action required

