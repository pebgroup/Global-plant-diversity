# Data prep for SEM
# produces data_for_SEM.rds

library(sf)
library(ggcorrplot)
library(caret)
library(gbm)
library(cowplot)
library(ggplot2)
library(data.table)
library(parallel)
theme_set(theme_bw())
source("0_functions.R")

rm(list = setdiff(ls(), lsf.str()))  
shp <- readRDS("../processed_data/shp_object_fin_analysis.RDS")


org <- st_drop_geometry(shp)

# remove not needed columns 
## We exclude mangroves and flooded grasslands because they are hierarchical
## very fine constructs compared to other biome types and therefore pretty rare

rownames(org) <- org$LEVEL_3_CO
dat <- org[,-grep("_n|mangroves|flooded|LEVEL_3_CO|ID|LEVEL_NAME|REGION_NAM|CONTINENT",
                  names(org))]


dat_no.na <- na.omit(dat)
dat_no.na$level3 <- row.names(dat_no.na)
dim(dat_no.na)



#### CORRELATION ##############################################################

cor.dat_no.na<- cor(dat_no.na[,-grep("level3",names(dat_no.na))])
p.dat_no.na <- cor_pmat(dat_no.na[,-grep("level3",names(dat_no.na))])

my.corrplot(cor.dat_no.na, lab=F, p.mat = p.dat_no.na, insig = "blank",
            tl.cex = 8, digits = 1, type = "lower", hc.order = FALSE,
            lab_size = 4, highlight = TRUE, method="square", lab_col = "grey20")+
  scale_y_discrete(position = "right")+
  theme(axis.text.x = element_text(margin=margin(0,0,0,0)),  # Order: top, right, bottom, left
        axis.text.y = element_text(margin=margin(0,0,0,0)),
        plot.margin = margin(0, 0, 0, 0, "cm"),
        panel.grid = element_blank(),
        legend.position = c(.15, .8), legend.text.align = 1)
ggsave("../figures/correlation_2022.png", width=7, height=6, units = "in", dpi = 600, bg = "white")




# Multicollinearity ########################################################
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
lm_sr3 <- lm(sr_trans ~ . -mat_m -pre_lgm_ano_m, data=temp)
sort(car::vif(lm_sr3))

summary(lm_sr)
summary(lm_sr1)
summary(lm_sr2)
summary(lm_sr3)
# removing mat_m is more efficient in reducing VIFs, removing pet_m leaves more explanatory power.
# If necessary, remove pet_m before mat_m
# After removing either pet or mat, precipitation LGM anomaly stays a problem. Removing solves this.




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
# removing mat_m is more efficient in reducing VIFs, removing pet_m leaves slightly less explanatory power.
# No preference. See which one performs better, avoid having both in same regression
# pre_lgm_ano_m sis a problem here too


lm_mdr <- lm(mdr ~ . -sr_trans, data=temp)
sort(car::vif(lm_mdr))
# same dynamic as for MRD


# indirect effects:
lm_subtropmbf <- lm(sub_trop_mbf ~ mat_m + pre_m + tra_m, data=dat_no.na)
summary(lm_subtropmbf)
sort(car::vif(lm_subtropmbf))
# no action required

# remove pet and pre_lgm_ano_m
dat_no.na <- dat_no.na[,-grep("pre_lgm|pet",
                  names(dat_no.na))]


dat_no.na <- na.omit(dat_no.na)
dat_no.na$level3 <- row.names(dat_no.na)
dim(dat_no.na)


#### GBM CONTROLS #############################################################

# parameter tuning using caret
gbmControl <- trainControl(method = "repeatedcv", number = 10,
                           repeats = 3, savePredictions = "final",
                           returnResamp ="final")

gbmGrid <-  expand.grid(interaction.depth = c(1, 2, 3), 
                        n.trees = (1:10)*50, 
                        shrinkage = c(0.1, 0.01, 0.001),
                        n.minobsinnode = c(5, 10, 15))



# SR ----------------------------------------------------------------------

# define function
run.gbm <- function(i, seeds, data, trControl=gbmControl,
                    tuneGrid=gbmGrid){
  set.seed(s[i]) 
  if(!i%%1)cat(i,"\r")
  gbm_temp <- train(sr ~ ., data, method = "gbm", 
                  trControl = gbmControl, verbose=FALSE,
                  tuneGrid = gbmGrid)
  temp <- list(summary(gbm_temp)$var[1:16], s[i])
  return(temp)
}

n_cores=4
s <- seq(540,645,1)
system.time(
  sr_list <- mclapply(1:length(s), seeds = s, run.gbm,
                      data = dat_no.na[,-grep("level3", names(dat_no.na))],
                      #response=sr,
                      mc.cores=n_cores)
)
#saveRDS(sr_list, "sr_list100_lgm.rds")
#saveRDS(sr_list, "sr_list100.rds")

# check
#sr_list <- readRDS("sr_list100_lgm.rds")
filnames <- dir("../processed_data/gbm", full.names = T)
filnames <- filnames[grep("sr_list100_lgm5", filnames)]
length(filnames)
sr_list <- sapply(filnames, readRDS, simplify = F) 
# combine files
sr_list <- unlist(sr_list, recursive = F)


eq <- matrix(sapply(sr_list, "[[", 1), ncol=31, byrow = T)
eq <- melt(eq)
ggplot(eq, aes(value))+
  geom_bar()+
  facet_wrap(~Var2, scales = "free")+
  theme(axis.text.x = element_text(size=6, angle = 45, hjust=1),
        strip.background = element_blank())
#ggplot(eq, aes(x=Var2, group=value, fill=value))+
#  geom_bar(position = "fill")
ggsave("../figures/varImp_SR_gbm100_runs.png", width=7, height=7, units = "in", dpi = 600)

eqr <-tapply(eq$value, eq$Var2, function(x){names(sort(table(x), decreasing = T))})
fin_sr <- sort(tapply(eq$Var2, eq$value, function(x){sum(x)/length(x)}), decreasing=F)

# select variables that score lowest position relative to their number of
# appearances in different positions
fin_sr <- data.frame(fin_sr, variable=names(fin_sr))
fin_sr <- fin_sr[order(fin_sr$fin_sr, decreasing = T),]
fin_sr$variable <- factor(fin_sr$variable, levels = fin_sr$variable)
tot.position.sr <- ggplot(fin_sr[c((nrow(fin_sr)-30):nrow(fin_sr)),], aes(y=variable, x=fin_sr))+
  #xlab("Average variable position")+
  geom_bar(stat = "identity")+
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size=9))



# MRD ---------------------------------------------------------------------

# redefine function for MRD
run.gbm <- function(i, seeds, data, trControl=gbmControl,
                    tuneGrid=gbmGrid){
  set.seed(s[i]) 
  if(!i%%1)cat(i,"\r")
  gbm_temp <- train(mrd ~ ., data, method = "gbm", 
                    trControl = gbmControl, verbose=FALSE,
                    tuneGrid = gbmGrid)
  temp <- list(summary(gbm_temp)$var[1:16], s[i])
  return(temp)
}

n_cores=4
s <- seq(670,770,1)
system.time(
  mrd_list <- mclapply(1:length(s), seeds = s, run.gbm,
                       data = dat_no.na[,-grep("sr|level3|mdr",names(dat_no.na))],
                       mc.cores=n_cores)
)
saveRDS(mrd_list, "mrd_list100_lgm.rds")
#saveRDS(mrd_list, "mrd_list100.rds")

# check
#mrd_list <- readRDS("mrd_list100_lgm.rds")
filnames <- dir("../processed_data/gbm", full.names = T)
filnames <- filnames[grep("mrd_list100_lgm5_", filnames)]
length(filnames)
mrd_list <- sapply(filnames, readRDS, simplify = F) 
# combine files
mrd_list <- unlist(mrd_list, recursive = F)

eq <- matrix(sapply(mrd_list, "[[", 1), ncol=29, byrow = T)
eq <- melt(eq)
ggplot(eq, aes(value))+
  geom_bar()+
  facet_wrap(~Var2, scales = "free")+
  theme(axis.text.x = element_text(size=7,angle = 45, hjust=1),
        strip.background = element_blank())
#ggplot(eq, aes(x=Var2, group=value, fill=value))+
#  geom_bar(position = "fill")
ggsave("../figures/varImp_MRD_gbm100_runs.png", width=7, height=7, units = "in", dpi = 600)

eqr <-tapply(eq$value, eq$Var2, function(x){names(sort(table(x), decreasing = T))})
fin_mrd <- sort(tapply(eq$Var2, eq$value, function(x){sum(x)/length(x)}), decreasing=F)

# select variables that score lowest position relative to their number of
# appearances in different positions
names(fin_mrd[1:14])
fin_mrd <- data.frame(fin_mrd, variable=names(fin_mrd))
fin_mrd <- fin_mrd[order(fin_mrd$fin_mrd, decreasing = T),]
fin_mrd$variable <- factor(fin_mrd$variable, levels = fin_mrd$variable)
tot.position.mrd <- ggplot(fin_mrd[c((nrow(fin_mrd)-28):nrow(fin_mrd)),], aes(y=variable, x=fin_mrd))+
  #xlab("Average variable position")+
  geom_bar(stat = "identity")+
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size=9),
        axis.title.y = element_blank())

plot_grid(tot.position.sr, tot.position.mrd, ncol=2, labels = c("SR", "MRD"))

# DivRate -----------------------------------------------------------------

# redefine function for MDR
run.gbm <- function(i, seeds, data, trControl=gbmControl,
                    tuneGrid=gbmGrid){
  set.seed(s[i]) 
  if(!i%%1)cat(i,"\r")
  gbm_temp <- train(mdr ~ ., data, method = "gbm", 
                    trControl = gbmControl, verbose=FALSE,
                    tuneGrid = gbmGrid)
  temp <- list(summary(gbm_temp)$var[1:16], s[i])
  return(temp)
}

n_cores=4
s <- seq(770,870,1)
system.time(
  mdr_list <- mclapply(1:length(s), seeds = s, run.gbm,
                       data = dat_no.na[,-grep("sr|level3|mrd",names(dat_no.na))],
                       mc.cores=n_cores)
)
#saveRDS(mdr_list, "mdr_list100.rds")

# check
#mdr_list <- readRDS("mdr_list100.rds")
filnames <- dir("../processed_data/gbm", full.names = T)
filnames <- filnames[grep("mdr_list100_lgm5_", filnames)]
length(filnames)
mdr_list <- sapply(filnames, readRDS, simplify = F) 
# combine files
mdr_list <- unlist(mdr_list, recursive = F)
eq <- matrix(sapply(mdr_list, "[[", 1), ncol=29, byrow = T)
eq <- melt(eq)
ggplot(eq, aes(value))+
  geom_bar()+
  facet_wrap(~Var2, scales = "free")+
  theme(axis.text.x = element_text(size=7,angle = 45, hjust=1),
        strip.background = element_blank())
#ggplot(eq, aes(x=Var2, group=value, fill=value))+
#  geom_bar(position = "fill")
ggsave("../figures/varImp_MDR_gbm100_runs.png", width=8, height=7, units = "in", dpi = 600)

eqr <-tapply(eq$value, eq$Var2, function(x){names(sort(table(x), decreasing = T))})
fin_mdr <- sort(tapply(eq$Var2, eq$value, function(x){sum(x)/length(x)}), decreasing=F)

# select variables that score lowest position relative to their number of
# appearances in different positions
names(fin_mdr[1:13])
fin_mdr <- data.frame(fin_mdr, variable=names(fin_mdr))
fin_mdr <- fin_mdr[order(fin_mdr$fin_mdr, decreasing = T),]
fin_mdr$variable <- factor(fin_mdr$variable, levels = fin_mdr$variable)
(tot.position.mdr <- ggplot(fin_mdr[c((nrow(fin_mdr)-28):nrow(fin_mdr)),], aes(y=variable, x=fin_mdr))+
  #xlab("Average variable position")+
  geom_bar(stat = "identity")+
  theme(axis.title.x = element_blank(), axis.text.y = element_text(size=9),
  axis.title.y = element_blank()))


temp <- plot_grid(ncol = 3,
          labels = c("A", "B", "C"), label_fontface = "plain",
          tot.position.sr, tot.position.mrd, tot.position.mdr)
ggdraw(add_sub(temp, "Average variable position", vpadding=grid::unit(0.5,"lines"),y=3, x=0.5, vjust=4.5, 
               size = 11))
ggsave("../figures/total_var_position_GBM.png", width=7, height=5, units = "in", dpi = 600, bg = "white")



selected_variables <- list(fin_sr, fin_mrd, fin_mdr)
saveRDS(selected_variables, "../processed_data/selected_variables.rds")
saveRDS(dat_no.na, "../processed_data/data_for_SEM.rds")






# load data ##################################################################
rm(list=ls())
library(modEvA)
library(fitdistrplus)

dat_no.na <-  readRDS("../processed_data/data_for_SEM.rds")
names(dat_no.na) <- sub("_mean$", "_m", names(dat_no.na)) # shorten mean to _m


# test distribution of potential response variables ############################################################
## distributions SR and MRD 
# sqrt transform species richness for better normal distribution
dat_no.na$sr_trans <- sqrt(dat_no.na$sr)
hist(dat_no.na$sr_trans)

hist(dat_no.na$sr)
hist(dat_no.na$sr_trans)
hist(dat_no.na$mrd)
hist(dat_no.na$mdr)

# mrd seems normally distributed. Lets test:
descdist(dat_no.na$mrd, discrete = FALSE, boot=500)
shapiro.test(dat_no.na$mrd) # almost normal distribution
mrd_norm_fit <- fitdist(dat_no.na$mrd, "norm",method="mle")
summary(mrd_norm_fit)
plot(mrd_norm_fit)

descdist(dat_no.na$sr_trans, discrete = TRUE, boot=500)
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

hist(dat_no.na$mdr, breaks=20)
hist(log1p(dat_no.na$mdr), breaks=20)

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
dat_no.na$mdr <- log1p(dat_no.na$mdr)

# Scale data ################################################################
# z-transformation to avoid variance issues in SEM
dat_no.na <- dat_no.na[,-grep("level3",names(dat_no.na))]
dat_no.na <- apply(dat_no.na, 2, ztrans)
dat_no.na <- as.data.frame(dat_no.na)

saveRDS(dat_no.na, "../processed_data/sem_input_data.rds")
# feed this into model selection


