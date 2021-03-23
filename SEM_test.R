# sem test ### https://nmmichalak.github.io/nicholas_michalak/blog_entries/2018/nrg01/nrg01.html

rm(list=ls())
#library(tidyverse)
#library(knitr)
library(lavaan)
#library(psych)
library(semPlot)
#library(MBESS)
#library(ggcorrplot)

#### tests ####
thirst_dat <- "mackinnon_2008_t3.1.csv" %>% read_csv()

thirst_dat %>% 
  headTail() %>% 
  kable()

thirst_dat %>% 
  select(room_temp, thirst, consume) %>% 
  pairs.panels()

# Idea: room_temp is associated with consume indirectly through thirst

mod1 <- "# a path
         thirst ~ a * room_temp

         # b path
         consume ~ b * thirst

         # c prime path 
         consume ~ cp * room_temp

         # indirect and total effects
         ab := a * b
         total := cp + ab"
set.seed(1234)
fsem1 <- sem(mod1, data = thirst_dat, se = "bootstrap", bootstrap = 100)
summary(fsem1, standardized = TRUE, fit.measures=TRUE)

parameterestimates(fsem1, boot.ci.type = "bca.simple", standardized = TRUE)
semPaths(object = fsem1,
         layout = "tree",
         rotation = 2,
         whatLabels="estimate",
         label.cex = 1.5,
         edge.label.cex = 1.5,
         what="std",
         edge.color="blue")
# Interpretation:Every 1&degF increase in room temperature was associated with an a = 0.339 (estimates, not std) (S.E. = 0.101) increase in thirstiness units. Adjusting for room temperature, every 1-unit increase in thirstiness was associated with drinking b = 0.451 (S.E. = 0.149) more deciliters of water. Increases in room temperature were associated with increases in water drinking indirectly through increases in thirstiness. Specifically, for every a = 0.339 unit increase in the association between room temperature and thirstiness, there was an ab = 0.153 (estimate from combined paths) (S.E. = 0.064) increase in deciliters of water people drank. Importantly, a bias-corrected bootstrapped confidence interval with 10,000 samples was above zero, 95% CI [0.06, 0.32]. Last, there was no sufficient evidence that room temperature was associated with how many deciliters of water people drank independent of its association with thirstiness, c’ = 0.208 (S.E. = 0.130).
#--> What this model shows is that increases in consume are not caused directly by room temperature, but indirectly via thirstiness.

#mediate(consume ~ room_temp + thirst, data = thirst_dat, n.iter = 10000) %>% print(short = FALSE)


#### my data ####
shp <- readRDS("data/shp_object_fin_analysis.RDS")
library(sf)
dat <- st_drop_geometry(shp)
# remove not needed columns 
dat <- dat[,-grep("trops|_n|lng|lat|z_score", names(dat))]
rownames(dat) <- dat$LEVEL_3_CO
dat <- dat[,c(grep("sr", names(dat)):ncol(dat))]

# two datasets: complete, resulting in less data (310) because of more NAs in sd. means only, resulting in 334 data points

dat_means <- dat[,-grep("lng|lat_bin|z_score|trops|binary|_n|_sd|_even|_simp", names(dat))]
dat_means <- na.omit(dat_means)
dat_no.na <- na.omit(dat)

# standardize all variables with z-transformation to avoid variance issues
ztrans <- function(x){(x - mean(x)) / sd(x)}

dat_means_scaled <- apply(dat_means[-grep("predominant_biome", names(dat_means))], 2, ztrans)
dat_means_scaled <- as.data.frame(dat_means_scaled)
dat_means_scaled$predominant_biome <- dat_means$predominant_biome

dat_scaled <- apply(dat_no.na[-grep("predominant_biome", names(dat_no.na))], 2, ztrans)
dat_scaled <- as.data.frame(dat_scaled)
dat_scaled$predominant_biome <- dat_no.na$predominant_biome


# Idea: mat_mean is associated with species richness (sr) indirectly through mrd. Temeprature increases speciation, therefore increases species richness
mod1 <- "# a path
         mrd ~ a * mat_mean
         # b path
         sr ~ b * mrd
         # c prime path 
         sr ~ cp * mat_mean 
         # indirect and total effects
         ab := a * b
         total := cp + ab"
set.seed(1234)
fsem1 <- sem(mod1, data = dat_means_scaled)#, se = "bootstrap", bootstrap = 1000)
summary(fsem1, standardized = TRUE, fit.measures=TRUE)
parameterestimates(fsem1, boot.ci.type = "bca.simple", standardized = TRUE)
semPaths(object = fsem1, layout = "tree", rotation = 2,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1.5, what="std")
# Interpret: no indirect effect of mat on sr through mrd. Direct effect is significant though - used to be not in the first run.

# test for variables that are strongly correlated:
# Idea: area is associated with species richness (sr) indirectly through soil
mod2 <- "# a path
         soil ~ a * area
         # b path
         sr ~ b * soil
         # c prime path 
         sr ~ cp * area
         # indirect and total effects
         ab := a * b
         total := cp + ab"
set.seed(1234)
fsem2 <- sem(mod2, data = dat_means_scaled) #, se = "bootstrap", bootstrap = 1000
summary(fsem2, standardized = TRUE, fit.measures=TRUE)
parameterestimates(fsem2, boot.ci.type = "bca.simple", standardized = TRUE)
semPaths(object = fsem2, layout = "tree", rotation = 2,
         whatLabels="std", label.cex = 1.5, edge.label.cex = 1.5, what="std")

# interesting, area has no direct effect on SR when you account for soil!


# Adding another variable to the regression
mod2.2 <- '
         soil ~ a * area
         sr ~ b * soil + cp * area + d* mat_mean
         
         # indirect and total effects
         ab := a * b
         total := cp + ab'
set.seed(1234)
fsem2.2 <- sem(mod2.2, data = dat_means_scaled, se = "bootstrap", bootstrap = 100) 
summary(fsem2.2, standardized = TRUE, fit.measures=TRUE)
parameterestimates(fsem2.2, boot.ci.type = "bca.simple", standardized = TRUE)
fitmeasures(fsem2.2, c("cfi", "tli", "rmsea", "srmr"))
semPaths(object = fsem2.2, layout = "tree", rotation = 2,
         whatLabels="std", label.cex = 1, edge.label.cex = 1, what="std")
# all values are more or less stable compared to the models before, plus this model has okish fit


### CFA _ Heterogeneity as a latent variable ####
hetero.model1 <- "
        heterogeneity =~ soil + tri + tra_mean  + sea_mean  + number_biomes + area
"
hetero1.fit <- cfa(hetero.model1, data=dat_means_scaled)
summary(hetero1.fit, standardized=TRUE, fit.measures=TRUE)
modificationindices(hetero1.fit, sort=TRUE)

hetero.model2 <- "
        heterogeneity =~ soil + tri + tra_mean  + sea_mean  + number_biomes + area
        # mods
        tra_mean  ~~ area
"
hetero2.fit <- cfa(hetero.model2, data=dat_means_scaled)
fitmeasures(hetero1.fit, c("aic", "cfi", "rmsea", "srmr"))
fitmeasures(hetero2.fit, c("aic", "cfi", "rmsea", "srmr"))
anova(hetero1.fit, hetero2.fit)

modificationindices(hetero2.fit, sort=TRUE)
hetero.model3 <- "
        heterogeneity =~ soil + tri + tra_mean  + sea_mean  + number_biomes + area
        # mods
        tra_mean  ~~ area
        tri ~~ tra_mean 
"
hetero3.fit <- cfa(hetero.model3, data=dat_means_scaled)
fitmeasures(hetero2.fit, c("aic", "cfi", "rmsea", "srmr"))
fitmeasures(hetero3.fit, c("aic", "cfi", "tli", "rmsea", "srmr"))
anova(hetero2.fit, hetero3.fit)

# removing tra ?
hetero.model_no_tra <- "
        heterogeneity =~ tri + sea_mean  + number_biomes + area + soil
        tri ~~ area
"
hetero_no_tra.fit <- cfa(hetero.model_no_tra, data=dat_means_scaled)
summary(hetero_no_tra.fit, standardized=TRUE, fit.measures=TRUE)
parameterestimates(hetero_no_tra.fit, boot.ci.type = "bca.simple", standardized = TRUE)
modificationindices(hetero_no_tra.fit, sort=TRUE)
fitmeasures(hetero_no_tra.fit, c("aic", "cfi", "tli", "rmsea", "srmr"))

# Result: tra  as manifest variable is not predicted well by heterogeneity latent variable,
# the model performs better without it.
semPaths(object = hetero_no_tra.fit, layout = "tree", rotation = 2,
         whatLabels="std", label.cex = 1, edge.label.cex = 1, what="std")

# comments on the manifest variable order:
# by default the scaling of the latent variable is achieved by fixing the loading of the first indicator (manifest variable) for a certain latent variable to the value of 1. The order of the manifest variables is relevant only for fixing one loading per factor. Loading is fixed to 1.
# When I put soil first, it becomes non-significant. If I put it last, it`s fine and significant. don`T really understand how to interpret this yet


hetero.model4 <- "
        heterogeneity =~ elev_sd+ tri+ soil+ tra_sd+ sea_sd+
        tra_mean+ sea_mean+ pre_sd+ pet_sd+ mat_sd + number_biomes
        # mods
        elev_sd ~~ mat_sd
        tra_mean ~~ mat_sd
        sea_sd ~~ sea_mean
        tra_sd ~~tra_mean
        tra_sd ~~ mat_sd
        elev_sd ~~ tri
        tri ~~ mat_sd
        tra_mean ~~ pre_sd
        sea_mean ~~ pre_sd
"
hetero4.fit <- cfa(hetero.model4, data=dat_scaled)
fitmeasures(hetero4.fit, c("aic", "cfi", "tli", "rmsea", "srmr"))
fitmeasures(hetero3.fit, c("aic", "cfi", "tli", "rmsea", "srmr"))
summary(hetero4.fit, standardized=TRUE, fit.measures=TRUE)

modificationindices(hetero4.fit, sort=TRUE)
semPaths(hetero4.fit, edge.label.cex = 1, what="std")
semPaths(hetero3.fit, edge.label.cex = 1, what="std")

# this model includes all variables that make sense to put into heterogeneity, but it
# fits worse than the original and has plenty of correlation
# both have in common that tri and tra_mean are badly represented by heterogeneity, suggesting that they should not be included here but as 

#Try splitting up latent variable:
# DOES NOT WORK
# hetero.model5 <- "
#         heterogeneity =~ elev_sd+ soil+ tra_sd+ sea_sd+
#         sea_mean+ pet_sd+ mat_sd + number_biomes
#         hetero2 =~ tri + tra_mean + pre_sd
#         # mods
# "
# hetero5.fit <- cfa(hetero.model5, data=dat_scaled, check.gradient=FALSE)
# fitmeasures(hetero5.fit, c("aic", "cfi", "tli", "rmsea", "srmr"))
# fitmeasures(hetero3.fit, c("aic", "cfi", "tli", "rmsea", "srmr"))
# summary(hetero4.fit, standardized=TRUE, fit.measures=TRUE)


# adding heterogeneity as a latent variable
mod3 <- "# latent variable
         heterogeneity =~ number_biomes + sea_mean  + tri + area + soil
         tri ~~ area
        # regressions
          sr ~ heterogeneity + area + mrd
          mrd ~ heterogeneity
        # mods
        sea_mean  ~~ sr"

set.seed(1234)
fsem3 <- sem(mod3, data = dat_means_scaled)#, se = "bootstrap", bootstrap = 100) 
summary(fsem3, standardized = TRUE, fit.measures=TRUE)
# this model is really good!

parameterestimates(fsem3, boot.ci.type = "bca.simple", standardized = TRUE)
fitmeasures(fsem3, c("aic", "cfi", "tli", "rmsea", "srmr"))
modificationindices(fsem3, sort = TRUE)
semPaths(object = fsem3, layout = "tree", rotation = 2,
         whatLabels="std", label.cex = 1, edge.label.cex = 1, what="std")



#### climate as latent variable? ####
climate.model1 <- "
        climate =~ pet_mean  + mat_mean  + pre_mean 
"
climate1.fit <- cfa(climate.model1, data=dat_means_scaled)
summary(climate1.fit, standardized=TRUE, fit.measures=TRUE)
modificationindices(climate1.fit)

# This model is not right, something`s off. Also the inital correlation are low between 2 variables



# adding mat + precipitation (wolfs former model):
mod.wolf <- "# latent variable
          heterogeneity =~ number_biomes + sea_mean  + tri + area + soil
          tri ~~ area
        # regressions
          sr ~ heterogeneity +  area + mrd + mat_mean + pre_mean   
          mrd ~ heterogeneity + mat_mean  + pre_mean  
        # mods
          sea_mean  ~~ sr"
set.seed(1234)
fsem.wolf <- sem(mod.wolf, data = dat_means_scaled)#, se = "bootstrap", bootstrap = 100) 
summary(fsem.wolf, standardized = TRUE, fit.measures=TRUE)
fitmeasures(fsem.wolf, c("aic", "cfi", "tli", "rmsea", "srmr"))
# bad model fit

# Examine modification indices 
modificationindices(fsem.wolf, sort = TRUE)

mod.wolf2 <- "# latent variable
          heterogeneity =~ number_biomes + sea  + tri + area + soil
          tri ~~ area
        # regressions
          sr ~ heterogeneity +  area + mrd + mat  + pre
          mrd ~ heterogeneity + mat  + pre  
        # mods
          sea  ~~ sr
         area ~~ mrd
"
fsem.wolf2 <- sem(mod.wolf2, data = dat_means_scaled)#, se = "bootstrap", bootstrap = 100) 
fitmeasures(fsem.wolf2, c("aic", "cfi", "tli", "rmsea", "srmr"))
anova(fsem.wolf, fsem.wolf2)


mod.no_latent <-"
        # regressions
          sr ~ number_biomes + sea  + tri + area + soil + mat  + pre + pet + tra + elev
          mrd ~ number_biomes + sea  + tri + area + soil  + mat  + pre   + pet + tra + elev
        # mods
"
fit.no_latent <- sem(mod.no_latent, data = dat_means_scaled)#, se = "bootstrap", bootstrap = 100) 
summary(fit.no_latent, standardized = TRUE, fit.measures=TRUE)
semPaths(object = fit.no_latent, layout = "tree", rotation = 2,
         whatLabels="std", label.cex = 1, edge.label.cex = 1, what="std")

mod.latent <-"
          heterogeneity =~ number_biomes + sea  + tri + area + soil
          tri ~~ area
        # regressions
          sr ~  heterogeneity + mrd + mat  + pre + pet + tra + elev
          mrd ~ heterogeneity  + mat  + pre   + pet + tra + elev
        # mods
"
fit.latent <- sem(mod.latent, data = dat_means_scaled)
summary(fit.latent, standardized = TRUE, fit.measures=TRUE)
semPaths(object = fit.latent, layout = "tree", rotation = 2,
         whatLabels="std", label.cex = 1, edge.label.cex = 1, what="std")
anova(fit.no_latent, fit.latent)


mod.no_latent2 <-"
        # regressions
          number_biomes ~ area
          soil ~ area
          sr ~ number_biomes + sea  + tri + area + soil + mat  + pre + pet + tra + elev
          mrd ~ number_biomes + sea  + tri + area + soil  + mat  + pre   + pet + tra + elev
        # mods
"

fit.no_latent2 <- sem(mod.no_latent2, data = dat_means_scaled)#, se = "bootstrap", bootstrap = 100) 
summary(fit.no_latent2, standardized = TRUE, fit.measures=TRUE)
semPaths(object = fit.no_latent2, layout = "tree", rotation = 2,
         whatLabels="std", label.cex = 1, edge.label.cex = 1, what="std")


fitmeasures(fit.no_latent2, c("aic", "cfi", "tli", "rmsea", "srmr"))



# MEAN FULL MODEL ###########################################################################
mod.first_idea <-"
          spatial =~ tri + soil + number_biomes
          temporal =~ tra + sea
        # regressions
          spatial ~ area
          sr ~ spatial + temporal + mat  + pre + pet + elev + mrd
          mrd ~ spatial + temporal + mat  + pre + pet + elev
        # mods
"
fit.first_idea <- sem(mod.first_idea, data = dat_means_scaled)#, se = "bootstrap", bootstrap = 100) 
summary(fit.first_idea, standardized = TRUE, fit.measures=TRUE)
fitmeasures(fit.first_idea, c("aic", "cfi", "tli", "rmsea", "srmr"))
# bad fit
## temporal is not a great latent variable since the loadings are <0.3
## spatial is a strong latent variable
## temporal, precipitation and mrd have no significant direct influence on sr
## direct influeces on mrd: spatial,precipitation and elevation
modificationindices(fit.first_idea, sort = TRUE)

# remove temporal as latent variable
mod.first_idea2 <-"
          spatial =~ tri + soil + number_biomes
          #temporal =~ tra + sea
        # regressions
          spatial ~ area
          sr ~ spatial + mat  + pre + pet + elev + mrd + tra + sea
          mrd ~ spatial + mat  + pre + pet + elev + tra + sea
"
fit.first_idea2 <- sem(mod.first_idea2, data = dat_means_scaled)#, se = "bootstrap", bootstrap = 100) 
summary(fit.first_idea2, standardized = TRUE, fit.measures=TRUE)
fitmeasures(fit.first_idea2, c("aic", "cfi", "tli", "rmsea", "srmr"))
anova(fit.first_idea, fit.first_idea2)
# removing temporal as latent variable helps a lot
modificationindices(fit.first_idea2, sort = TRUE)

mod.first_idea3 <-"
          spatial =~ tri + soil + number_biomes
          #temporal =~ tra + sea
        # regressions
          spatial ~ area
          sr ~ spatial + mat  + pre + pet + elev + mrd + tra + sea
          mrd ~ spatial + mat  + pre + pet + elev + tra + sea
        # mod
          spatial ~ elev
"
fit.first_idea3 <- sem(mod.first_idea3, data = dat_means_scaled)
fitmeasures(fit.first_idea3, c("aic", "cfi", "tli", "rmsea", "srmr"))
anova(fit.first_idea2, fit.first_idea3) # adding influence of elevation on spatial helps the model
modificationindices(fit.first_idea3, sort = TRUE) # further suggested modifications dont make sense ecologically
summary(fit.first_idea3, standardized = TRUE, fit.measures=TRUE)

# Remove tra_mean from sr regression and move to spatial
mod.first_idea4 <-"
          spatial =~ tri + soil + number_biomes
          #temporal =~ tra + sea
        # regressions
          spatial ~ area
          sr ~ spatial + mat  + pre + pet + elev + mrd  + sea
          mrd ~ spatial + mat  + pre + pet + elev + tra + sea
        # mod
          spatial ~ elev + tra
"
fit.first_idea4 <- sem(mod.first_idea4, data = dat_means_scaled)
fitmeasures(fit.first_idea4, c("aic", "cfi", "tli", "rmsea", "srmr"))
anova(fit.first_idea3, fit.first_idea4) # adding influence of elevation on spatial helps the model
modificationindices(fit.first_idea4, sort = TRUE) # further suggested modifications dont make sense ecologically
summary(fit.first_idea4, standardized = TRUE, fit.measures=TRUE)



# SD FULL MODEL ########################################################################
#  # model does not converge, besides correlations of variables with sr /mrd actually being more meaningful.
# adding mat + precipitation (wolfs former model):
mod.wolf_sd <- "# latent variable
          heterogeneity =~ number_biomes + sea_mean  + tri + area + soil
          tri ~~ area
        # regressions
          sr ~ heterogeneity +  area + mrd + mat_sd + pre_sd  
          mrd ~ heterogeneity + mat_sd  + pre_sd
        # mods
          sea_mean  ~~ sr"
fsem.wolf_sd <- sem(mod.wolf_sd, data = dat_scaled) 
summary(fsem.wolf_sd, standardized = TRUE, fit.measures=TRUE)
fitmeasures(fsem.wolf_sd, c("aic", "cfi", "tli", "rmsea", "srmr"))
# model fit is actually worse although correlation with mat_sd and pre_sd is higher.




mod.first_idea_sd <-"
          spatial =~ tri + soil + number_biomes
          temporal =~ tra_sd + sea_sd
          
        # regressions
          spatial ~ area
          sr ~ spatial + temporal + mat_sd  + pre_sd + pet_sd + elev_sd + mrd
          mrd ~ spatial + temporal + mat_sd  + pre_sd + pet_sd + elev_sd
        # mods
"
fit.first_idea_sd <- sem(mod.first_idea_sd, data = dat_scaled)
summary(fit.first_idea_sd, standardized = TRUE, fit.measures=TRUE)
fitmeasures(fit.first_idea_sd, c("aic", "cfi", "tli", "rmsea", "srmr"))

# some SD trials
mod3_sd <- "# latent variable
         heterogeneity =~ number_biomes  + area + soil + sea_mean
         tri ~~ area
        # regressions
          sr ~ heterogeneity + area + mrd
          mrd ~ heterogeneity
        # mods
     sea_mean  ~~ sr"
fsem3_sd <- sem(mod3_sd, data = dat_scaled)#, se = "bootstrap", bootstrap = 100) 
fitmeasures(fsem3_sd, c("aic", "cfi", "tli", "rmsea", "srmr"))



# AREa predicting heterogeneity, not being part of it
mod4 <- "# latent variable
         heterogeneity =~ number_biomes + sea_mean  + tri + soil
         tri ~~ area
        # regressions
          sr ~ heterogeneity + area + mrd
          mrd ~ heterogeneity
          heterogeneity ~ area"

set.seed(1234)
fsem4 <- sem(mod4, data = dat_means_scaled) 
summary(fsem4, standardized = TRUE, fit.measures=TRUE)
# this model is really good

parameterestimates(fsem4, boot.ci.type = "bca.simple", standardized = TRUE)
fitmeasures(fsem4, c("aic", "cfi", "tli", "rmsea", "srmr"))

semPaths(fsem3, what="std", edge.label.cex=1)
semPaths(fsem4, what="std", edge.label.cex=1)
# Does not change the model, but area as part of heterogeneity is supported slightly better




#### basing the heterogeneity variable on all possible manifest ones   
mod3.2 <- "# latent variable
         heterogeneity =~ elev_sd+ tri+ soil+ tra_sd+ sea_sd+
          tra_mean+ sea_mean+ pre_sd+ pet_sd+ mat_sd + number_biomes
          # mods
          elev_sd ~~ mat_sd
          tra_mean ~~ mat_sd
          sea_sd ~~ sea_mean
          tra_sd ~~tra_mean
          tra_sd ~~ mat_sd
          elev_sd ~~ tri
          tri ~~ mat_sd
          tra_mean ~~ pre_sd
          sea_mean ~~ pre_sd
        # regressions
          sr ~ heterogeneity + area + mrd
          mrd ~ heterogeneity
          heterogeneity ~ area
"
fsem3.2 <- sem(mod3.2, data = na.omit(dat_scaled))
summary(fsem3.2, standardized = TRUE, fit.measures=TRUE)
# this model is not so great

fitmeasures(fsem3.2, c("aic", "cfi", "tli", "rmsea", "srmr"))
modificationindices(fsem3.2, sort = TRUE)
semPaths(object = fsem3.2, layout = "tree", rotation = 1,
         whatLabels="std", label.cex = 1, edge.label.cex = 1, what="std")

# figure out removing the tri + tra_mean from heterogeneity and putting as direct manifest variables



# model with theoretically informed heterogeneity variable and GBM results for SR and MRD
# (chosing most important variables from the GBM while separating them into manifest + latent)
mod_gbm <- "# latent variable
         heterogeneity =~ area + elev_sd + mat_sd + number_biomes + pre_sd + sea_sd + soil + tra_sd 
        # regressions
          sr ~ heterogeneity + area + mrd + '(sub)tropical moist broadleaf forest' + mangroves 
            + '(sub)tropical dry broadleaf forest' + 'montane grasslands and shrublands'
          mrd ~ heterogeneity + elev_mean + pre_mean + sea_mean + tra_mean
        # mods
          elev_sd ~~        mat_sd
          mat_sd ~~        tra_sd
          heterogeneity ~ elev_mean
"

fsem_gbm <- sem(mod_gbm, data = dat_scaled) 
#summary(fsem_gbm, standardized = TRUE, fit.measures=TRUE)
# this model is really bad!

fitmeasures(fsem_gbm, c("aic", "cfi", "tli", "rmsea", "srmr"))
modificationindices(fsem_gbm, sort = TRUE)
semPaths(object = fsem_gbm, layout = "tree", rotation = 2,
         whatLabels="std", label.cex = 1, edge.label.cex = 1, what="std")




# reduced dataset, following variable importance from GBMs ####
dat_red <- readRDS("dat_red.rds")
ztrans <- function(x){(x - mean(x)) / sd(x)}
dat_red <- apply(dat_red, 2, ztrans)
dat_red <- as.data.frame(dat_red)
names(dat_red)[18:29] <- paste0("biome",1:12)
  
mod_red_simple <-"
        # regressions
          sr ~ soil + 
          biome1 +
          pre_sd + 
          area +
          biome2 + 
          pre_mean + 
          elev_sd + 
          biome9 + 
          tra_sd + 
          mrd + 
          sea_mean+
          mat_mean+
          tra_mean+
          sea_sd+
          elev_mean+
          mat_sd+
          pet_sd
           mrd ~ soil + 
          biome1 +
          pre_sd + 
          area +
          biome2 + 
          pre_mean + 
          elev_sd + 
          biome9 + 
          tra_sd + 
#          mrd + 
          sea_mean+
          mat_mean+
          tra_mean+
          sea_sd+
          elev_mean+
          mat_sd+
          pet_sd
          # mods
#          pet_mean ~~ mat_mean
#	        pre_mean ~~ '(sub)tropical moist broadleaf forest'
#	        mat_sd ~~ elev_sd
#	        tra_sd ~~ area
"
# mod_red_simple <- "# latent variable
#          heterogeneity =~ number_biomes + sea_mean  + tri + area + soil
#          tri ~~ area
#         # regressions
#           sr ~ heterogeneity + area + mrd
#           mrd ~ heterogeneity
#         # mods
#         sea_mean  ~~ sr"

fsem_red_simp <- sem(mod_red_simple, data = dat_red)#, se = "bootstrap", bootstrap = 100) 
summary(fsem_red_simp, standardized = TRUE, fit.measures=TRUE)
# this model is weird

fitmeasures(fsem_red_simp, c("aic", "cfi", "tli", "rmsea", "srmr"))
modificationindices(fsem_red_simp, sort = TRUE)
semPaths(object = fsem_red_simp, layout = "tree", rotation = 2,
         whatLabels="std", label.cex = 1, edge.label.cex = 1, what="std")

mod7 <- "# latent variable
         heterogeneity =~ area + soil + pre_sd  + elev_sd + tra_mean
        # regressions
          sr ~ heterogeneity + area + mrd
          mrd ~ heterogeneity
        # mods
#        sea_mean  ~~ sr
"
fsem_7 <- sem(mod7, data = dat_red)#, se = "bootstrap", bootstrap = 100) 
summary(fsem_7, standardized = TRUE, fit.measures=TRUE)
