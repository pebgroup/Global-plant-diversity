# Other model tests
# produces data_for_SEM.rds
# next script: scaling_distribution_multicol-test.R

rm(list=ls())
library(sf)
library(caret)
library(gbm)
library(ggsci)
library(cowplot) 
library(ggcorrplot)
library(tidyverse)
library(data.table)
library(ppsr)
source("scripts/my_corrplot_func.R")

#### DATA PREP ##################################################################################################


shp <- readRDS("data/shp_object_fin_analysis.RDS")
org <- st_drop_geometry(shp)

# remove not needed columns 
# --> remove predomiant biomes as well as its a rougher and duplicate assessment compared to the biome percentages. 
#  --> mangroves and flooded grasslands are excluded because they are hierarchical very fine constructs compared to other biome types and therefore pretty rare
# --> remove elev_mean, that`s just bycatch from getting measn for other variables, we were interested in topgraphy structure, so SD or TRI

dat <- org[,-grep("trops|_n|lng|lat|z_score|soil_even|soil_simp|mangroves|flooded|predom|elev_mean|MRD.sd|MDR.sd", names(org))]
rownames(dat) <- dat$LEVEL_3_CO
dat <- dat[,c(grep("sr", names(dat)):ncol(dat))]

# two datasets: complete, resulting in less data (310) because of more NAs in sd. means only, resulting in 334 data points
# not working with the one without sds, uncomment
# 
# dat_means <- dat[,-grep("lng|lat_bin|z_score|trops|binary|_n|_sd|_even|_simp", names(dat))]
# dat_means <- na.omit(dat_means)
# 
# # standardize all variables with z-transformation to avoid variance issues - only necessary for SEM later
# ztrans <- function(x){(x - mean(x)) / sd(x)}
# 
# dat_means_scaled <- apply(dat_means[-grep("predominant_biome", names(dat_means))], 2, ztrans)
# dat_means_scaled <- as.data.frame(dat_means_scaled)
# dat_means_scaled$predominant_biome <- dat_means$predominant_biome
# 
# except paleoclimate from NA excloding and merge it back in again
temp <- dat[,grep("ano", names(dat))]
temp$level3 <- row.names(temp)
dat <- dat[,-grep("ano", names(dat))]
#dat_no.na <- dat[!is.na(c(paste0("dat$", names(dat)[-grep("ano", names(dat))]))), ] # keep paleoclim NAs
dat_no.na <- na.omit(dat)
dat_no.na$level3 <- row.names(dat_no.na)
dat_no.na <- merge(dat_no.na, temp, all.x=TRUE)
# dat_scaled <- apply(dat_no.na[-grep("predominant_biome", names(dat_no.na))], 2, ztrans)
# dat_scaled <- as.data.frame(dat_scaled)
# dat_scaled$predominant_biome <- dat_no.na$predominant_biome
dat_no.na





#### CORRELATION #############################################################################
# paleoclim-extrawurst
cor.test(dat_no.na$sr, dat_no.na$Mio_ano_AP)
cor.test(dat_no.na$sr, dat_no.na$Mio_ano_MAT)
cor.test(dat_no.na$sr, dat_no.na$Plio_ano_AP)
cor.test(dat_no.na$sr, dat_no.na$LGM_ano_MAT)
rownames(dat_no.na) <- dat_no.na$level3
dat_no.na <- dat_no.na[,-grep("ano|level3",names(dat_no.na))]
dat_no.na <- na.omit(dat_no.na)
dim(dat_no.na)
saveRDS(dat_no.na, "data/data_for_SEM.rds")

# rename sea to prs
names(dat_no.na) <- gsub("sea_","prs_", names(dat_no.na))
cor.dat_no.na<- cor(dat_no.na[,-grep("mdr|number|elev",names(dat_no.na))])
p.dat_no.na <- cor_pmat(dat_no.na[,-grep("mdr|number|elev",names(dat_no.na))])

my_corrplot(cor.dat_no.na, lab=TRUE, p.mat = p.dat_no.na, insig = "blank",
            tl.cex = 8, digits = 1, type = "lower", hc.order = FALSE,
            lab_size = 3, highlight = TRUE)+
  theme(axis.text.x = element_text(margin=margin(0,0,0,0)),  # Order: top, right, bottom, left
        axis.text.y = element_text(margin=margin(0,0,0,0)),
        plot.margin = margin(0, 0, 0, 0, "cm"))
ggsave("results/correlation_april2021.png", width=7, height=6, units = "in", dpi = 600)




# monitor pet_m + mat_m, mat_m + tra_m for multicollinearity
# chose tri over elev_sd, they monitor a similar thing and elev_sd has higher correlation with other variables
# chose soil over number of biomes, they monitor a similar thing and biomes types are covered in the proportions as well
# --> decision: remove elev_sd and number biomes from further analysis



#### GBM CONTROLS ###################################################################################################

# parameter tuning using caret
gbmControl <- trainControl(method = "repeatedcv", number = 10,
                           repeats = 3, savePredictions = "final",
                           returnResamp ="final")

gbmGrid <-  expand.grid(interaction.depth = c(1, 2, 3), 
                        n.trees = (1:10)*50, 
                        shrinkage = c(0.1, 0.01, 0.001),
                        n.minobsinnode = c(5, 10, 15))

# Full model SR ##################################################################################################
set.seed(519)
gbmFit4 <- train(sr ~ ., data = dat_no.na, method = "gbm", 
                 trControl = gbmControl, verbose=FALSE,
                 tuneGrid=gbmGrid)
gbmFit4$results[row.names(gbmFit4$bestTune),]
(full_sr <- ggplot(summary(gbmFit4), aes(x=reorder(var, rel.inf), y=rel.inf))+
  geom_bar(stat="identity")+
  coord_flip()+
  xlab("")+ylab("Relative influence"))
#ggsave(plot = full_sr, "results/varImp_full_gbm_SR_feb2021.png", width=4, height=4, units = "in", dpi = 600)


# Full model MRD ##################################################################################################
set.seed(973)
gbmFit4_mrd <- train(mrd ~ ., data = dat_no.na[,-grep("sr",names(dat_no.na))], method = "gbm", 
                     trControl = gbmControl, verbose=FALSE,
                     tuneGrid=gbmGrid)
gbmFit4_mrd$results[row.names(gbmFit4_mrd$bestTune),]
full_mrd <- ggplot(summary(gbmFit4_mrd), aes(x=reorder(var, rel.inf), y=rel.inf))+
  geom_bar(stat="identity")+
  coord_flip()+
  xlab("")+ylab("Relative influence")
#ggsave(plot=full_mrd, "results/varImp_full_gbm_mrd_feb2021.png", width=4, height=4, units = "in", dpi = 600)

plot_grid(full_sr, full_mrd, labels = c("A", "B"), ncol = 2)
#ggsave("results/varImp_full_gbm_feb2021.png", width=8, height=4, units = "in", dpi = 600)

# alternative back to back plot
df <- rbind(summary(gbmFit4), summary(gbmFit4_mrd))
df$response <- c(rep("SR", nrow(summary(gbmFit4))), rep("MRD", nrow(summary(gbmFit4_mrd))))
df$rel.inf[df$response=="MRD"] <- df$rel.inf[df$response=="MRD"]*-1
ggplot(df, aes(x=reorder(var, abs(rel.inf)), y= rel.inf, fill=response)) +
  facet_wrap(~ response, scales = "free_x") +
  geom_col() +
  coord_flip() +
  scale_y_continuous("Relative influence", expand = c(0, 0), labels = function(x) signif(abs(x), 3)) +
  xlab("")+
  theme(panel.spacing.x = unit(0, "mm"), legend.position = "none")+
  scale_fill_startrek()
#ggsave("results/varImp_full_gbm_btb_feb2021.png", width=4, height=4, units = "in", dpi = 600)


# model fit:
gbmFit4$results[row.names(gbmFit4$bestTune),]
gbmFit4_mrd$results[row.names(gbmFit4_mrd$bestTune),]



# model without number of biomes and elevationSD variables ########################
dat_no.na <- dat_no.na[,-grep("mdr", names(dat_no.na))]
set.seed(591)
gbmFit4_red <- train(sr ~ ., data = dat_no.na[,-grep("elev_sd|number", names(dat_no.na))], method = "gbm", 
                     trControl = gbmControl, verbose=FALSE,
                     tuneGrid=gbmGrid)
gbmFit4_red$results[row.names(gbmFit4_red$bestTune),]
(red_sr <- ggplot(summary(gbmFit4_red), aes(x=reorder(var, rel.inf), y=rel.inf))+
  geom_bar(stat="identity")+
  coord_flip()+
  xlab("")+ylab("Relative influence"))

set.seed(466) 
gbmFit4_red_mrd <- train(mrd ~ ., data = dat_no.na[,-grep("sr|elev_|number",names(dat_no.na))], method = "gbm", 
                     trControl = gbmControl, verbose=FALSE,
                     tuneGrid=gbmGrid)
gbmFit4_red_mrd$results[row.names(gbmFit4_red_mrd$bestTune),]
red_mrd <- ggplot(summary(gbmFit4_red_mrd), aes(x=reorder(var, rel.inf), y=rel.inf))+
  geom_bar(stat="identity")+
  coord_flip()+
  xlab("")+ylab("Relative influence")

plot_grid(red_sr, red_mrd, labels = c("A", "B"), ncol = 2)
ggsave("results/varImp_full_gbm_feb2021_no_elev+nobiomes.png", width=7, height=4, units = "in", dpi = 600)

# alternative back to back plot
df <- rbind(summary(gbmFit4_red), summary(gbmFit4_red_mrd))
df$response <- c(rep("SR", nrow(summary(gbmFit4_red))), rep("MRD", nrow(summary(gbmFit4_red_mrd))))
df$rel.inf[df$response=="MRD"] <- df$rel.inf[df$response=="MRD"]*-1
ggplot(df, aes(x=reorder(var, abs(rel.inf)), y= rel.inf, fill=response)) +
  facet_wrap(~ response, scales = "free_x") +
  geom_col() +
  coord_flip() +
  scale_y_continuous("Relative influence", expand = c(0, 0), labels = function(x) signif(abs(x), 3)) +
  xlab("")+
  theme(panel.spacing.x = unit(0, "mm"), legend.position = "none")+
  scale_fill_startrek()
ggsave("results/varImp_full_gbm_btb_feb2021_no_elev+nobiomes.png", width=4, height=4, units = "in", dpi = 600)


