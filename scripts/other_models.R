# Other model tests
# produces data_for_SEM.rds
# next script: scaling_distribution_multicol-test.R

rm(list=ls())
library(sf)
library(caret)
library(gbm)
library(ggsci)
#library(ggpubr) # this package has broom version issues
library(cowplot) 
library(ggcorrplot)
library(tidyverse)
library(data.table)
theme_set(theme_bw())
source("scripts/my_corrplot_func.R")

#### DATA PREP ##################################################################################################


shp <- readRDS("data/shp_object_fin_analysis.RDS")
org <- st_drop_geometry(shp)

# remove not needed columns 
# --> remove predomiant biomes as well as its a rougher and duplicate assessment compared to the biome percentages. 
#  --> mangroves and flooded grasslands are excluded because they are hierarchical very fine constructs compared to other biome types and therefore pretty rare
# --> remove elev_mean, that`s just bycatch from getting measn for other variables, we were interested in topgraphy structure, so SD or TRI
# --> remove soil eveness and simpson diversity because different approach
# --> remove mrd_z_scores since we do not scale MRD to a null model. 
# --> remove lng + lat + trops, never mind the gradients

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



# get catchy names
## get first col starting with (sub)
first <- grep("(sub)", names(dat_no.na))[1]
last <- grep("xeric", names(dat_no.na))
names(dat_no.na)[first:last] <- c("sub_trop_mbf", "sub_trop_dbf", "sub_trop_cf", "temp_bmf", "temp_cf",
                                    "boreal_f_taiga", "_sub_trop_gss", "temp_gss", "mont_gs", "tundra",
                                    "medit_fws", "deserts_x_shrub")



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

cor.dat_no.na<- cor(dat_no.na[,-grep("sr|mrd|mdr|number|elev",names(dat_no.na))])
p.dat_no.na <- cor_pmat(dat_no.na[,-grep("sr|mrd|mdr|number|elev",names(dat_no.na))])

my_corrplot(cor.dat_no.na, lab=TRUE, p.mat = p.dat_no.na, insig = "blank",
            tl.cex = 8, digits = 1, type = "lower", hc.order = FALSE, lab_size = 3, highlight = TRUE)+
  theme(axis.text.x = element_text(margin=margin(0,0,0,0)),  # Order: top, right, bottom, left
        axis.text.y = element_text(margin=margin(0,0,0,0)))
ggsave("results/correlation_march2021.png", width=6, height=5, units = "in", dpi = 600)

# correlations with SR and MRD
cor.dat_no.na<- cor(dat_no.na[,-grep("number|elev",names(dat_no.na))])
p.dat_no.na <- cor_pmat(dat_no.na[,-grep("number|elev",names(dat_no.na))])

my_corrplot(cor.dat_no.na[,c(1,2,3)], lab=TRUE, p.mat = p.dat_no.na[,c(1,2,3)], insig = "blank",
           tl.cex = 8, digits = 1, type = "lower", hc.order = FALSE, lab_size = 2.5)+
  theme(axis.text.x = element_text(margin=margin(0,0,0,0)),  # Order: top, right, bottom, left
        axis.text.y = element_text(margin=margin(0,0,0,0)))
ggsave("results/correlation_mrd_sr_march2021.png", width=6, height=2.1, units = "in", dpi = 600)


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
ggsave("results/varImp_full_gbm_feb2021_no_elev+nobiomes.png", width=8, height=4, units = "in", dpi = 600)

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





##### MODELS FOR SINGLE CATEGORIES #################################################

# CLIMATE MODELS ######

gbm_clim_sr <- train(sr ~ ., data = dat_no.na[,grep("sr|pet|mat|tra|pre_|sea", names(dat_no.na))], method = "gbm", 
                 trControl = gbmControl, verbose=FALSE,
                 tuneGrid=gbmGrid)
gbm_clim_sr$results[row.names(gbm_clim_sr$bestTune),]
gbm_clim_sr_plot <- ggplot(summary(gbm_clim_sr), aes(x=reorder(var, rel.inf), y=rel.inf))+
  geom_bar(stat="identity")+
  coord_flip()+
  xlab("")+ylab("Relative influence")

gbm_clim_mrd <- train(mrd ~ ., data = dat_no.na[,grep("mrd|pet|mat|tra|pre_|sea", names(dat_no.na))], method = "gbm", 
                     trControl = gbmControl, verbose=FALSE,
                     tuneGrid=gbmGrid)
gbm_clim_mrd$results[row.names(gbm_clim_mrd$bestTune),]
gbm_clim_mrd_plot <- ggplot(summary(gbm_clim_mrd), aes(x=reorder(var, rel.inf), y=rel.inf))+
  geom_bar(stat="identity")+
  coord_flip()+
  xlab("")+ylab("Relative influence")

plot_grid(gbm_clim_sr_plot, gbm_clim_mrd_plot, labels = c("A", "B"), ncol = 2)
ggsave("results/varImp_clim_gbm_feb2021.png", width=4, height=2, units = "in", dpi = 600)

# model fit:
gbm_clim_sr$results[row.names(gbm_clim_sr$bestTune),]
gbm_clim_mrd$results[row.names(gbm_clim_mrd$bestTune),]



# BIOME MODELS ######

gbm_bio_sr <- train(sr ~ ., data = dat_no.na[,grep("sr|sub|temp|taiga|tundra|mont|medit|desert", names(dat_no.na))], method = "gbm", 
                     trControl = gbmControl, verbose=FALSE,
                     tuneGrid=gbmGrid)
gbm_bio_sr$results[row.names(gbm_bio_sr$bestTune),]
gbm_bio_sr_plot <- ggplot(summary(gbm_bio_sr), aes(x=reorder(var, rel.inf), y=rel.inf))+
  geom_bar(stat="identity")+
  coord_flip()+
  xlab("")+ylab("Relative influence")

gbm_bio_mrd <- train(mrd ~ ., data = dat_no.na[,grep("mrd|sub|temp|taiga|tundra|mont|medit|desert", names(dat_no.na))], method = "gbm", 
                     trControl = gbmControl, verbose=FALSE,
                     tuneGrid=gbmGrid)
gbm_bio_mrd$results[row.names(gbm_bio_mrd$bestTune),]
gbm_bio_mrd_plot <- ggplot(summary(gbm_bio_mrd), aes(x=reorder(var, rel.inf), y=rel.inf))+
  geom_bar(stat="identity")+
  coord_flip()+
  xlab("")+ylab("Relative influence")

# model fit:
gbm_bio_sr$results[row.names(gbm_bio_sr$bestTune),]
gbm_bio_mrd$results[row.names(gbm_bio_mrd$bestTune),]



# terrain MODELS ######

gbm_geo_sr <- train(sr ~ ., data = dat_no.na[,grep("sr|area|soil|tri|elev", names(dat_no.na))], method = "gbm", 
                    trControl = gbmControl, verbose=FALSE,
                    tuneGrid=gbmGrid)
gbm_geo_sr$results[row.names(gbm_geo_sr$bestTune),]
gbm_geo_sr_plot <- ggplot(summary(gbm_geo_sr), aes(x=reorder(var, rel.inf), y=rel.inf))+
  geom_bar(stat="identity")+
  coord_flip()+
  xlab("")+ylab("Relative influence")

gbm_geo_mrd <- train(mrd ~ ., data = dat_no.na[,grep("mrd|area|soil|tri|elev", names(dat_no.na))], method = "gbm", 
                     trControl = gbmControl, verbose=FALSE,
                     tuneGrid=gbmGrid)
gbm_geo_mrd$results[row.names(gbm_geo_mrd$bestTune),]
gbm_geo_mrd_plot <- ggplot(summary(gbm_geo_mrd), aes(x=reorder(var, rel.inf), y=rel.inf))+
  geom_bar(stat="identity")+
  coord_flip()+
  xlab("")+ylab("Relative influence")

# model fit:
gbm_geo_sr$results[row.names(gbm_geo_sr$bestTune),]
gbm_geo_mrd$results[row.names(gbm_geo_mrd$bestTune),]


# PLOTs ####


# left_col <- plot_grid(gbm_clim_sr_plot, gbm_bio_sr_plot, gbm_geo_sr_plot, labels = c('A','C','E'),
#                       label_size = 10, ncol=1, label_y = c(0.9,1,1))
# right_col <- plot_grid(gbm_clim_mrd_plot, gbm_bio_mrd_plot, gbm_geo_mrd_plot, labels = c('B','D','F'),
#                        label_size = 10, ncol=1)
# plot_grid(left_col, right_col, labels = c('SR', 'MRD'), label_size = 12, ncol = 2, label_y = 1.01)

first_row <- plot_grid(gbm_clim_sr_plot, gbm_clim_mrd_plot, labels = c('SR','MRD'),
                      label_size = 10, ncol=2, label_x = c(0.5), label_y = 1.07,  label_fontface = "plain")
sec_row <- plot_grid(gbm_bio_sr_plot, gbm_bio_mrd_plot, labels = c('SR','MRD'),
                       label_size = 10, ncol=2, label_x = c(0.5), label_y = 1.07, label_fontface = "plain")
third_row <- plot_grid(gbm_geo_sr_plot, gbm_geo_mrd_plot, labels = c('SR','MRD'),
                     label_size = 10, ncol=2, label_x = c(0.5), label_y = 1.07, label_fontface = "plain")
#plot_grid(first_row, sec_row, third_row, labels = c('climate', 'biomes', 'aerial'), label_size = 12, ncol = 1)

# now add the title
title1 <- ggdraw() + 
  draw_label("A", fontface = 'bold',x = 0,hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))
title2 <- ggdraw() + 
  draw_label("B", fontface = 'bold',x = 0,hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))
title3 <- ggdraw() + 
  draw_label("C", fontface = 'bold',x = 0,hjust = 0) +
  theme(plot.margin = margin(0, 0, 0, 7))
plot_grid(title1, first_row, title2, sec_row, title3, third_row, ncol = 1, rel_heights = c(0.1, 1))
ggsave("results/varImp_subs_gbm_feb2021.png", width=8, height=8, units = "in", dpi = 600)

dfclim <- rbind(summary(gbm_clim_sr), summary(gbm_clim_mrd))
dfclim$response <- c(rep("SR", nrow(summary(gbm_clim_sr))), rep("MRD", nrow(summary(gbm_clim_mrd))))
dfclim$rel.inf[dfclim$response=="MRD"] <- dfclim$rel.inf[dfclim$response=="MRD"]*-1
clim_plot <- ggplot(dfclim, aes(x=reorder(var, rel.inf), y= rel.inf, fill=response)) +
  facet_wrap(~ response, scales = "free_x") +
  geom_col() +
  coord_flip() +
  scale_y_continuous("Relative influence", expand = c(0, 0), labels = function(x) signif(abs(x), 3)) +
  xlab("")+
  theme(panel.spacing.x = unit(0, "mm"), legend.position = "none")+
  scale_fill_startrek()

dfbio <- rbind(summary(gbm_bio_sr), summary(gbm_bio_mrd))
dfbio$response <- c(rep("SR", nrow(summary(gbm_bio_sr))), rep("MRD", nrow(summary(gbm_bio_mrd))))
dfbio$rel.inf[dfbio$response=="MRD"] <- dfbio$rel.inf[dfbio$response=="MRD"]*-1
bio_plot <- ggplot(dfbio, aes(x=reorder(var, abs(rel.inf)), y= rel.inf, fill=response)) +
  facet_wrap(~ response, scales = "free_x") +
  geom_col() +
  coord_flip() +
  scale_y_continuous("Relative influence", expand = c(0, 0), labels = function(x) signif(abs(x), 3)) +
  xlab("")+
  theme(panel.spacing.x = unit(0, "mm"), legend.position = "none")+
  scale_fill_startrek()

dfgeo <- rbind(summary(gbm_geo_sr), summary(gbm_geo_mrd))
dfgeo$response <- c(rep("SR", nrow(summary(gbm_geo_sr))), rep("MRD", nrow(summary(gbm_geo_mrd))))
dfgeo$rel.inf[dfgeo$response=="MRD"] <- dfgeo$rel.inf[dfgeo$response=="MRD"]*-1
geo_plot <- ggplot(dfgeo, aes(x=reorder(var, rel.inf), y= rel.inf, fill=response)) +
  facet_wrap(~ response, scales = "free_x") +
  geom_col() +
  coord_flip() +
  scale_y_continuous("Relative influence", expand = c(0, 0), labels = function(x) signif(abs(x), 3)) +
  xlab("")+
  theme(panel.spacing.x = unit(0, "mm"), legend.position = "none")+
  scale_fill_startrek()

plot_grid(clim_plot, bio_plot, geo_plot, ncol = 1, labels="AUTO", rel_heights = c(10,11,7.5))
ggsave("results/varImp_sub_gbm_btb_feb2021.png", width=4, height=6, units = "in", dpi = 600)

# plot_grid(gbm_clim_sr_plot, gbm_clim_mrd_plot, 
#           gbm_bio_sr_plot, gbm_bio_mrd_plot,
#           gbm_geo_sr_plot, gbm_geo_mrd_plot, labels = "AUTO", ncol = 2)

















































# TRASH
# 
# # FOLLOWING PART IS OLD STUFF ##################################################################################
# 
# ##  tri + elev_sd ##############################################################################################
# # correlation with SR and MRD
# cor.test(dat_no.na$sr, dat_no.na$elev_sd, method="s") #0.45
# cor.test(dat_no.na$sr, dat_no.na$tri, method="s") #0.15
# cor.test(dat_no.na$mrd, dat_no.na$elev_sd, method="s") #0.29
# cor.test(dat_no.na$mrd, dat_no.na$tri, method="s") # NS
# 
# cor.test(dat_no.na$elev_sd, dat_no.na$tri, method="s") # 0.71
# 
# # TRI in the full model has importance for MRD, less SR. Subset model supports this. 
# # elev_sd has less importance than TRI in the subsets, same for full MRD model, but no SR (TRI less important)
# 
# # SR TEST"
# gbmFit5_no_elev_sd <- train(sr ~ ., data = dat_no.na[,-grep("elev_sd",names(dat_no.na))], method = "gbm", 
#                             trControl = gbmControl, verbose=FALSE,
#                             tuneGrid=gbmGrid)
# gbmFit5_no_tri <- train(sr ~ ., data = dat_no.na[,-grep("tri",names(dat_no.na))], method = "gbm", 
#                         trControl = gbmControl, verbose=FALSE,
#                         tuneGrid=gbmGrid)
# 
# gbmFit4$results[row.names(gbmFit4$bestTune),]
# gbmFit5_no_tri$results[row.names(gbmFit5_no_tri$bestTune),] # better
# gbmFit5_no_elev_sd$results[row.names(gbmFit5_no_elev_sd$bestTune),] # worse
# # models are almost the same
# 
# # variable order
# ord <- as.data.frame(summary(gbmFit4))
# res2 <- as.data.frame(summary(gbmFit5_no_elev_sd))
# res3 <- as.data.frame(summary(gbmFit5_no_tri))
# ord_sum <- merge(ord, res2, by="var", all=TRUE)
# ord_sum <- merge(ord_sum, res3, by="var", all=TRUE)
# 
# # plot relative influence in each model instead of order?
# test <- ord_sum[,c(1:4)]
# names(test) <- c("var", "full_model", "no_elev_sd_model", "no_tri_model")
# test2 <- test %>% 
#   pivot_longer(cols=c(full_model, no_elev_sd_model, no_tri_model))
# test2$col <- "other"
# test2$col[test2$var=="tri"] <- "tri"
# test2$col[test2$var=="elev_sd"] <- "elev_sd"
# test2$name <- factor(test2$name, levels = c("no_elev_sd_model", "full_model", "no_tri_model"))
# ggplot(test2, aes(x=name, y=value, group=var, col=test2$col))+
#   geom_point()+
#   geom_line()+
#   scale_color_manual("Variable", values=c("blue", "grey", "red"))+
#   #geom_label(aes(label=var), size=3, col=test2$col)+
#   scale_y_log10("Relative influence")+
#   theme(legend.position = "right")
# 
# ggplot(temp, aes(x=name, y=value, group=var, col=temp$col))+ 
#   geom_line()+
#   geom_label(aes(label=var), size=3)+
#   scale_y_continuous("Variable rank")+
#   scale_color_manual(values=c("blue", "grey", "red"))+
#   facet_wrap(~ response)+
#   theme(legend.position = "none", 
#         panel.spacing.x = unit(0, "mm"))
# 
# 
# ####
# 
# 
# ord_sum$no_tri_model <- rank(ord_sum$rel.inf, ties.method = "last", na.last = "keep")
# ord_sum$full_model <- rank(ord_sum$rel.inf.x, ties.method = "last", na.last = "keep")
# ord_sum$no_elev_sd_model <- rank(ord_sum$rel.inf.y, ties.method = "last", na.last = "keep")
# 
# ord_sum2 <- ord_sum %>% 
#   pivot_longer(cols=c(full_model, no_elev_sd_model, no_tri_model))
# #ord_sum2$var2 <- abbreviate(ord_sum2$var, minlength = 7, named = F)
# ord_sum2$name <- factor(ord_sum2$name, levels = c("no_elev_sd_model", "full_model", "no_tri_model"))
# 
# ord_sum2$col <- "grey"
# ord_sum2$col[ord_sum2$var=="tri"] <- "red"
# ord_sum2$col[ord_sum2$var=="elev_sd"] <- "blue"
# ggplot(ord_sum2, aes(x=name, y=value, col=col))+
# #  geom_point()+
#   geom_line(aes(group=var), col=ord_sum2$col)+
#   geom_label(aes(label=var), size=3, col=ord_sum2$col)+
#   scale_y_continuous("Variable rank")+
#   theme(legend.position = "none")
# ggsave("results/tri_vs_elev_sd_SR_feb2021.png", width=6, height=5, units = "in", dpi = 600)
# 
# # MRD TEST"
# gbm_mrd_no_elev_sd <- train(mrd ~ ., data = dat_no.na[,-grep("elev_sd|sr",names(dat_no.na))], method = "gbm", 
#                             trControl = gbmControl, verbose=FALSE,
#                             tuneGrid=gbmGrid)
# gbm_mrd_no_tri <- train(mrd ~ ., data = dat_no.na[,-grep("tri|sr",names(dat_no.na))], method = "gbm", 
#                         trControl = gbmControl, verbose=FALSE,
#                         tuneGrid=gbmGrid)
# 
# gbmFit4_mrd$results[row.names(gbmFit4_mrd$bestTune),]
# gbm_mrd_no_elev_sd$results[row.names(gbm_mrd_no_elev_sd$bestTune),] # worse
# gbm_mrd_no_tri$results[row.names(gbm_mrd_no_tri$bestTune),] # better
# 
# # variable order
# ord <- as.data.frame(summary(gbmFit4_mrd))
# res2 <- as.data.frame(summary(gbm_mrd_no_elev_sd))
# res3 <- as.data.frame(summary(gbm_mrd_no_tri))
# ord_sum <- merge(ord, res2, by="var", all=TRUE)
# ord_sum <- merge(ord_sum, res3, by="var", all=TRUE)
# ord_sum$no_tri_model <- rank(ord_sum$rel.inf, ties.method = "last", na.last = "keep")
# ord_sum$full_model <- rank(ord_sum$rel.inf.x, ties.method = "last", na.last = "keep")
# ord_sum$no_elev_sd_model <- rank(ord_sum$rel.inf.y, ties.method = "last", na.last = "keep")
# 
# ord_sum2_mrd <- ord_sum %>% 
#   pivot_longer(cols=c(full_model, no_elev_sd_model, no_tri_model))
# #ord_sum2$var2 <- abbreviate(ord_sum2$var, minlength = 7, named = F)
# ord_sum2_mrd$name <- factor(ord_sum2_mrd$name, levels = c("no_elev_sd_model", "full_model", "no_tri_model"))
# 
# ord_sum2_mrd$col <- "grey"
# ord_sum2_mrd$col[ord_sum2_mrd$var=="tri"] <- "red"
# ord_sum2_mrd$col[ord_sum2_mrd$var=="elev_sd"] <- "blue"
# ggplot(ord_sum2_mrd, aes(x=name, y=value, col=col))+
#   #  geom_point()+
#   geom_line(aes(group=var), col=ord_sum2_mrd$col)+
#   geom_label(aes(label=var), size=3, col=ord_sum2_mrd$col)+
#   scale_y_continuous("Variable rank")+
#   theme(legend.position = "none")
# ggsave("results/tri_vs_elev_sd_MRD_feb2021.png", width=6, height=5, units = "in", dpi = 600)
# 
# # FACET PLOT
# ord_sum2$response <- "SR"
# ord_sum2_mrd$response <- "MRD"
# temp <- rbind(ord_sum2, ord_sum2_mrd)
# ggplot(temp, aes(x=name, y=value, group=var, col=temp$col))+ 
#   geom_line()+
#   geom_label(aes(label=var), size=3)+
#   scale_y_continuous("Variable rank")+
#   scale_color_manual(values=c("blue", "grey", "red"))+
#   facet_wrap(~ response)+
#   theme(legend.position = "none", 
#         panel.spacing.x = unit(0, "mm"))
# ggsave("results/tri_vs_elev_sd_feb2021.png", width=9, height=5, units = "in", dpi = 600)
# 
# 
# # decision: 
# ## - TRI small correlation with SR and no correlation with MRD
# ## - removing TRI: reduced RMSE and increased Rsquared compared to full model for both SR and MRD - thats for the GBM, does not say much about SEM
# ## does not seem to cause much trouble and can probably be combined in one latent variable
# 
# 
# 
# 
# 
# 
# ##  elev_mean + elev_sd ##############################################################################################
# # correlation with SR and MRD
# cor.test(dat_no.na$sr, dat_no.na$elev_sd, method="s") #0.45
# cor.test(dat_no.na$sr, dat_no.na$elev_mean, method="s") #0.34
# cor.test(dat_no.na$mrd, dat_no.na$elev_sd, method="s") #0.29
# cor.test(dat_no.na$mrd, dat_no.na$elev_mean, method="s") # 0.34
# 
# cor.test(dat_no.na$elev_sd, dat_no.na$elev_mean, method="s") # 0.83
# 
# 
# # SR TEST"
# gbmFit5_no_elev_mean <- train(sr ~ ., data = dat_no.na[,-grep("elev_mean",names(dat_no.na))], method = "gbm", 
#                             trControl = gbmControl, verbose=FALSE,
#                             tuneGrid=gbmGrid)
# 
# gbmFit4$results[row.names(gbmFit4$bestTune),]
# gbmFit5_no_elev_mean$results[row.names(gbmFit5_no_elev_mean$bestTune),] # better
# gbmFit5_no_elev_sd$results[row.names(gbmFit5_no_elev_sd$bestTune),] # worse
# # removing elev_mean improved model minimally
# 
# # variable order
# ord <- as.data.frame(summary(gbmFit4))
# res2 <- as.data.frame(summary(gbmFit5_no_elev_sd))
# res3 <- as.data.frame(summary(gbmFit5_no_elev_mean))
# ord_sum <- merge(ord, res2, by="var", all=TRUE)
# ord_sum <- merge(ord_sum, res3, by="var", all=TRUE)
# 
# # plot relative influence in each model instead of order?
# test <- ord_sum[,c(1:4)]
# names(test) <- c("var", "full_model", "no_elev_sd_model", "no_elev_mean_model")
# test2 <- test %>% 
#   pivot_longer(cols=c(full_model, no_elev_sd_model, no_elev_mean_model))
# test2$col <- "other"
# test2$col[test2$var=="elev_mean"] <- "elev_mean"
# test2$col[test2$var=="elev_sd"] <- "elev_sd"
# test2$name <- factor(test2$name, levels = c("no_elev_sd_model", "full_model", "no_elev_mean_model"))
# ggplot(test2, aes(x=name, y=value, group=var, col=test2$col))+
#   geom_point()+
#   geom_line()+
#   scale_color_manual("Variable", values=c("blue", "red", "grey"))+
#   #geom_label(aes(label=var), size=3, col=test2$col)+
#   scale_y_log10("Relative influence")+
#   theme(legend.position = "right")
# 
# 
# ord_sum$no_elev_mean_model <- rank(ord_sum$rel.inf, ties.method = "last", na.last = "keep")
# ord_sum$full_model <- rank(ord_sum$rel.inf.x, ties.method = "last", na.last = "keep")
# ord_sum$no_elev_sd_model <- rank(ord_sum$rel.inf.y, ties.method = "last", na.last = "keep")
# 
# ord_sum2 <- ord_sum %>% 
#   pivot_longer(cols=c(full_model, no_elev_sd_model, no_elev_mean_model))
# #ord_sum2$var2 <- abbreviate(ord_sum2$var, minlength = 7, named = F)
# ord_sum2$name <- factor(ord_sum2$name, levels = c("no_elev_sd_model", "full_model", "no_elev_mean_model"))
# 
# ord_sum2$col <- "grey"
# ord_sum2$col[ord_sum2$var=="elev_mean"] <- "red"
# ord_sum2$col[ord_sum2$var=="elev_sd"] <- "blue"
# ggplot(ord_sum2, aes(x=name, y=value, col=col))+
#   #  geom_point()+
#   geom_line(aes(group=var), col=ord_sum2$col)+
#   geom_label(aes(label=var), size=3, col=ord_sum2$col)+
#   scale_y_continuous("Variable rank")+
#   theme(legend.position = "none")
# #ggsave("results/elev_mean_vs_elev_sd_SR.png", width=6, height=5, units = "in", dpi = 600)
# 
# # MRD TEST"
# gbm_mrd_no_elev_mean <- train(mrd ~ ., data = dat_no.na[,-grep("elev_mean|sr",names(dat_no.na))], method = "gbm", 
#                         trControl = gbmControl, verbose=FALSE,
#                         tuneGrid=gbmGrid)
# 
# gbmFit4_mrd$results[row.names(gbmFit4_mrd$bestTune),]
# gbm_mrd_no_elev_sd$results[row.names(gbm_mrd_no_elev_sd$bestTune),] # worse
# gbm_mrd_no_elev_mean$results[row.names(gbm_mrd_no_elev_mean$bestTune),] # better
# 
# 
# # variable order
# ord <- as.data.frame(summary(gbmFit4_mrd))
# res2 <- as.data.frame(summary(gbm_mrd_no_elev_sd))
# res3 <- as.data.frame(summary(gbm_mrd_no_elev_mean))
# ord_sum <- merge(ord, res2, by="var", all=TRUE)
# ord_sum <- merge(ord_sum, res3, by="var", all=TRUE)
# ord_sum$no_elev_mean_model <- rank(ord_sum$rel.inf, ties.method = "last", na.last = "keep")
# ord_sum$full_model <- rank(ord_sum$rel.inf.x, ties.method = "last", na.last = "keep")
# ord_sum$no_elev_sd_model <- rank(ord_sum$rel.inf.y, ties.method = "last", na.last = "keep")
# 
# ord_sum2_mrd <- ord_sum %>% 
#   pivot_longer(cols=c(full_model, no_elev_sd_model, no_elev_mean_model))
# #ord_sum2$var2 <- abbreviate(ord_sum2$var, minlength = 7, named = F)
# ord_sum2_mrd$name <- factor(ord_sum2_mrd$name, levels = c("no_elev_sd_model", "full_model", "no_elev_mean_model"))
# 
# ord_sum2_mrd$col <- "grey"
# ord_sum2_mrd$col[ord_sum2_mrd$var=="elev_mean"] <- "red"
# ord_sum2_mrd$col[ord_sum2_mrd$var=="elev_sd"] <- "blue"
# ggplot(ord_sum2_mrd, aes(x=name, y=value, col=col))+
#   #  geom_point()+
#   geom_line(aes(group=var), col=ord_sum2_mrd$col)+
#   geom_label(aes(label=var), size=3, col=ord_sum2_mrd$col)+
#   scale_y_continuous("Variable rank")+
#   theme(legend.position = "none")
# #ggsave("results/elev_mean_vs_elev_sd_MRD.png", width=6, height=5, units = "in", dpi = 600)
# 
# # FACET PLOT
# ord_sum2$response <- "SR"
# ord_sum2_mrd$response <- "MRD"
# temp <- rbind(ord_sum2, ord_sum2_mrd)
# ggplot(temp, aes(x=name, y=value, group=var, col=temp$col))+ 
#   geom_line()+
#   geom_label(aes(label=var), size=3)+
#   scale_y_continuous("Variable rank")+
#   scale_color_manual(values=c("blue", "grey", "red"))+
#   facet_wrap(~ response)+
#   theme(legend.position = "none", 
#         panel.spacing.x = unit(0, "mm"))
# ggsave("results/elev_mean_vs_elev_sd.png", width=9, height=5, units = "in", dpi = 600)
# 
# 
# # decision: 
# # elev_mean and elev_tri are less correlated than elev_mean and sd, so removing elev_sd would kil two birds with one stone.
# ## - elev_mean small correlation with SR and no correlation with MRD
# ## - removing elev_mean: reduced RMSE and increased Rsqured compared to full model for both SR and MRD - thats for the GBM, does not say much about SEM
# 
# 
# 
# 
# 
# 
# 
# 
# ##  mat_sd + elev_sd ##########################################################################################
# # correlation with SR and MRD
# cor.test(dat_no.na$sr, dat_no.na$elev_sd, method="s") #0.45
# cor.test(dat_no.na$sr, dat_no.na$mat_sd, method="s") #0.43
# cor.test(dat_no.na$mrd, dat_no.na$elev_sd, method="s") #0.29
# cor.test(dat_no.na$mrd, dat_no.na$mat_sd, method="s") # 0.29
# 
# cor.test(dat_no.na$elev_sd, dat_no.na$mat_sd, method="s") # 0.69
# 
# 
# # SR TEST"
# gbmFit5_no_mat_sd <- train(sr ~ ., data = dat_no.na[,-grep("mat_sd",names(dat_no.na))], method = "gbm", 
#                               trControl = gbmControl, verbose=FALSE,
#                               tuneGrid=gbmGrid)
# 
# gbmFit4$results[row.names(gbmFit4$bestTune),]
# gbmFit5_no_mat_sd$results[row.names(gbmFit5_no_mat_sd$bestTune),] # better
# gbmFit5_no_elev_sd$results[row.names(gbmFit5_no_elev_sd$bestTune),] # worse
# 
# # variable order
# ord <- as.data.frame(summary(gbmFit4))
# res2 <- as.data.frame(summary(gbmFit5_no_elev_sd))
# res3 <- as.data.frame(summary(gbmFit5_no_mat_sd))
# ord_sum <- merge(ord, res2, by="var", all=TRUE)
# ord_sum <- merge(ord_sum, res3, by="var", all=TRUE)
# 
# # plot relative influence in each model instead of order?
# test <- ord_sum[,c(1:4)]
# names(test) <- c("var", "full_model", "no_elev_sd_model", "no_mat_sd_model")
# test2 <- test %>% 
#   pivot_longer(cols=c(full_model, no_elev_sd_model, no_mat_sd_model))
# test2$col <- "other"
# test2$col[test2$var=="mat_sd"] <- "mat_sd"
# test2$col[test2$var=="elev_sd"] <- "elev_sd"
# test2$name <- factor(test2$name, levels = c("no_elev_sd_model", "full_model", "no_mat_sd_model"))
# ggplot(test2, aes(x=name, y=value, group=var, col=test2$col))+
#   geom_point()+
#   geom_line()+
#   scale_color_manual("Variable", values=c("blue", "red", "grey"))+
#   #geom_label(aes(label=var), size=3, col=test2$col)+
#   scale_y_log10("Relative influence")+
#   theme(legend.position = "right")
# 
# 
# ord_sum$no_mat_sd_model <- rank(ord_sum$rel.inf, ties.method = "last", na.last = "keep")
# ord_sum$full_model <- rank(ord_sum$rel.inf.x, ties.method = "last", na.last = "keep")
# ord_sum$no_elev_sd_model <- rank(ord_sum$rel.inf.y, ties.method = "last", na.last = "keep")
# 
# ord_sum2 <- ord_sum %>% 
#   pivot_longer(cols=c(full_model, no_elev_sd_model, no_mat_sd_model))
# #ord_sum2$var2 <- abbreviate(ord_sum2$var, minlength = 7, named = F)
# ord_sum2$name <- factor(ord_sum2$name, levels = c("no_elev_sd_model", "full_model", "no_mat_sd_model"))
# 
# ord_sum2$col <- "grey"
# ord_sum2$col[ord_sum2$var=="mat_sd"] <- "red"
# ord_sum2$col[ord_sum2$var=="elev_sd"] <- "blue"
# ggplot(ord_sum2, aes(x=name, y=value, col=col))+
#   #  geom_point()+
#   geom_line(aes(group=var), col=ord_sum2$col)+
#   geom_label(aes(label=var), size=3, col=ord_sum2$col)+
#   scale_y_continuous("Variable rank")+
#   theme(legend.position = "none")
# #ggsave("results/mat_sd_vs_elev_sd_SR.png", width=6, height=5, units = "in", dpi = 600)
# 
# # MRD TEST"
# gbm_mrd_no_mat_sd <- train(mrd ~ ., data = dat_no.na[,-grep("mat_sd|sr",names(dat_no.na))], method = "gbm", 
#                               trControl = gbmControl, verbose=FALSE,
#                               tuneGrid=gbmGrid)
# 
# gbmFit4_mrd$results[row.names(gbmFit4_mrd$bestTune),]
# gbm_mrd_no_elev_sd$results[row.names(gbm_mrd_no_elev_sd$bestTune),] # worse
# gbm_mrd_no_mat_sd$results[row.names(gbm_mrd_no_mat_sd$bestTune),] # better
# 
# 
# # variable order
# ord <- as.data.frame(summary(gbmFit4_mrd))
# res2 <- as.data.frame(summary(gbm_mrd_no_elev_sd))
# res3 <- as.data.frame(summary(gbm_mrd_no_mat_sd))
# ord_sum <- merge(ord, res2, by="var", all=TRUE)
# ord_sum <- merge(ord_sum, res3, by="var", all=TRUE)
# ord_sum$no_mat_sd_model <- rank(ord_sum$rel.inf, ties.method = "last", na.last = "keep")
# ord_sum$full_model <- rank(ord_sum$rel.inf.x, ties.method = "last", na.last = "keep")
# ord_sum$no_elev_sd_model <- rank(ord_sum$rel.inf.y, ties.method = "last", na.last = "keep")
# 
# ord_sum2_mrd <- ord_sum %>% 
#   pivot_longer(cols=c(full_model, no_elev_sd_model, no_mat_sd_model))
# #ord_sum2$var2 <- abbreviate(ord_sum2$var, minlength = 7, named = F)
# ord_sum2_mrd$name <- factor(ord_sum2_mrd$name, levels = c("no_elev_sd_model", "full_model", "no_mat_sd_model"))
# 
# ord_sum2_mrd$col <- "grey"
# ord_sum2_mrd$col[ord_sum2_mrd$var=="mat_sd"] <- "red"
# ord_sum2_mrd$col[ord_sum2_mrd$var=="elev_sd"] <- "blue"
# ggplot(ord_sum2_mrd, aes(x=name, y=value, col=col))+
#   #  geom_point()+
#   geom_line(aes(group=var), col=ord_sum2_mrd$col)+
#   geom_label(aes(label=var), size=3, col=ord_sum2_mrd$col)+
#   scale_y_continuous("Variable rank")+
#   theme(legend.position = "none")
# #ggsave("results/mat_sd_vs_elev_sd_MRD.png", width=6, height=5, units = "in", dpi = 600)
# 
# # FACET PLOT
# ord_sum2$response <- "SR"
# ord_sum2_mrd$response <- "MRD"
# temp <- rbind(ord_sum2, ord_sum2_mrd)
# ggplot(temp, aes(x=name, y=value, group=var, col=temp$col))+ 
#   geom_line()+
#   geom_label(aes(label=var), size=3)+
#   scale_y_continuous("Variable rank")+
#   scale_color_manual(values=c("blue", "grey", "red"))+
#   facet_wrap(~ response)+
#   theme(legend.position = "none", 
#         panel.spacing.x = unit(0, "mm"))
# ggsave("results/mat_sd_vs_elev_sd.png", width=9, height=5, units = "in", dpi = 600)
# 
# 
# # decision: 
# # removing mat_sd seems justified due to better model performance and lower Relative importance in both models
# 
# 
# 
# 
# 
# ##  pet_sd + sea_sd ##########################################################################################
# # correlation with SR and MRD
# cor.test(dat_no.na$sr, dat_no.na$pet_sd, method="s") #0.51
# cor.test(dat_no.na$sr, dat_no.na$sea_sd, method="s") #0.45
# cor.test(dat_no.na$mrd, dat_no.na$pet_sd, method="s") #0.26
# cor.test(dat_no.na$mrd, dat_no.na$sea_sd, method="s") # NS
# 
# cor.test(dat_no.na$pet_sd, dat_no.na$sea_sd, method="s") # 0.69
# 
# 
# # SR TEST"
# gbmFit5_no_pet_sd <- train(sr ~ ., data = dat_no.na[,-grep("pet_sd",names(dat_no.na))], method = "gbm", 
#                            trControl = gbmControl, verbose=FALSE,
#                            tuneGrid=gbmGrid)
# gbmFit5_no_sea_sd <- train(sr ~ ., data = dat_no.na[,-grep("sea_sd",names(dat_no.na))], method = "gbm", 
#                            trControl = gbmControl, verbose=FALSE,
#                            tuneGrid=gbmGrid)
# 
# gbmFit4$results[row.names(gbmFit4$bestTune),]
# gbmFit5_no_pet_sd$results[row.names(gbmFit5_no_pet_sd$bestTune),] # worse
# gbmFit5_no_sea_sd$results[row.names(gbmFit5_no_sea_sd$bestTune),] # little better
# 
# # variable order
# ord <- as.data.frame(summary(gbmFit4))
# res2 <- as.data.frame(summary(gbmFit5_no_sea_sd))
# res3 <- as.data.frame(summary(gbmFit5_no_pet_sd))
# ord_sum <- merge(ord, res2, by="var", all=TRUE)
# ord_sum <- merge(ord_sum, res3, by="var", all=TRUE)
# 
# # plot relative influence in each model instead of order?
# test <- ord_sum[,c(1:4)]
# names(test) <- c("var", "full_model", "no_sea_sd_model", "no_pet_sd_model")
# test2 <- test %>% 
#   pivot_longer(cols=c(full_model, no_sea_sd_model, no_pet_sd_model))
# test2$col <- "other"
# test2$col[test2$var=="pet_sd"] <- "pet_sd"
# test2$col[test2$var=="sea_sd"] <- "sea_sd"
# test2$name <- factor(test2$name, levels = c("no_sea_sd_model", "full_model", "no_pet_sd_model"))
# ggplot(test2, aes(x=name, y=value, group=var, col=test2$col))+
#   geom_point()+
#   geom_line()+
#   scale_color_manual("Variable", values=c("grey", "red", "blue"))+
#   #geom_label(aes(label=var), size=3, col=test2$col)+
#   scale_y_log10("Relative influence")+
#   theme(legend.position = "right")
# 
# 
# ord_sum$no_pet_sd_model <- rank(ord_sum$rel.inf, ties.method = "last", na.last = "keep")
# ord_sum$full_model <- rank(ord_sum$rel.inf.x, ties.method = "last", na.last = "keep")
# ord_sum$no_sea_sd_model <- rank(ord_sum$rel.inf.y, ties.method = "last", na.last = "keep")
# 
# ord_sum2 <- ord_sum %>% 
#   pivot_longer(cols=c(full_model, no_sea_sd_model, no_pet_sd_model))
# #ord_sum2$var2 <- abbreviate(ord_sum2$var, minlength = 7, named = F)
# ord_sum2$name <- factor(ord_sum2$name, levels = c("no_sea_sd_model", "full_model", "no_pet_sd_model"))
# 
# ord_sum2$col <- "grey"
# ord_sum2$col[ord_sum2$var=="pet_sd"] <- "red"
# ord_sum2$col[ord_sum2$var=="sea_sd"] <- "blue"
# ggplot(ord_sum2, aes(x=name, y=value, col=col))+
#   #  geom_point()+
#   geom_line(aes(group=var), col=ord_sum2$col)+
#   geom_label(aes(label=var), size=3, col=ord_sum2$col)+
#   scale_y_continuous("Variable rank")+
#   theme(legend.position = "none")
# #ggsave("results/pet_sd_vs_sea_sd_SR.png", width=6, height=5, units = "in", dpi = 600)
# 
# # MRD TEST"
# gbm_mrd_no_pet_sd <- train(mrd ~ ., data = dat_no.na[,-grep("pet_sd|sr",names(dat_no.na))], method = "gbm", 
#                            trControl = gbmControl, verbose=FALSE,
#                            tuneGrid=gbmGrid)
# gbm_mrd_no_sea_sd <- train(mrd ~ ., data = dat_no.na[,-grep("sea_sd|sr",names(dat_no.na))], method = "gbm", 
#                            trControl = gbmControl, verbose=FALSE,
#                            tuneGrid=gbmGrid)
# 
# 
# gbmFit4_mrd$results[row.names(gbmFit4_mrd$bestTune),]
# gbm_mrd_no_sea_sd$results[row.names(gbm_mrd_no_sea_sd$bestTune),] # better
# gbm_mrd_no_pet_sd$results[row.names(gbm_mrd_no_pet_sd$bestTune),] # better
# 
# 
# # variable order
# ord <- as.data.frame(summary(gbmFit4_mrd))
# res2 <- as.data.frame(summary(gbm_mrd_no_sea_sd))
# res3 <- as.data.frame(summary(gbm_mrd_no_pet_sd))
# ord_sum <- merge(ord, res2, by="var", all=TRUE)
# ord_sum <- merge(ord_sum, res3, by="var", all=TRUE)
# ord_sum$no_pet_sd_model <- rank(ord_sum$rel.inf, ties.method = "last", na.last = "keep")
# ord_sum$full_model <- rank(ord_sum$rel.inf.x, ties.method = "last", na.last = "keep")
# ord_sum$no_sea_sd_model <- rank(ord_sum$rel.inf.y, ties.method = "last", na.last = "keep")
# 
# ord_sum2_mrd <- ord_sum %>% 
#   pivot_longer(cols=c(full_model, no_sea_sd_model, no_pet_sd_model))
# #ord_sum2$var2 <- abbreviate(ord_sum2$var, minlength = 7, named = F)
# ord_sum2_mrd$name <- factor(ord_sum2_mrd$name, levels = c("no_sea_sd_model", "full_model", "no_pet_sd_model"))
# 
# ord_sum2_mrd$col <- "grey"
# ord_sum2_mrd$col[ord_sum2_mrd$var=="pet_sd"] <- "red"
# ord_sum2_mrd$col[ord_sum2_mrd$var=="sea_sd"] <- "blue"
# ggplot(ord_sum2_mrd, aes(x=name, y=value, col=col))+
#   #  geom_point()+
#   geom_line(aes(group=var), col=ord_sum2_mrd$col)+
#   geom_label(aes(label=var), size=3, col=ord_sum2_mrd$col)+
#   scale_y_continuous("Variable rank")+
#   theme(legend.position = "none")
# #ggsave("results/pet_sd_vs_sea_sd_MRD.png", width=6, height=5, units = "in", dpi = 600)
# 
# # FACET PLOT
# ord_sum2$response <- "SR"
# ord_sum2_mrd$response <- "MRD"
# temp <- rbind(ord_sum2, ord_sum2_mrd)
# ggplot(temp, aes(x=name, y=value, group=var, col=temp$col))+ 
#   geom_line()+
#   geom_label(aes(label=var), size=3)+
#   scale_y_continuous("Variable rank")+
#   scale_color_manual(values=c("blue", "grey", "red"))+
#   facet_wrap(~ response)+
#   theme(legend.position = "none", 
#         panel.spacing.x = unit(0, "mm"))
# ggsave("results/pet_sd_vs_sea_sd.png", width=9, height=5, units = "in", dpi = 600)
# 
# 
# # decision: 
# # removing one variable makes tiny diffs. pet_sd has more correlation with SR and MRD, but sea_sd is more important in both models
# 
# 
# 
# 
# 
# 
# ##  soil + number biomes ##########################################################################################
# # correlation with SR and MRD
# cor.test(dat_no.na$sr, dat_no.na$soil, method="s") #0.72
# cor.test(dat_no.na$sr, dat_no.na$number_biomes, method="s") #0.64
# cor.test(dat_no.na$mrd, dat_no.na$soil, method="s") #0.32
# cor.test(dat_no.na$mrd, dat_no.na$number_biomes, method="s") # 0.25
# 
# cor.test(dat_no.na$soil, dat_no.na$number_biomes, method="s") # 0.699
# 
# 
# # SR TEST"
# gbmFit5_no_soil <- train(sr ~ ., data = dat_no.na[,-grep("soil",names(dat_no.na))], method = "gbm", 
#                            trControl = gbmControl, verbose=FALSE,
#                            tuneGrid=gbmGrid)
# gbmFit5_no_number_biomes <- train(sr ~ ., data = dat_no.na[,-grep("number_biomes",names(dat_no.na))], method = "gbm", 
#                            trControl = gbmControl, verbose=FALSE,
#                            tuneGrid=gbmGrid)
# 
# gbmFit4$results[row.names(gbmFit4$bestTune),]
# gbmFit5_no_soil$results[row.names(gbmFit5_no_soil$bestTune),] # worse
# gbmFit5_no_number_biomes$results[row.names(gbmFit5_no_number_biomes$bestTune),] # worse
# 
# # variable order
# ord <- as.data.frame(summary(gbmFit4))
# res2 <- as.data.frame(summary(gbmFit5_no_number_biomes))
# res3 <- as.data.frame(summary(gbmFit5_no_soil))
# ord_sum <- merge(ord, res2, by="var", all=TRUE)
# ord_sum <- merge(ord_sum, res3, by="var", all=TRUE)
# 
# # plot relative influence in each model instead of order?
# test <- ord_sum[,c(1:4)]
# names(test) <- c("var", "full_model", "no_number_biomes_model", "no_soil_model")
# test2 <- test %>% 
#   pivot_longer(cols=c(full_model, no_number_biomes_model, no_soil_model))
# test2$col <- "other"
# test2$col[test2$var=="soil"] <- "soil"
# test2$col[test2$var=="number_biomes"] <- "number_biomes"
# test2$name <- factor(test2$name, levels = c("no_number_biomes_model", "full_model", "no_soil_model"))
# ggplot(test2, aes(x=name, y=value, group=var, col=test2$col))+
#   geom_point()+
#   geom_line()+
#   scale_color_manual("Variable", values=c("red", "grey", "blue"))+
#   #geom_label(aes(label=var), size=3, col=test2$col)+
#   scale_y_log10("Relative influence")+
#   theme(legend.position = "right")
# 
# 
# ord_sum$no_soil_model <- rank(ord_sum$rel.inf, ties.method = "last", na.last = "keep")
# ord_sum$full_model <- rank(ord_sum$rel.inf.x, ties.method = "last", na.last = "keep")
# ord_sum$no_number_biomes_model <- rank(ord_sum$rel.inf.y, ties.method = "last", na.last = "keep")
# 
# ord_sum2 <- ord_sum %>% 
#   pivot_longer(cols=c(full_model, no_number_biomes_model, no_soil_model))
# #ord_sum2$var2 <- abbreviate(ord_sum2$var, minlength = 7, named = F)
# ord_sum2$name <- factor(ord_sum2$name, levels = c("no_number_biomes_model", "full_model", "no_soil_model"))
# 
# ord_sum2$col <- "grey"
# ord_sum2$col[ord_sum2$var=="soil"] <- "red"
# ord_sum2$col[ord_sum2$var=="number_biomes"] <- "blue"
# ggplot(ord_sum2, aes(x=name, y=value, col=col))+
#   #  geom_point()+
#   geom_line(aes(group=var), col=ord_sum2$col)+
#   geom_label(aes(label=var), size=3, col=ord_sum2$col)+
#   scale_y_continuous("Variable rank")+
#   theme(legend.position = "none")
# #ggsave("results/soil_vs_number_biomes_SR.png", width=6, height=5, units = "in", dpi = 600)
# 
# # MRD TEST"
# gbm_mrd_no_soil <- train(mrd ~ ., data = dat_no.na[,-grep("soil|sr",names(dat_no.na))], method = "gbm", 
#                            trControl = gbmControl, verbose=FALSE,
#                            tuneGrid=gbmGrid)
# gbm_mrd_no_number_biomes <- train(mrd ~ ., data = dat_no.na[,-grep("number_biomes|sr",names(dat_no.na))], method = "gbm", 
#                            trControl = gbmControl, verbose=FALSE,
#                            tuneGrid=gbmGrid)
# 
# 
# gbmFit4_mrd$results[row.names(gbmFit4_mrd$bestTune),]
# gbm_mrd_no_number_biomes$results[row.names(gbm_mrd_no_number_biomes$bestTune),] # better
# gbm_mrd_no_soil$results[row.names(gbm_mrd_no_soil$bestTune),] # better
# 
# 
# # variable order
# ord <- as.data.frame(summary(gbmFit4_mrd))
# res2 <- as.data.frame(summary(gbm_mrd_no_number_biomes))
# res3 <- as.data.frame(summary(gbm_mrd_no_soil))
# ord_sum <- merge(ord, res2, by="var", all=TRUE)
# ord_sum <- merge(ord_sum, res3, by="var", all=TRUE)
# ord_sum$no_soil_model <- rank(ord_sum$rel.inf, ties.method = "last", na.last = "keep")
# ord_sum$full_model <- rank(ord_sum$rel.inf.x, ties.method = "last", na.last = "keep")
# ord_sum$no_number_biomes_model <- rank(ord_sum$rel.inf.y, ties.method = "last", na.last = "keep")
# 
# ord_sum2_mrd <- ord_sum %>% 
#   pivot_longer(cols=c(full_model, no_number_biomes_model, no_soil_model))
# #ord_sum2$var2 <- abbreviate(ord_sum2$var, minlength = 7, named = F)
# ord_sum2_mrd$name <- factor(ord_sum2_mrd$name, levels = c("no_number_biomes_model", "full_model", "no_soil_model"))
# 
# ord_sum2_mrd$col <- "grey"
# ord_sum2_mrd$col[ord_sum2_mrd$var=="soil"] <- "red"
# ord_sum2_mrd$col[ord_sum2_mrd$var=="number_biomes"] <- "blue"
# ggplot(ord_sum2_mrd, aes(x=name, y=value, col=col))+
#   #  geom_point()+
#   geom_line(aes(group=var), col=ord_sum2_mrd$col)+
#   geom_label(aes(label=var), size=3, col=ord_sum2_mrd$col)+
#   scale_y_continuous("Variable rank")+
#   theme(legend.position = "none")
# #ggsave("results/soil_vs_number_biomes_MRD.png", width=6, height=5, units = "in", dpi = 600)
# 
# # FACET PLOT
# ord_sum2$response <- "SR"
# ord_sum2_mrd$response <- "MRD"
# temp <- rbind(ord_sum2, ord_sum2_mrd)
# ggplot(temp, aes(x=name, y=value, group=var, col=col))+ 
#   geom_line()+
#   geom_label(aes(label=var), size=3)+
#   scale_y_continuous("Variable rank")+
#   scale_color_manual(values=c("blue", "grey", "red"))+
#   facet_wrap(~ response)+
#   theme(legend.position = "none", 
#         panel.spacing.x = unit(0, "mm"))
# ggsave("results/soil_vs_number_biomes.png", width=9, height=5, units = "in", dpi = 600)
# 
# # I would argue to exlclude number of biomes since its a rougher resilution than number of soil types, and it seems to code a similar thing, but not in as much detail as soil since the habitat types become more important when removing soil, indicating soil representing physical characteristics of those biomes
# 
# 
# 
# ## pre_mean + subtrop mbf ################################################################################33
# cor.test(dat_no.na$pre_mean, dat_no.na$`(sub)tropical moist broadleaf forest`, method="s")
# gbmFit5_no_pre_mean <- train(sr ~ ., data = dat_no.na[,-grep("pre_mean",names(dat_no.na))], method = "gbm", 
#                              trControl = gbmControl, verbose=FALSE,
#                              tuneGrid=gbmGrid)
# gbmFit5_no_trf <- train(sr ~ ., data = dat_no.na[,-grep("moist",names(dat_no.na))], method = "gbm", 
#                         trControl = gbmControl, verbose=FALSE,
#                         tuneGrid=gbmGrid)
# 
# gbmFit5_no_pre_mean$results[row.names(gbmFit5_no_pre_mean$bestTune),]
# gbmFit5_no_trf$results[row.names(gbmFit5_no_trf$bestTune),]
# #plot(varImp(gbmFit5_no_soil, scale = TRUE), main="gbmFit5_no_pet_sd")
# #plot(varImp(gbmFit5_no_nb, scale = TRUE), main="gbmFit5_no_sea_sd")
# # models are almost the same
# 
# cor.test(dat_no.na$sr, dat_no.na$pre_mean, method="s") #0.21
# cor.test(dat_no.na$sr, dat_no.na$`(sub)tropical moist broadleaf forest`, method="s") #0.47
# cor.test(dat_no.na$mrd, dat_no.na$pre_mean, method="s") # -0.56
# cor.test(dat_no.na$mrd, dat_no.na$`(sub)tropical moist broadleaf forest`, method="s") # -0.39
# 
# # variable order
# ord <- as.data.frame(summary(gbmFit4))
# res2 <- as.data.frame(summary(gbmFit5_no_pre_mean))
# res3 <- as.data.frame(summary(gbmFit5_no_trf))
# ord_sum <- merge(ord, res2, by="var", all=TRUE)
# ord_sum <- merge(ord_sum, res3, by="var", all=TRUE)
# ord_sum$full_model <- rank(ord_sum$rel.inf.x, ties.method = "last", na.last = "keep")
# ord_sum$no_pre_mean_model <- rank(ord_sum$rel.inf.y, ties.method = "last", na.last = "keep")
# ord_sum$no_trf_model <- rank(ord_sum$rel.inf, ties.method = "last", na.last = "keep")
# ord_sum2 <- ord_sum %>% 
#   pivot_longer(cols=c(full_model, no_pre_mean_model, no_trf_model))
# ord_sum2$var2 <- abbreviate(ord_sum2$var, minlength = 7, named = F)
# ord_sum2$name <- factor(ord_sum2$name, levels = c("no_pre_mean_model", "full_model", "no_trf_model"))
# ggplot(ord_sum2, aes(x=name, y=value, col=var2))+
#   geom_point()+
#   geom_line(aes(group=var2))+
#   geom_label(aes(label=var2), size=3)+
#   scale_y_continuous("importance")+
#   theme(legend.position = "none")
# 
# par(mfrow=c(2,2))
# plot(dat_no.na$pre_mean, dat_no.na$`(sub)tropical moist broadleaf forest`)
# plot(dat_no.na$sr, dat_no.na$pre_mean)
# plot(dat_no.na$sr, dat_no.na$`(sub)tropical moist broadleaf forest`)
# 
# # I would argue since TRF is more important for SR and pre_mean for MRD, and they do not change their positions much in the model, we keep both
# 
# 
# 
# ##  tra_sd + area ##########################################################################################
# # correlation with SR and MRD
# cor.test(dat_no.na$sr, dat_no.na$tra_sd, method="s") #0.18
# cor.test(dat_no.na$sr, dat_no.na$area, method="s") #0.46
# cor.test(dat_no.na$mrd, dat_no.na$tra_sd, method="s") #0.24
# cor.test(dat_no.na$mrd, dat_no.na$area, method="s") # 0.21
# 
# cor.test(dat_no.na$tra_sd, dat_no.na$area, method="s") # 0.73
# 
# 
# # SR TEST"
# gbmFit5_no_tra_sd <- train(sr ~ ., data = dat_no.na[,-grep("tra_sd",names(dat_no.na))], method = "gbm", 
#                          trControl = gbmControl, verbose=FALSE,
#                          tuneGrid=gbmGrid)
# gbmFit5_no_area <- train(sr ~ ., data = dat_no.na[,-grep("area",names(dat_no.na))], method = "gbm", 
#                                   trControl = gbmControl, verbose=FALSE,
#                                   tuneGrid=gbmGrid)
# 
# gbmFit4$results[row.names(gbmFit4$bestTune),]
# gbmFit5_no_tra_sd$results[row.names(gbmFit5_no_tra_sd$bestTune),] # worse
# gbmFit5_no_area$results[row.names(gbmFit5_no_area$bestTune),] # worse
# 
# # variable order
# ord <- as.data.frame(summary(gbmFit4))
# res2 <- as.data.frame(summary(gbmFit5_no_area))
# res3 <- as.data.frame(summary(gbmFit5_no_tra_sd))
# ord_sum <- merge(ord, res2, by="var", all=TRUE)
# ord_sum <- merge(ord_sum, res3, by="var", all=TRUE)
# 
# # plot relative influence in each model instead of order?
# test <- ord_sum[,c(1:4)]
# names(test) <- c("var", "full_model", "no_area_model", "no_tra_sd_model")
# test2 <- test %>% 
#   pivot_longer(cols=c(full_model, no_area_model, no_tra_sd_model))
# test2$col <- "other"
# test2$col[test2$var=="tra_sd"] <- "tra_sd"
# test2$col[test2$var=="area"] <- "area"
# test2$name <- factor(test2$name, levels = c("no_area_model", "full_model", "no_tra_sd_model"))
# ggplot(test2, aes(x=name, y=value, group=var, col=test2$col))+
#   geom_point()+
#   geom_line()+
#   scale_color_manual("Variable", values=c("red", "grey", "blue"))+
#   #geom_label(aes(label=var), size=3, col=test2$col)+
#   scale_y_log10("Relative influence")+
#   theme(legend.position = "right")
# 
# 
# ord_sum$no_tra_sd_model <- rank(ord_sum$rel.inf, ties.method = "last", na.last = "keep")
# ord_sum$full_model <- rank(ord_sum$rel.inf.x, ties.method = "last", na.last = "keep")
# ord_sum$no_area_model <- rank(ord_sum$rel.inf.y, ties.method = "last", na.last = "keep")
# 
# ord_sum2 <- ord_sum %>% 
#   pivot_longer(cols=c(full_model, no_area_model, no_tra_sd_model))
# #ord_sum2$var2 <- abbreviate(ord_sum2$var, minlength = 7, named = F)
# ord_sum2$name <- factor(ord_sum2$name, levels = c("no_area_model", "full_model", "no_tra_sd_model"))
# 
# ord_sum2$col <- "grey"
# ord_sum2$col[ord_sum2$var=="tra_sd"] <- "red"
# ord_sum2$col[ord_sum2$var=="area"] <- "blue"
# ggplot(ord_sum2, aes(x=name, y=value, col=col))+
#   #  geom_point()+
#   geom_line(aes(group=var), col=ord_sum2$col)+
#   geom_label(aes(label=var), size=3, col=ord_sum2$col)+
#   scale_y_continuous("Variable rank")+
#   theme(legend.position = "none")
# #ggsave("results/tra_sd_vs_area_SR.png", width=6, height=5, units = "in", dpi = 600)
# 
# # MRD TEST"
# gbm_mrd_no_tra_sd <- train(mrd ~ ., data = dat_no.na[,-grep("tra_sd|sr",names(dat_no.na))], method = "gbm", 
#                          trControl = gbmControl, verbose=FALSE,
#                          tuneGrid=gbmGrid)
# gbm_mrd_no_area <- train(mrd ~ ., data = dat_no.na[,-grep("area|sr",names(dat_no.na))], method = "gbm", 
#                                   trControl = gbmControl, verbose=FALSE,
#                                   tuneGrid=gbmGrid)
# 
# 
# gbmFit4_mrd$results[row.names(gbmFit4_mrd$bestTune),]
# gbm_mrd_no_area$results[row.names(gbm_mrd_no_area$bestTune),] # better
# gbm_mrd_no_tra_sd$results[row.names(gbm_mrd_no_tra_sd$bestTune),] # better
# 
# 
# # variable order
# ord <- as.data.frame(summary(gbmFit4_mrd))
# res2 <- as.data.frame(summary(gbm_mrd_no_area))
# res3 <- as.data.frame(summary(gbm_mrd_no_tra_sd))
# ord_sum <- merge(ord, res2, by="var", all=TRUE)
# ord_sum <- merge(ord_sum, res3, by="var", all=TRUE)
# ord_sum$no_tra_sd_model <- rank(ord_sum$rel.inf, ties.method = "last", na.last = "keep")
# ord_sum$full_model <- rank(ord_sum$rel.inf.x, ties.method = "last", na.last = "keep")
# ord_sum$no_area_model <- rank(ord_sum$rel.inf.y, ties.method = "last", na.last = "keep")
# 
# ord_sum2_mrd <- ord_sum %>% 
#   pivot_longer(cols=c(full_model, no_area_model, no_tra_sd_model))
# #ord_sum2$var2 <- abbreviate(ord_sum2$var, minlength = 7, named = F)
# ord_sum2_mrd$name <- factor(ord_sum2_mrd$name, levels = c("no_area_model", "full_model", "no_tra_sd_model"))
# 
# ord_sum2_mrd$col <- "grey"
# ord_sum2_mrd$col[ord_sum2_mrd$var=="tra_sd"] <- "red"
# ord_sum2_mrd$col[ord_sum2_mrd$var=="area"] <- "blue"
# ggplot(ord_sum2_mrd, aes(x=name, y=value, col=col))+
#   #  geom_point()+
#   geom_line(aes(group=var), col=ord_sum2_mrd$col)+
#   geom_label(aes(label=var), size=3, col=ord_sum2_mrd$col)+
#   scale_y_continuous("Variable rank")+
#   theme(legend.position = "none")
# #ggsave("results/tra_sd_vs_area_MRD.png", width=6, height=5, units = "in", dpi = 600)
# 
# # FACET PLOT
# ord_sum2$response <- "SR"
# ord_sum2_mrd$response <- "MRD"
# temp <- rbind(ord_sum2, ord_sum2_mrd)
# ggplot(temp, aes(x=name, y=value, group=var, col=col))+ 
#   geom_line()+
#   geom_label(aes(label=var), size=3)+
#   scale_y_continuous("Variable rank")+
#   scale_color_manual(values=c("blue", "grey", "red"))+
#   facet_wrap(~ response)+
#   theme(legend.position = "none", 
#         panel.spacing.x = unit(0, "mm"))
# ggsave("results/tra_sd_vs_area.png", width=9, height=5, units = "in", dpi = 600)
# 
# # worse models for SR, better models for MRD. would not exclude any of them since correlation with area is expected and the strength is not super high. correlations with SR and MRD also work in the same direction
# 
# 
# ## pre_mean + subtrop mbf ################################################################################33
# cor.test(dat_no.na$pre_mean, dat_no.na$`(sub)tropical moist broadleaf forest`, method="s")
# gbmFit5_no_pre_mean <- train(sr ~ ., data = dat_no.na[,-grep("pre_mean",names(dat_no.na))], method = "gbm", 
#                              trControl = gbmControl, verbose=FALSE,
#                              tuneGrid=gbmGrid)
# gbmFit5_no_trf <- train(sr ~ ., data = dat_no.na[,-grep("moist",names(dat_no.na))], method = "gbm", 
#                         trControl = gbmControl, verbose=FALSE,
#                         tuneGrid=gbmGrid)
# 
# gbmFit5_no_pre_mean$results[row.names(gbmFit5_no_pre_mean$bestTune),]
# gbmFit5_no_trf$results[row.names(gbmFit5_no_trf$bestTune),]
# #plot(varImp(gbmFit5_no_soil, scale = TRUE), main="gbmFit5_no_pet_sd")
# #plot(varImp(gbmFit5_no_nb, scale = TRUE), main="gbmFit5_no_sea_sd")
# # models are almost the same
# 
# cor.test(dat_no.na$sr, dat_no.na$pre_mean, method="s") #0.21
# cor.test(dat_no.na$sr, dat_no.na$`(sub)tropical moist broadleaf forest`, method="s") #0.47
# cor.test(dat_no.na$mrd, dat_no.na$pre_mean, method="s") # -0.56
# cor.test(dat_no.na$mrd, dat_no.na$`(sub)tropical moist broadleaf forest`, method="s") # -0.39
# 
# # variable order
# ord <- as.data.frame(summary(gbmFit4))
# res2 <- as.data.frame(summary(gbmFit5_no_pre_mean))
# res3 <- as.data.frame(summary(gbmFit5_no_trf))
# ord_sum <- merge(ord, res2, by="var", all=TRUE)
# ord_sum <- merge(ord_sum, res3, by="var", all=TRUE)
# ord_sum$full_model <- rank(ord_sum$rel.inf.x, ties.method = "last", na.last = "keep")
# ord_sum$no_pre_mean_model <- rank(ord_sum$rel.inf.y, ties.method = "last", na.last = "keep")
# ord_sum$no_trf_model <- rank(ord_sum$rel.inf, ties.method = "last", na.last = "keep")
# ord_sum2 <- ord_sum %>% 
#   pivot_longer(cols=c(full_model, no_pre_mean_model, no_trf_model))
# ord_sum2$var2 <- abbreviate(ord_sum2$var, minlength = 7, named = F)
# ord_sum2$name <- factor(ord_sum2$name, levels = c("no_pre_mean_model", "full_model", "no_trf_model"))
# ggplot(ord_sum2, aes(x=name, y=value, col=var2))+
#   geom_point()+
#   geom_line(aes(group=var2))+
#   geom_label(aes(label=var2), size=3)+
#   scale_y_continuous("importance")+
#   theme(legend.position = "none")
# 
# par(mfrow=c(2,2))
# plot(dat_no.na$pre_mean, dat_no.na$`(sub)tropical moist broadleaf forest`)
# plot(dat_no.na$sr, dat_no.na$pre_mean)
# plot(dat_no.na$sr, dat_no.na$`(sub)tropical moist broadleaf forest`)
# 
# # I would argue since TRF is more important for SR and pre_mean for MRD, and they do not change their positions much in the model, we keep both
# 
# 
# 
# 
# ##  pre_mean + subtrop mbf ##########################################################################################
# # correlation with SR and MRD
# cor.test(dat_no.na$sr, dat_no.na$pre_mean, method="s") #0.20
# cor.test(dat_no.na$sr, dat_no.na$sub_trop_mbf, method="s") #0.48
# cor.test(dat_no.na$mrd, dat_no.na$pre_mean, method="s") # - 0.56
# cor.test(dat_no.na$mrd, dat_no.na$sub_trop_mbf, method="s") # -0.39
# 
# cor.test(dat_no.na$pre_mean, dat_no.na$sub_trop_mbf, method="s") # 0.69
# 
# 
# # SR TEST"
# gbmFit5_no_pre_mean <- train(sr ~ ., data = dat_no.na[,-grep("pre_mean",names(dat_no.na))], method = "gbm", 
#                            trControl = gbmControl, verbose=FALSE,
#                            tuneGrid=gbmGrid)
# gbmFit5_no_trop_mbf <- train(sr ~ ., data = dat_no.na[,-grep("mbf",names(dat_no.na))], method = "gbm", 
#                          trControl = gbmControl, verbose=FALSE,
#                          tuneGrid=gbmGrid)
# 
# gbmFit4$results[row.names(gbmFit4$bestTune),]
# gbmFit5_no_pre_mean$results[row.names(gbmFit5_no_pre_mean$bestTune),] # worse
# gbmFit5_no_trop_mbf$results[row.names(gbmFit5_no_trop_mbf$bestTune),] # worse
# 
# # variable order
# ord <- as.data.frame(summary(gbmFit4))
# res2 <- as.data.frame(summary(gbmFit5_no_trop_mbf))
# res3 <- as.data.frame(summary(gbmFit5_no_pre_mean))
# ord_sum <- merge(ord, res2, by="var", all=TRUE)
# ord_sum <- merge(ord_sum, res3, by="var", all=TRUE)
# 
# # plot relative influence in each model instead of order?
# test <- ord_sum[,c(1:4)]
# names(test) <- c("var", "full_model", "no_trop_mbf_model", "no_pre_mean_model")
# test2 <- test %>% 
#   pivot_longer(cols=c(full_model, no_trop_mbf_model, no_pre_mean_model))
# test2$col <- "other"
# test2$col[test2$var=="pre_mean"] <- "pre_mean"
# test2$col[test2$var=="sub_trop_mbf"] <- "trop_mbf"
# test2$name <- factor(test2$name, levels = c("no_trop_mbf_model", "full_model", "no_pre_mean_model"))
# ggplot(test2, aes(x=name, y=value, group=var, col=col))+
#   geom_point()+
#   geom_line()+
#   scale_color_manual("Variable", values=c("grey", "red", "blue"))+
#   #geom_label(aes(label=var), size=3, col=test2$col)+
#   scale_y_log10("Relative influence")+
#   theme(legend.position = "right")
# 
# 
# ord_sum$no_pre_mean_model <- rank(ord_sum$rel.inf, ties.method = "last", na.last = "keep")
# ord_sum$full_model <- rank(ord_sum$rel.inf.x, ties.method = "last", na.last = "keep")
# ord_sum$no_trop_mbf_model <- rank(ord_sum$rel.inf.y, ties.method = "last", na.last = "keep")
# 
# ord_sum2 <- ord_sum %>% 
#   pivot_longer(cols=c(full_model, no_trop_mbf_model, no_pre_mean_model))
# #ord_sum2$var2 <- abbreviate(ord_sum2$var, minlength = 7, named = F)
# ord_sum2$name <- factor(ord_sum2$name, levels = c("no_trop_mbf_model", "full_model", "no_pre_mean_model"))
# 
# ord_sum2$col <- "grey"
# ord_sum2$col[ord_sum2$var=="pre_mean"] <- "red"
# ord_sum2$col[ord_sum2$var=="sub_trop_mbf"] <- "blue"
# ggplot(ord_sum2, aes(x=name, y=value, col=col))+
#   #  geom_point()+
#   geom_line(aes(group=var), col=ord_sum2$col)+
#   geom_label(aes(label=var), size=3, col=ord_sum2$col)+
#   scale_y_continuous("Variable rank")+
#   theme(legend.position = "none")
# #ggsave("results/pre_mean_vs_trop_mbf_SR.png", width=6, height=5, units = "in", dpi = 600)
# 
# # MRD TEST"
# gbm_mrd_no_pre_mean <- train(mrd ~ ., data = dat_no.na[,-grep("pre_mean|sr",names(dat_no.na))], method = "gbm", 
#                            trControl = gbmControl, verbose=FALSE,
#                            tuneGrid=gbmGrid)
# gbm_mrd_no_trop_mbf <- train(mrd ~ ., data = dat_no.na[,-grep("mbf|sr",names(dat_no.na))], method = "gbm", 
#                          trControl = gbmControl, verbose=FALSE,
#                          tuneGrid=gbmGrid)
# 
# 
# gbmFit4_mrd$results[row.names(gbmFit4_mrd$bestTune),]
# gbm_mrd_no_trop_mbf$results[row.names(gbm_mrd_no_trop_mbf$bestTune),] # better
# gbm_mrd_no_pre_mean$results[row.names(gbm_mrd_no_pre_mean$bestTune),] # worse
# 
# 
# # variable order
# ord <- as.data.frame(summary(gbmFit4_mrd))
# res2 <- as.data.frame(summary(gbm_mrd_no_trop_mbf))
# res3 <- as.data.frame(summary(gbm_mrd_no_pre_mean))
# ord_sum <- merge(ord, res2, by="var", all=TRUE)
# ord_sum <- merge(ord_sum, res3, by="var", all=TRUE)
# ord_sum$no_pre_mean_model <- rank(ord_sum$rel.inf, ties.method = "last", na.last = "keep")
# ord_sum$full_model <- rank(ord_sum$rel.inf.x, ties.method = "last", na.last = "keep")
# ord_sum$no_trop_mbf_model <- rank(ord_sum$rel.inf.y, ties.method = "last", na.last = "keep")
# 
# ord_sum2_mrd <- ord_sum %>% 
#   pivot_longer(cols=c(full_model, no_trop_mbf_model, no_pre_mean_model))
# #ord_sum2$var2 <- abbreviate(ord_sum2$var, minlength = 7, named = F)
# ord_sum2_mrd$name <- factor(ord_sum2_mrd$name, levels = c("no_trop_mbf_model", "full_model", "no_pre_mean_model"))
# 
# ord_sum2_mrd$col <- "grey"
# ord_sum2_mrd$col[ord_sum2_mrd$var=="pre_mean"] <- "red"
# ord_sum2_mrd$col[ord_sum2_mrd$var=="sub_trop_mbf"] <- "blue"
# ggplot(ord_sum2_mrd, aes(x=name, y=value, col=col))+
#   #  geom_point()+
#   geom_line(aes(group=var), col=ord_sum2_mrd$col)+
#   geom_label(aes(label=var), size=3, col=ord_sum2_mrd$col)+
#   scale_y_continuous("Variable rank")+
#   theme(legend.position = "none")
# #ggsave("results/pre_mean_vs_trop_mbf_MRD.png", width=6, height=5, units = "in", dpi = 600)
# 
# # FACET PLOT
# ord_sum2$response <- "SR"
# ord_sum2_mrd$response <- "MRD"
# temp <- rbind(ord_sum2, ord_sum2_mrd)
# ggplot(temp, aes(x=name, y=value, group=var, col=col))+ 
#   geom_line()+
#   geom_label(aes(label=var), size=3)+
#   scale_y_continuous("Variable rank")+
#   scale_color_manual(values=c("blue", "grey", "red"))+
#   facet_wrap(~ response)+
#   theme(legend.position = "none", 
#         panel.spacing.x = unit(0, "mm"))
# ggsave("results/pre_mean_vs_trop_mbf.png", width=9, height=5, units = "in", dpi = 600)
# 
# # positive correlation with SR, negative (and rather strong) with MRD, with differing importances in the model. 
# # Removing subtrop mbf in the SR model causes pre_mean to become much more important, suggesting they account for similar characteristics. However, subtrop mbf does not change in its influence when removing pre_mean.
# # REmoving subtrop mbf in the MRD model does not change pre_mean position, removing pre_mean makes subtrop mbf more important. 
# # --> opposing effects and importances: keep both, account for correlation, however they cannot be part of one latent variable
# 
# 
# 
# 
# 
# ##  pet_mean + mat_mean ##########################################################################################
# # correlation with SR and MRD
# cor.test(dat_no.na$sr, dat_no.na$mat_mean, method="s") #0.28
# cor.test(dat_no.na$sr, dat_no.na$pet_mean, method="s") #0.31
# cor.test(dat_no.na$mrd, dat_no.na$mat_mean, method="s") # -0.22
# cor.test(dat_no.na$mrd, dat_no.na$pet_mean, method="s") # NS
# 
# cor.test(dat_no.na$mat_mean, dat_no.na$pet_mean, method="s") # 0.80
# 
# 
# # SR TEST"
# gbmFit5_no_mat_mean <- train(sr ~ ., data = dat_no.na[,-grep("mat_mean",names(dat_no.na))], method = "gbm", 
#                            trControl = gbmControl, verbose=FALSE,
#                            tuneGrid=gbmGrid)
# gbmFit5_no_pet_mean <- train(sr ~ ., data = dat_no.na[,-grep("pet_mean",names(dat_no.na))], method = "gbm", 
#                            trControl = gbmControl, verbose=FALSE,
#                            tuneGrid=gbmGrid)
# 
# gbmFit4$results[row.names(gbmFit4$bestTune),]
# gbmFit5_no_mat_mean$results[row.names(gbmFit5_no_mat_mean$bestTune),] # worse
# gbmFit5_no_pet_mean$results[row.names(gbmFit5_no_pet_mean$bestTune),] # worse
# 
# # variable order
# ord <- as.data.frame(summary(gbmFit4))
# res2 <- as.data.frame(summary(gbmFit5_no_pet_mean))
# res3 <- as.data.frame(summary(gbmFit5_no_mat_mean))
# ord_sum <- merge(ord, res2, by="var", all=TRUE)
# ord_sum <- merge(ord_sum, res3, by="var", all=TRUE)
# 
# # plot relative influence in each model instead of order?
# test <- ord_sum[,c(1:4)]
# names(test) <- c("var", "full_model", "no_pet_mean_model", "no_mat_mean_model")
# test2 <- test %>% 
#   pivot_longer(cols=c(full_model, no_pet_mean_model, no_mat_mean_model))
# test2$col <- "other"
# test2$col[test2$var=="mat_mean"] <- "mat_mean"
# test2$col[test2$var=="pet_mean"] <- "pet_mean"
# test2$name <- factor(test2$name, levels = c("no_pet_mean_model", "full_model", "no_mat_mean_model"))
# ggplot(test2, aes(x=name, y=value, group=var, col=test2$col))+
#   geom_point()+
#   geom_line()+
#   scale_color_manual("Variable", values=c("red", "grey", "blue"))+
#   #geom_label(aes(label=var), size=3, col=test2$col)+
#   scale_y_log10("Relative influence")+
#   theme(legend.position = "right")
# # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# ord_sum$no_mat_mean_model <- rank(ord_sum$rel.inf, ties.method = "last", na.last = "keep")
# ord_sum$full_model <- rank(ord_sum$rel.inf.x, ties.method = "last", na.last = "keep")
# ord_sum$no_pet_mean_model <- rank(ord_sum$rel.inf.y, ties.method = "last", na.last = "keep")
# 
# ord_sum2 <- ord_sum %>% 
#   pivot_longer(cols=c(full_model, no_pet_mean_model, no_mat_mean_model))
# #ord_sum2$var2 <- abbreviate(ord_sum2$var, minlength = 7, named = F)
# ord_sum2$name <- factor(ord_sum2$name, levels = c("no_pet_mean_model", "full_model", "no_mat_mean_model"))
# 
# ord_sum2$col <- "grey"
# ord_sum2$col[ord_sum2$var=="mat_mean"] <- "red"
# ord_sum2$col[ord_sum2$var=="pet_mean"] <- "blue"
# ggplot(ord_sum2, aes(x=name, y=value, col=col))+
#   #  geom_point()+
#   geom_line(aes(group=var), col=ord_sum2$col)+
#   geom_label(aes(label=var), size=3, col=ord_sum2$col)+
#   scale_y_continuous("Variable rank")+
#   theme(legend.position = "none")
# #ggsave("results/mat_mean_vs_pet_mean_SR.png", width=6, height=5, units = "in", dpi = 600)
# 
# # MRD TEST"
# gbm_mrd_no_mat_mean <- train(mrd ~ ., data = dat_no.na[,-grep("mat_mean|sr",names(dat_no.na))], method = "gbm", 
#                            trControl = gbmControl, verbose=FALSE,
#                            tuneGrid=gbmGrid)
# gbm_mrd_no_pet_mean <- train(mrd ~ ., data = dat_no.na[,-grep("pet_mean|sr",names(dat_no.na))], method = "gbm", 
#                            trControl = gbmControl, verbose=FALSE,
#                            tuneGrid=gbmGrid)
# 
# 
# gbmFit4_mrd$results[row.names(gbmFit4_mrd$bestTune),]
# gbm_mrd_no_pet_mean$results[row.names(gbm_mrd_no_pet_mean$bestTune),] # better
# gbm_mrd_no_mat_mean$results[row.names(gbm_mrd_no_mat_mean$bestTune),] # worse
# 
# 
# # variable order
# ord <- as.data.frame(summary(gbmFit4_mrd))
# res2 <- as.data.frame(summary(gbm_mrd_no_pet_mean))
# res3 <- as.data.frame(summary(gbm_mrd_no_mat_mean))
# ord_sum <- merge(ord, res2, by="var", all=TRUE)
# ord_sum <- merge(ord_sum, res3, by="var", all=TRUE)
# ord_sum$no_mat_mean_model <- rank(ord_sum$rel.inf, ties.method = "last", na.last = "keep")
# ord_sum$full_model <- rank(ord_sum$rel.inf.x, ties.method = "last", na.last = "keep")
# ord_sum$no_pet_mean_model <- rank(ord_sum$rel.inf.y, ties.method = "last", na.last = "keep")
# 
# ord_sum2_mrd <- ord_sum %>% 
#   pivot_longer(cols=c(full_model, no_pet_mean_model, no_mat_mean_model))
# #ord_sum2$var2 <- abbreviate(ord_sum2$var, minlength = 7, named = F)
# ord_sum2_mrd$name <- factor(ord_sum2_mrd$name, levels = c("no_pet_mean_model", "full_model", "no_mat_mean_model"))
# 
# ord_sum2_mrd$col <- "grey"
# ord_sum2_mrd$col[ord_sum2_mrd$var=="mat_mean"] <- "red"
# ord_sum2_mrd$col[ord_sum2_mrd$var=="pet_mean"] <- "blue"
# ggplot(ord_sum2_mrd, aes(x=name, y=value, col=col))+
#   #  geom_point()+
#   geom_line(aes(group=var), col=ord_sum2_mrd$col)+
#   geom_label(aes(label=var), size=3, col=ord_sum2_mrd$col)+
#   scale_y_continuous("Variable rank")+
#   theme(legend.position = "none")
# #ggsave("results/mat_mean_vs_pet_mean_MRD.png", width=6, height=5, units = "in", dpi = 600)
# 
# # FACET PLOT
# ord_sum2$response <- "SR"
# ord_sum2_mrd$response <- "MRD"
# temp <- rbind(ord_sum2, ord_sum2_mrd)
# ggplot(temp, aes(x=name, y=value, group=var, col=col))+ 
#   geom_line()+
#   geom_label(aes(label=var), size=3)+
#   scale_y_continuous("Variable rank")+
#   scale_color_manual(values=c("blue", "grey", "red"))+
#   facet_wrap(~ response)+
#   theme(legend.position = "none", 
#         panel.spacing.x = unit(0, "mm"))
# ggsave("results/mat_mean_vs_pet_mean.png", width=9, height=5, units = "in", dpi = 600)
# 
# 
# # decision: 
# # removing one variable makes tiny diffs. mat_mean has more correlation with SR and importance in SR model, else no diffs.
# 
# 
# ##  tra_mean + mat_mean ##########################################################################################
# # correlation with SR and MRD
# cor.test(dat_no.na$sr, dat_no.na$mat_mean, method="s") # 0.28
# cor.test(dat_no.na$sr, dat_no.na$tra_mean, method="s") # -0.28
# cor.test(dat_no.na$mrd, dat_no.na$mat_mean, method="s") # -0.22
# cor.test(dat_no.na$mrd, dat_no.na$tra_mean, method="s") # 0.30
# 
# cor.test(dat_no.na$mat_mean, dat_no.na$tra_mean, method="s") # - 0.81
# 
# 
# # SR TEST"
# gbmFit5_no_tra_mean <- train(sr ~ ., data = dat_no.na[,-grep("tra_mean",names(dat_no.na))], method = "gbm", 
#                              trControl = gbmControl, verbose=FALSE,
#                              tuneGrid=gbmGrid)
# 
# gbmFit4$results[row.names(gbmFit4$bestTune),]
# gbmFit5_no_mat_mean$results[row.names(gbmFit5_no_mat_mean$bestTune),] # worse
# gbmFit5_no_tra_mean$results[row.names(gbmFit5_no_tra_mean$bestTune),] # better
# 
# # variable order
# ord <- as.data.frame(summary(gbmFit4))
# res2 <- as.data.frame(summary(gbmFit5_no_tra_mean))
# res3 <- as.data.frame(summary(gbmFit5_no_mat_mean))
# ord_sum <- merge(ord, res2, by="var", all=TRUE)
# ord_sum <- merge(ord_sum, res3, by="var", all=TRUE)
# 
# # plot relative influence in each model instead of order?
# test <- ord_sum[,c(1:4)]
# names(test) <- c("var", "full_model", "no_tra_mean_model", "no_mat_mean_model")
# test2 <- test %>% 
#   pivot_longer(cols=c(full_model, no_tra_mean_model, no_mat_mean_model))
# test2$col <- "other"
# test2$col[test2$var=="mat_mean"] <- "mat_mean"
# test2$col[test2$var=="tra_mean"] <- "tra_mean"
# test2$name <- factor(test2$name, levels = c("no_tra_mean_model", "full_model", "no_mat_mean_model"))
# ggplot(test2, aes(x=name, y=value, group=var, col=test2$col))+
#   geom_point()+
#   geom_line()+
#   scale_color_manual("Variable", values=c("red", "grey", "blue"))+
#   #geom_label(aes(label=var), size=3, col=test2$col)+
#   scale_y_log10("Relative influence")+
#   theme(legend.position = "right")
# # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# ord_sum$no_mat_mean_model <- rank(ord_sum$rel.inf, ties.method = "last", na.last = "keep")
# ord_sum$full_model <- rank(ord_sum$rel.inf.x, ties.method = "last", na.last = "keep")
# ord_sum$no_tra_mean_model <- rank(ord_sum$rel.inf.y, ties.method = "last", na.last = "keep")
# 
# ord_sum2 <- ord_sum %>% 
#   pivot_longer(cols=c(full_model, no_tra_mean_model, no_mat_mean_model))
# #ord_sum2$var2 <- abbreviate(ord_sum2$var, minlength = 7, named = F)
# ord_sum2$name <- factor(ord_sum2$name, levels = c("no_tra_mean_model", "full_model", "no_mat_mean_model"))
# 
# ord_sum2$col <- "grey"
# ord_sum2$col[ord_sum2$var=="mat_mean"] <- "red"
# ord_sum2$col[ord_sum2$var=="tra_mean"] <- "blue"
# ggplot(ord_sum2, aes(x=name, y=value, col=col))+
#   #  geom_point()+
#   geom_line(aes(group=var), col=ord_sum2$col)+
#   geom_label(aes(label=var), size=3, col=ord_sum2$col)+
#   scale_y_continuous("Variable rank")+
#   theme(legend.position = "none")
# #ggsave("results/mat_mean_vs_tra_mean_SR.png", width=6, height=5, units = "in", dpi = 600)
# 
# # MRD TEST"
# gbm_mrd_no_mat_mean <- train(mrd ~ ., data = dat_no.na[,-grep("mat_mean|sr",names(dat_no.na))], method = "gbm", 
#                              trControl = gbmControl, verbose=FALSE,
#                              tuneGrid=gbmGrid)
# gbm_mrd_no_tra_mean <- train(mrd ~ ., data = dat_no.na[,-grep("tra_mean|sr",names(dat_no.na))], method = "gbm", 
#                              trControl = gbmControl, verbose=FALSE,
#                              tuneGrid=gbmGrid)
# 
# 
# gbmFit4_mrd$results[row.names(gbmFit4_mrd$bestTune),]
# gbm_mrd_no_tra_mean$results[row.names(gbm_mrd_no_tra_mean$bestTune),] # better
# gbm_mrd_no_mat_mean$results[row.names(gbm_mrd_no_mat_mean$bestTune),] # better
# 
# 
# # variable order
# ord <- as.data.frame(summary(gbmFit4_mrd))
# res2 <- as.data.frame(summary(gbm_mrd_no_tra_mean))
# res3 <- as.data.frame(summary(gbm_mrd_no_mat_mean))
# ord_sum <- merge(ord, res2, by="var", all=TRUE)
# ord_sum <- merge(ord_sum, res3, by="var", all=TRUE)
# ord_sum$no_mat_mean_model <- rank(ord_sum$rel.inf, ties.method = "last", na.last = "keep")
# ord_sum$full_model <- rank(ord_sum$rel.inf.x, ties.method = "last", na.last = "keep")
# ord_sum$no_tra_mean_model <- rank(ord_sum$rel.inf.y, ties.method = "last", na.last = "keep")
# 
# ord_sum2_mrd <- ord_sum %>% 
#   pivot_longer(cols=c(full_model, no_tra_mean_model, no_mat_mean_model))
# #ord_sum2$var2 <- abbreviate(ord_sum2$var, minlength = 7, named = F)
# ord_sum2_mrd$name <- factor(ord_sum2_mrd$name, levels = c("no_tra_mean_model", "full_model", "no_mat_mean_model"))
# 
# ord_sum2_mrd$col <- "grey"
# ord_sum2_mrd$col[ord_sum2_mrd$var=="mat_mean"] <- "red"
# ord_sum2_mrd$col[ord_sum2_mrd$var=="tra_mean"] <- "blue"
# ggplot(ord_sum2_mrd, aes(x=name, y=value, col=col))+
#   #  geom_point()+
#   geom_line(aes(group=var), col=ord_sum2_mrd$col)+
#   geom_label(aes(label=var), size=3, col=ord_sum2_mrd$col)+
#   scale_y_continuous("Variable rank")+
#   theme(legend.position = "none")
# #ggsave("results/mat_mean_vs_tra_mean_MRD.png", width=6, height=5, units = "in", dpi = 600)
# 
# # FACET PLOT
# ord_sum2$response <- "SR"
# ord_sum2_mrd$response <- "MRD"
# temp <- rbind(ord_sum2, ord_sum2_mrd)
# ggplot(temp, aes(x=name, y=value, group=var, col=col))+ 
#   geom_line()+
#   geom_label(aes(label=var), size=3)+
#   scale_y_continuous("Variable rank")+
#   scale_color_manual(values=c("blue", "grey", "red"))+
#   facet_wrap(~ response)+
#   theme(legend.position = "none", 
#         panel.spacing.x = unit(0, "mm"))
# ggsave("results/mat_mean_vs_tra_mean.png", width=9, height=5, units = "in", dpi = 600)
# 
# 
# # decision: 
# # opposed correlation between variables and SR/MRD
# # tiny changes for SR model. both models better with removing tra mean. highly correlated but cannot be coded in a latent variable, difficult
# 
# 
# 
# 
# 
# ##  tra_mean + subtrop mbf ##########################################################################################
# # correlation with SR and MRD
# cor.test(dat_no.na$sr, dat_no.na$sub_trop_mbf, method="s") # 0.47
# cor.test(dat_no.na$sr, dat_no.na$tra_mean, method="s") # -0.28
# cor.test(dat_no.na$mrd, dat_no.na$sub_trop_mbf, method="s") # -0.40
# cor.test(dat_no.na$mrd, dat_no.na$tra_mean, method="s") # 0.30
# 
# cor.test(dat_no.na$sub_trop_mbf, dat_no.na$tra_mean, method="s") # - 0.73
# 
# 
# # SR TEST"
# gbmFit4$results[row.names(gbmFit4$bestTune),]
# gbmFit5_no_trop_mbf$results[row.names(gbmFit5_no_trop_mbf$bestTune),] # worse
# gbmFit5_no_tra_mean$results[row.names(gbmFit5_no_tra_mean$bestTune),] # better
# 
# # variable order
# ord <- as.data.frame(summary(gbmFit4))
# res2 <- as.data.frame(summary(gbmFit5_no_tra_mean))
# res3 <- as.data.frame(summary(gbmFit5_no_trop_mbf))
# ord_sum <- merge(ord, res2, by="var", all=TRUE)
# ord_sum <- merge(ord_sum, res3, by="var", all=TRUE)
# 
# # plot relative influence in each model instead of order?
# test <- ord_sum[,c(1:4)]
# names(test) <- c("var", "full_model", "no_tra_mean_model", "no_trop_mbf_model")
# test2 <- test %>% 
#   pivot_longer(cols=c(full_model, no_tra_mean_model, no_trop_mbf_model))
# test2$col <- "other"
# test2$col[test2$var=="sub_trop_mbf"] <- "sub_trop_mbf"
# test2$col[test2$var=="tra_mean"] <- "tra_mean"
# test2$name <- factor(test2$name, levels = c("no_tra_mean_model", "full_model", "no_trop_mbf_model"))
# ggplot(test2, aes(x=name, y=value, group=var, col=col))+
#   geom_point()+
#   geom_line()+
#   scale_color_manual("Variable", values=c("red", "grey", "blue"))+
#   #geom_label(aes(label=var), size=3, col=test2$col)+
#   scale_y_log10("Relative influence")+
#   theme(legend.position = "right")
# # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 
# ord_sum$no_trop_mbf_model <- rank(ord_sum$rel.inf, ties.method = "last", na.last = "keep")
# ord_sum$full_model <- rank(ord_sum$rel.inf.x, ties.method = "last", na.last = "keep")
# ord_sum$no_tra_mean_model <- rank(ord_sum$rel.inf.y, ties.method = "last", na.last = "keep")
# 
# ord_sum2 <- ord_sum %>% 
#   pivot_longer(cols=c(full_model, no_tra_mean_model, no_trop_mbf_model))
# #ord_sum2$var2 <- abbreviate(ord_sum2$var, minlength = 7, named = F)
# ord_sum2$name <- factor(ord_sum2$name, levels = c("no_tra_mean_model", "full_model", "no_trop_mbf_model"))
# 
# ord_sum2$col <- "grey"
# ord_sum2$col[ord_sum2$var=="sub_trop_mbf"] <- "red"
# ord_sum2$col[ord_sum2$var=="tra_mean"] <- "blue"
# ggplot(ord_sum2, aes(x=name, y=value, col=col))+
#   #  geom_point()+
#   geom_line(aes(group=var), col=ord_sum2$col)+
#   geom_label(aes(label=var), size=3, col=ord_sum2$col)+
#   scale_y_continuous("Variable rank")+
#   theme(legend.position = "none")
# #ggsave("results/mat_mean_vs_tra_mean_SR.png", width=6, height=5, units = "in", dpi = 600)
# 
# # MRD TEST"
# gbmFit4_mrd$results[row.names(gbmFit4_mrd$bestTune),]
# gbm_mrd_no_tra_mean$results[row.names(gbm_mrd_no_tra_mean$bestTune),] # better
# gbm_mrd_no_mat_mean$results[row.names(gbm_mrd_no_trop_mbf$bestTune),] # worse
# 
# 
# # variable order
# ord <- as.data.frame(summary(gbmFit4_mrd))
# res2 <- as.data.frame(summary(gbm_mrd_no_tra_mean))
# res3 <- as.data.frame(summary(gbm_mrd_no_trop_mbf))
# ord_sum <- merge(ord, res2, by="var", all=TRUE)
# ord_sum <- merge(ord_sum, res3, by="var", all=TRUE)
# ord_sum$no_trop_mbf_model <- rank(ord_sum$rel.inf, ties.method = "last", na.last = "keep")
# ord_sum$full_model <- rank(ord_sum$rel.inf.x, ties.method = "last", na.last = "keep")
# ord_sum$no_tra_mean_model <- rank(ord_sum$rel.inf.y, ties.method = "last", na.last = "keep")
# 
# ord_sum2_mrd <- ord_sum %>% 
#   pivot_longer(cols=c(full_model, no_tra_mean_model, no_trop_mbf_model))
# #ord_sum2$var2 <- abbreviate(ord_sum2$var, minlength = 7, named = F)
# ord_sum2_mrd$name <- factor(ord_sum2_mrd$name, levels = c("no_tra_mean_model", "full_model", "no_trop_mbf_model"))
# 
# ord_sum2_mrd$col <- "grey"
# ord_sum2_mrd$col[ord_sum2_mrd$var=="sub_trop_mbf"] <- "red"
# ord_sum2_mrd$col[ord_sum2_mrd$var=="tra_mean"] <- "blue"
# ggplot(ord_sum2_mrd, aes(x=name, y=value, col=col))+
#   #  geom_point()+
#   geom_line(aes(group=var), col=ord_sum2_mrd$col)+
#   geom_label(aes(label=var), size=3, col=ord_sum2_mrd$col)+
#   scale_y_continuous("Variable rank")+
#   theme(legend.position = "none")
# #ggsave("results/mat_mean_vs_tra_mean_MRD.png", width=6, height=5, units = "in", dpi = 600)
# 
# # FACET PLOT
# ord_sum2$response <- "SR"
# ord_sum2_mrd$response <- "MRD"
# temp <- rbind(ord_sum2, ord_sum2_mrd)
# ggplot(temp, aes(x=name, y=value, group=var, col=col))+ 
#   geom_line()+
#   geom_label(aes(label=var), size=3)+
#   scale_y_continuous("Variable rank")+
#   scale_color_manual(values=c("blue", "grey", "red"))+
#   facet_wrap(~ response)+
#   theme(legend.position = "none", 
#         panel.spacing.x = unit(0, "mm"))
# ggsave("results/trop_mbf_vs_tra_mean.png", width=9, height=5, units = "in", dpi = 600)
# 
# 
# # decision: 
# # opposed correlation between variables and SR/MRD
# # tiny changes for both models
# 
# 
# 
# 
# 
# 
# 
# 
# # model performance summary table compared to full model ####
# library(knitr)
# m_names <- c("gbmFit4", "gbmFit5_no_trop_mbf", "gbmFit5_no_tra_mean",
#               "gbmFit5_no_tri", "gbmFit5_no_area", "gbmFit5_no_elev_mean", "gbmFit5_no_elev_sd", "gbmFit5_no_mat_mean",
#               "gbmFit5_no_mat_sd", "gbmFit5_no_number_biomes", "gbmFit5_no_pet_mean", "gbmFit5_no_pet_sd", 
#               "gbmFit5_no_pre_mean", "gbmFit5_no_sea_sd", "gbmFit5_no_soil", "gbmFit5_no_tra_sd")
# #m_list <- list(mget(m_names))
# m_list <- list(gbmFit4, gbmFit5_no_trop_mbf, gbmFit5_no_tra_mean,
#                gbmFit5_no_tri, gbmFit5_no_area, gbmFit5_no_elev_mean, gbmFit5_no_elev_sd, gbmFit5_no_mat_mean,
#                gbmFit5_no_mat_sd, gbmFit5_no_number_biomes, gbmFit5_no_pet_mean, gbmFit5_no_pet_sd, 
#                gbmFit5_no_pre_mean, gbmFit5_no_sea_sd, gbmFit5_no_soil, gbmFit5_no_tra_sd)
# for(i in 1:length(m_names)){
#   if(i==1)
#     df <- m_list[[i]]$results[row.names(m_list[[i]]$bestTune),c(5:6)]
#   else 
#     df <- rbind(df, m_list[[i]]$results[row.names(m_list[[i]]$bestTune),c(5:6)])
#   print(m_list[[i]]$results[row.names(m_list[[i]]$bestTune),c(5:6)])
# }
# df$model_name <- (m_names)
# kable(df, format = "html", col.names = c("RMSE", "R squared", "Model name"))
# 
# 
# 
# 
# 
# # remove selected variables:
# #dat_no.na <- dat_no.na[,-grep("number_biomes", names(dat_no.na))]
# #dat_no.na <- dat_no.na[,-grep("elev_sd", names(dat_no.na))]
# 
# 
# 
# 
# 
# 
# 
# ##### SEM LAVAAN GLOBAL ESTIMATION ##########################################################
# library(lavaan)
# library(semPlot)
# library(tidySEM)
# library(fitdistrplus)
# library(modEvA)
# source("scripts/clean_semplot_functions.R") # script to remove unsignificant edges for plotting
# 
# ## distributions ##
# hist(dat_no.na[,c("sr", "mrd")])
# hist(dat_no.na$sr, breaks=100)
# dev.off()
# # mrd is kinda normally distributed, SR not. Lets test:
# descdist(dat_no.na$mrd, discrete = FALSE, boot=500)
# shapiro.test(dat_no.na$mrd) # not far away from normal distribution
# mrd_norm_fit <- fitdist(dat_no.na$mrd, "norm",method="mle")
# summary(mrd_norm_fit)
# plot(mrd_norm_fit)
# 
# descdist(dat_no.na$sr, discrete = TRUE, boot=500) # Poisson regression or negative binomial regression
# descdist(log(dat_no.na$sr), discrete = FALSE, boot=500) # doesnt make it better
# descdist(sqrt(dat_no.na$sr), discrete = FALSE, boot=500) # lognormal?
# fitp <- fitdist(dat_no.na$sr,"pois") # untransformed poisson
# fnbinom <- fitdist(dat_no.na$sr,"nbinom", discrete = TRUE) # untransformed neg.binomial
# fnorm_nt <- fitdist(dat_no.na$sr,"norm") # untransformed normal
# flnorm <- fitdist(sqrt(dat_no.na$sr),"lnorm")
# fnorm <- fitdist(sqrt(dat_no.na$sr),"norm")
# 
# summary(fitp)
# summary(fnbinom)
# summary(fnorm_nt)
# #_____________
# summary(flnorm)
# summary(fnorm)
# 
# plot(fitp) # numbers too large
# plot(fnbinom)
# plot(fnorm_nt)
# #_____________
# plot(flnorm)
# plot(fnorm)
# 
# # test with soil + area
# lmp <- lm(sr~soil + area, dat_no.na)
# lms <- lm(sqrt(sr)~soil + area, dat_no.na)
# m1 <- glm.nb(sr ~ soil + area, dat_no.na)
# hist(lmp$residuals, breaks=20)
# hist(lms$residuals, breaks=20)
# hist(m1$residuals, breaks=20)
# modEvA::RsqGLM(m1)
# 
# # Decision: Raw data has best fit for negative binomial distribution.
# ## AIC between sqrt() transformed lnorm and norm only show minor differences
# ## fitting a simple model for SR with soil and area shows best residual distribution for the sqrt() transformed noemal model, the Rsquared is comparable with the neg binomial model pseudo R squared.
# # --> sqrt() transform the species richness data, keep the global model estimation
# 
# dat_no.na$sr_trans <- sqrt(dat_no.na$sr)
# 
# # scale data?
# ## standardize all variables with z-transformation to avoid variance issues
# ztrans <- function(x){(x - mean(x)) / sd(x)}
# 
# dat_no.na <- apply(dat_no.na, 2, ztrans)
# dat_no.na <- as.data.frame(dat_no.na)
# 
# # rename all _mean to _m for better visibility in der model figures
# names(dat_no.na) <- sub("_mean$", "_m", names(dat_no.na))
# 
# # model building procedure############
# # Models were build starting off with no latent variables and 2 regression models including all variables 
# # with relative influence >5% as estimated in the exploratory analysis full GBMs (see figure S2).
# # We subsequently added more variables, following the their importance rank es estimated in the full GBMS, until the model stopped matching the acceptable parameters indicating a good fit (CFI > 0.9, RMSEA < 0.08, df>0, chisq p-value > 0.05). 
# 
# # simple model including all variables with importance > 5% in GBMs (and mrd) ####################################
# mod_base5 <- "
#           # regressions
#           sr_trans ~ soil 
#           + sub_trop_mbf 
#           + sub_trop_dbf
#           + area
#           + pre_sd 
#           + mrd 
#           
#           mrd ~ pre_m
#           + sea_m
#           + tra_m
#           + sea_sd
#           + soil 
#           + tri
# "
# fsem_base5 <- sem(mod_base5, data = dat_no.na)#, se = "bootstrap", bootstrap = 100) 
# summary(fsem_base5, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
# # adding significance to edges
# semPaths(object = fsem_base5, layout = "circle2", rotation = 1,
#          whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE,
#          edgeLabels = sem_sig_labels(fsem_base5))
# fitmeasures(fsem_base5, c("chisq", "df", "pvalue", "cfi", "rmsea", "aic"))
# modificationindices(fsem_base5, sort = TRUE)
# # fsem_base5_mod <- clean.semplot(fsem_base5)
# # semPaths(object = fsem_base5, layout = "circle2", rotation = 1,
# #          whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
# #          exoCov = FALSE, exoVar = FALSE,
# #          nCharNodes = 0, fade=FALSE)
# # good model. several non-significant edges
# # remove non signifant edges ####
# mod_base5_ns <- "
#           # regressions
#           sr_trans ~ soil 
#           + sub_trop_mbf 
#           #+ sub_trop_dbf
#           #+ area
#           + pre_sd 
#           + mrd 
#           
#           mrd ~ pre_m
#           #+ sea_m
#           #+ tri
#           + sea_sd 
#           + soil 
# "
# fsem_base5_ns <- sem(mod_base5_ns, data = dat_no.na)
# summary(fsem_base5_ns, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
# # adding significance to edges
# semPaths(object = fsem_base5_ns, layout = "circle2", rotation = 1,
#          whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE,
#          edgeLabels = sem_sig_labels(fsem_base5_ns))
# fitmeasures(fsem_base5_ns, c("chisq", "df", "pvalue", "cfi", "rmsea", "aic"))
# anova(fsem_base5_ns,fsem_base5)
# # slightly better, but as we are trying to show all variables influence, not worth removing them
# 
# 
# # adding more variables --> forward selection ###################################################################
# ## - forward selection based on AIC (while keeping the other goodness of fit values under control) 
# mod_base5_plus <- "
#           # regressions
#           sr_trans ~ soil 
#           + sub_trop_mbf 
#           + sub_trop_dbf
#           + area
#           + pre_sd 
#           + mrd
#           + mont_gs 
#           + pre_m
#           + sea_sd
#           + tra_m
#           + elev_m
#           + mat_m
#           #+ tra_sd
#           #+ mat_sd
# 
#           mrd ~ pre_m
#           + sea_m
#           + tra_m
#           + sea_sd
#           + soil 
#           + tri
#           + tra_sd 
#           + pre_sd
#           + elev_m
#           + pet_m
#           + mat_m
#           # + area
#           #+ mat_sd
#           
#           # modifications
#           mrd ~ sub_trop_mbf + mont_gs
#           sr_trans ~ pet_m + sea_m
#           
# "
# fsem_base5_plus <- sem(mod_base5_plus, data = dat_no.na)
# fitmeasures(fsem_base5_plus, c("chisq", "df", "pvalue", "cfi", "rmsea", "aic"))
# summary(fsem_base5_plus, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
# semPaths(object = fsem_base5_plus, layout = "circle2", rotation = 1,
#          whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE,
#          edgeLabels = sem_sig_labels(fsem_base5_plus))
# modificationindices(fsem_base5_plus, sort = TRUE)[1:10,]
# # good model, as many variables as possible in this configuration model fit becomes bad (AIC stays about the same)
# 
# # FULL MODEL #####
# mod_full <- "
#           # regressions # all variables with releative influence>1
#           sr_trans ~ soil + sub_trop_mbf + sub_trop_dbf+ area+ pre_sd+ mrd+ mont_gs + pre_m+ sea_sd+ tra_m
#           + elev_m+ mat_m +tra_sd+ mat_sd + sea_m
#             # all variables with releative influence>2
#           mrd ~ pre_m+ sea_m+ tra_m+ sea_sd+ soil + tri+ tra_sd + pre_sd+ elev_m+ pet_m+ mat_m + area
#           + mat_sd + temp_bmf + sub_trop_mbf + pet_sd
#           
#           # modifications
#           sr_trans ~ temp_bmf + pet_m
#           mrd ~ mont_gs
# "
# fsem_full <- sem(mod_full, data = dat_no.na,  se = "bootstrap", bootstrap = 500)
# summary(fsem_full, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
# 
# 
# # check validity of regressions - SR ####################
# lm_full <- lm(sr_trans ~ soil + sub_trop_mbf + sub_trop_dbf+ area+ pre_sd+ mrd+ mont_gs + pre_m+ sea_sd+ tra_m
# + elev_m+ mat_m +tra_sd+ mat_sd + sea_m, data=dat_no.na)
# summary(lm_full)
# car::vif(lm_full)
# step.sr <- step(lm_full)
# step.sr$coefficients
# lm_full_red <- lm(sr_trans ~ soil + sub_trop_mbf + area + pre_sd +  pre_m + sea_sd + tra_m + elev_m + mat_m + mat_sd
#                   + sea_m+  mrd, data = dat_no.na) # removed: mrd, mont_gs, tra_sd, sub_trop_dbf
# car::vif(lm_full_red)
# 
# # remove mat_m
# lm_full_nomat <- lm(sr_trans ~ soil + sub_trop_mbf + sub_trop_dbf+ area+ pre_sd+ mrd+ mont_gs + pre_m+ sea_sd+ tra_m
#               + elev_m+tra_sd+ mat_sd + sea_m, data=dat_no.na)
# summary(lm_full_nomat)
# car::vif(lm_full_nomat) # all issues solved
# step.sr.nomat <- step(lm_full_nomat)
# summary(step.sr.nomat)
# lm_full_red2 <- lm(formula = sr_trans ~ soil + sub_trop_mbf + area + pre_sd + mrd + sea_sd + tra_m + mat_sd, 
#                    data = dat_no.na)
# summary(lm_full_red2)
# 
# # check validity of regressions - MRD ####################
# lm_mrd <- lm(mrd ~ pre_m+ sea_m+ tra_m+ sea_sd+ soil + tri+ tra_sd + pre_sd+ elev_m+ pet_m + area
# + mat_sd + temp_bmf + sub_trop_mbf + pet_sd, data = dat_no.na)
# car::vif(lm_mrd) # mat_m + pet_m cause massive collinearity, removing mat solves the problem
# step(lm_mrd) # removes pre_sd, tri, pet_sd, sea_m, mad_sd, area, 
# lm_mrd_red <- lm(formula = mrd ~ pre_m + tra_m + sea_sd + soil + tra_sd + elev_m + 
#      pet_m + temp_bmf + sub_trop_mbf, data = dat_no.na)
# summary(lm_mrd_red)
# 
# # MODEL BASED ON STEPWISE AIC SELECTION #######
# mod_aic <- "
#             # regressions
#             sr_trans ~ soil + sub_trop_mbf + area + pre_sd + mrd + sea_sd + tra_m + mat_sd
#             mrd ~ pre_m + tra_m + sea_sd + soil + tra_sd + elev_m + pet_m + temp_bmf + sub_trop_mbf
#             # indirect effects
#             sub_trop_mbf ~ pre_m + tra_m + pet_m
# #            tra_sd ~ area # thats the strongest correlation
# #            sea_sd ~ area
# #            soil ~ area
# #            pre_sd ~ area
# #            mat_sd ~ area
#             # mods
#             sr_trans ~ temp_bmf
# "
# fsem_aic <- sem(mod_aic, data = dat_no.na)
# fitmeasures(fsem_aic, c("chisq", "df", "pvalue", "cfi", "rmsea", "aic"))
# summary(fsem_aic, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
# semPaths(object = fsem_aic, layout = "circle2", rotation = 1,
#          whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE,
#          edgeLabels = sem_sig_labels(fsem_aic))
# modificationindices(fsem_aic, sort = TRUE)[1:15,]
# 
# 
# # now backwards selection: remove non-significant variables from full model (AIC decreases) #####
# mod_full_bs <- "
#           # regressions
#           sr_trans ~ soil + sub_trop_mbf + pre_sd+ pre_m+ tra_m +
#            mat_m + mat_sd + sea_m + mrd + area
#             
#           mrd ~ pre_m+ tra_m+ sea_sd+ soil + tra_sd + mat_m
#            + temp_bmf + sub_trop_mbf
#           
#           # modifications
#           sr_trans ~ temp_bmf + pet_m
#           mrd ~ mont_gs
#           
#           #indirect area effects
# #          soil ~ area
# "
# fsem_full_bs <- sem(mod_full_bs, data = dat_no.na)
# fitmeasures(fsem_full_bs, c("chisq", "df", "pvalue", "cfi", "rmsea", "aic"))
# summary(fsem_full_bs, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
# semPaths(object = fsem_full_bs, layout = "circle2", rotation = 1,
#          whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE,
#          edgeLabels = sem_sig_labels(fsem_full_bs))
# modificationindices(fsem_full_bs, sort = TRUE)[1:15,]
# 
# mod_full_sr <- "
#           # regressions # all variables with releative influence>1
#           sr_trans ~ soil + sub_trop_mbf + sub_trop_dbf+ area+ pre_sd+ mrd+ mont_gs + sea_sd
#           + elev_m+tra_sd+ mat_sd
#           
#           # area influence on spatial heterogeneity
#           soil ~ area + elev_m
#           #pre_sd ~ area
#           #tra_sd ~ area
#           #sea_sd ~ area
#           #mat_sd ~ area
#           
#           # climate influence on habitat
#           #sub_trop_mbf ~ pre_m + pet_m + tra_m  + mat_m
#           #sub_trop_dbf ~ pre_m + pet_m + tra_m  + mat_m
# "
# fsem_full_sr <- sem(mod_full_sr, data = dat_no.na)
# fitmeasures(fsem_full_sr, c("chisq", "df", "pvalue", "cfi", "rmsea", "aic"))
# summary(fsem_full_sr, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
# modificationindices(fsem_full_sr, sort = TRUE)[1:10,]
# 
# # adding indirect effects  ######################################################
# # base model #####
# mod_base5_ind <- "
#           # regressions
#           sr_trans ~ a*soil 
#           + sub_trop_mbf 
#           + sub_trop_dbf
#           + area
#           + b*pre_sd 
#           + mrd 
#           
#           mrd ~ pre_m
#           + sea_m
#           + tri
#           + sea_sd 
#           + soil 
#           
#           soil ~ a1*area
#           #pre_sd ~ a2*area
#           #pre_sd ~ p1*pre_m
#           
#           area_ind := a1*a #+ a2*b #(via soil and pre_sd)
#           
# "
# fsem_base5_ind <- sem(mod_base5_ind, data = dat_no.na) 
# summary(fsem_base5_ind, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
# # adding significance to edges
# semPaths(object = fsem_base5_ind, layout = "circle2", rotation = 1,
#          whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE,
#          edgeLabels = sem_sig_labels(fsem_base5_ind))
# fitmeasures(fsem_base5_ind, c("chisq", "df", "pvalue", "cfi", "rmsea", "aic"))
# modificationindices(fsem_base5_ind, sort = TRUE)
# 
# # compare to without indirect effects:
# anova(fsem_base5, fsem_base5_ind)
# # adding area effect on soil makes model worse although it makes biologically sense?
# 
# 
# # advanced model with more variables ####
# # all indirect effects are based on theoretical expectations!
# mod_full_ind <- "
#           # regressions
#           sr_trans ~ soil + sub_trop_mbf + sub_trop_dbf+ area+ pre_sd+ mrd+ mont_gs + pre_m+ sea_sd+ tra_m
#           + elev_m+ mat_m +tra_sd+ mat_sd + sea_m 
#           + temp_bmf + pet_m
#           
#           mrd ~ pre_m+ sea_m+ tra_m+ sea_sd+ soil + tri+ tra_sd + pre_sd+ elev_m+ pet_m+ mat_m + area
#           + mat_sd + temp_bmf + sub_trop_mbf + pet_sd + 
#           mont_gs
#           
#           # area influence on spatial heterogeneity
#           soil ~ area
#           #pre_sd ~ area
#           #tra_sd ~ area
#           #sea_sd ~ area
#           #mat_sd ~ area
# 
#           # climate influence on habitat
#           #sub_trop_mbf ~ pre_m + pet_m + tra_m  + mat_m
#           #sub_trop_dbf ~ pre_m + pet_m + tra_m  + mat_m
#           
#           # modifications
# #          mat_sd~elev_m
# #          sea_sd ~ sea_m
# #          pre_sd ~ pre_m
# #          tra_sd ~ tra_m
# #          soil ~~ pre_sd
# #          pre_sd ~~ sea_sd
# #          sub_trop_mbf ~ temp_bmf
# "
# fsem_full_ind <- sem(mod_full_ind, data = dat_no.na)
# fitmeasures(fsem_full_ind, c("chisq", "df", "pvalue", "cfi", "rmsea", "aic"))
# modificationindices(fsem_full_ind, sort = TRUE)[1:10,]
# 
# # try adding to base model first 
# mod_base5_ind <- "
#           # regressions
#           sr_trans ~ soil + sub_trop_mbf + sub_trop_dbf+ area+ pre_sd + mrd 
#   
#           mrd ~ pre_m+ sea_m+ tra_m + sea_sd + soil + tri
#           
#           # area influence on spatial heterogeneity
# #          soil ~ area
#           #pre_sd ~ area
#           #tra_sd ~ area
#           #sea_sd ~ area
#           #mat_sd ~ area
# 
#           # climate influence on habitat
# #          sub_trop_mbf ~ pre_m + tra_m #+ pet_m+ mat_m 
# #          sub_trop_dbf ~ pre_m + pet_m + tra_m  + mat_m
# "
# fsem_base5_ind <- sem(mod_base5_ind, data = dat_no.na)#, se = "bootstrap", bootstrap = 100) 
# fitmeasures(fsem_base5_ind, c("chisq", "df", "pvalue", "cfi", "rmsea", "aic"))
# modificationindices(fsem_base5_ind, sort = TRUE)[1:10,]
# 
# # 
# # mod_base5_plus_ind <- "
# #           # regressions
# #           sr_trans ~ soil + sub_trop_mbf + sub_trop_dbf + area + pre_sd + mont_gs + pre_m + mrd + tri + tra_sd
# #           + elev_m + sea_m + tra_m + pet_m
# #           
# #           mrd ~ pre_m + sea_m + tri + sea_sd + soil + tra_sd + pre_sd + elev_m + tra_m + sub_trop_mbf + pet_m + mont_gs
# #           
# #           # area influence on spatial heterogeneity
# #           soil ~ area
# #           pre_sd ~ area
# #           tra_sd ~ area
# #           sea_sd ~ area
# #           
# #           # climate influence on habitat
# #           sub_trop_mbf ~ pre_m + pet_m + tra_m  # mat would go here but is not part of the model
# #           sub_trop_dbf ~ pre_m + pet_m + tra_m  # mat would go here but is not part of the model
# #           
# #           # modifications
# #           sea_sd ~ sea_m
# #           pre_sd ~ pre_m
# #           tra_sd ~ tra_m
# #           soil ~~ pre_sd    # uncausated correlation
# #           pre_sd ~~ sea_sd  # correlation with a cause?
# #           tra_sd ~ sea_sd
# #           sea_sd ~ tra_sd
# # "
# # fsem_ind2 <- sem(mod_base5_plus_ind, data = dat_no.na)
# # fitmeasures(fsem_ind2, c("chisq", "df", "pvalue", "cfi", "rmsea", "aic"))
# # modificationindices(fsem_ind2, sort = TRUE)[c(1:15),]
# # summary(fsem_ind2, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
# # 
# # # adding significance to edges
# # semPaths(object = fsem_ind2, layout = "circle2", rotation = 1, whatLabels="std", label.cex = 2, 
# #          edge.label.cex = 1.3, what = "std", exoCov = FALSE, exoVar = FALSE, nCharNodes = 0, 
# #          fade=FALSE, edgeLabels = sem_sig_labels(fsem_ind1))
# # 
# # 
# # mod_ind_area <- "
# #           # regressions
# #           sr_trans ~ a*soil 
# #           + b*sub_trop_mbf 
# #           + c*sub_trop_dbf
# #           + d*area
# #           + e*pre_sd 
# #           + f*mrd
# #           + l*pre_m
# #           pre_sd ~ a1*area
# #           soil ~ a2*area
# #           
# #           mrd ~ g*pre_m + 
# #           + h*sea_m 
# #           + i*tri 
# #           + j*sea_sd 
# #           + k*soil 
# # #          sub_trop_mbf ~ p1*pre_m + sea_m# + sea_sd
# # #          sub_trop_dbf ~ p2*pre_m
# #           tri ~ a3*area
# #           sea_sd ~ a4*area
# #           
# # #          pre_m_via_trop_dbf := p1*c
# # #          pre_m_via_trop_mbf := p2*b
# # #          sea_m_via_mrd := h*f
# # #          soil_via_mrd := k*f
# # #          tri_via_mrd := i*f
# # #          sea_sd_via_mrd := j*f
# #           
# #           ind_area_on_sr := a1*a + a2*e + a3*i*f + a4*j*f
# # #          ind_soil_on_sr := d*f
# #           
# #           ind_pre_m_on_sr := g*f + p1*c + p2*b
# #           area_on_sr_all := ind_area_on_sr + d
# #           
# #           # mods
# #           # sea_sd  ~        sea_m
# #           # pre_sd  ~        pre_m 
# #           # pre_sd ~~         soil # covariance is non causal
# #           # pre_sd ~~       sea_sd # covariance is non causal?
# #           # tri  ~       pre_sd
# #           # sub_trop_dbf  ~        sea_m
# #           # sub_trop_mbf  ~       sea_sd
# #           
# # "
# # 
# # fsem_ind1 <- sem(mod_ind1, data = dat_no.na)
# # fitmeasures(fsem_ind1, c("chisq", "df", "pvalue", "cfi", "rmsea", "aic"))
# # modificationindices(fsem_ind1, sort = TRUE)
# # summary(fsem_ind1, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
# # 
# # # adding significance to edges
# # semPaths(object = fsem_ind1, layout = "circle2", rotation = 1,
# #          whatLabels="std", label.cex = 2, edge.label.cex = 1.3, what = "std", 
# #          exoCov = FALSE, exoVar = FALSE,
# #          nCharNodes = 0, fade=FALSE
# #          , edgeLabels = sem_sig_labels(fsem_ind1))
# # 
# # # modifications make model acceptable. Rsquared SR is 0.65, thats rather good
# # ## modify model plot to only show significant paths
# # fsem_ind1_mod <- clean.semplot(fsem_ind1)
# # semPaths(object = fsem_ind1_mod, layout = "circle2", rotation = 1,
# #          whatLabels="std", label.cex = 2, edge.label.cex = 1, what = "std", 
# #          exoCov = FALSE, exoVar = FALSE,
# #          nCharNodes = 0, fade=FALSE)
# # #         , edgeLabels = b)
# 
# 
# 
# 
# # FULL VS PARTIAL MEDIATION
# # indirect effect of climate via sub_trop_mbf ####
# library(car) # check for multicollinearity
# 
# cor(dat_no.na[,c("sr_trans", "sub_trop_mbf", "tra_m", "pet_m", "pre_m", "mat_m")])
# summary(lm(sr_trans~ sub_trop_mbf, data=dat_no.na))
# summary(lm(sr_trans~ tra_m, data=dat_no.na))
# summary(lm_clim <- lm(sr_trans~ pre_m + sub_trop_mbf + mat_m + pet_m + tra_m, data=dat_no.na))
# # direct effect dominated by habitat and pre
# car::vif(lm_clim)
# summary(lm_clim)
# summary(step(lm_clim)) # keep only pre_m, pet_m, subtrop_mbf for SR direct regression
# 
# summary(lm_sub_trop_mbf <- lm(sub_trop_mbf~ pre_m  + mat_m + tra_m + pet_m, data=dat_no.na)) # tra_m is not significant despite high correlation, why?
# car::vif(lm_sub_trop_mbf)
# summary(lm(sub_trop_mbf~ pre_m  + tra_m + pet_m, data=dat_no.na))
# summary(lm(sub_trop_mbf~ pre_m  + tra_m, data=dat_no.na))
# summary(step(lm_sub_trop_mbf)) # stepwise AIC # keep mat, pre and pet for subtrop_mbf regression
# car::vif(lm(sub_trop_mbf~ pre_m + mat_m + pet_m, data=dat_no.na))
# # --> tra_m causes problems in the model, removing it is suggested via stepwiseAIC model selection and leaves the remaining model with better almost identical Rsquared while keeping all variables significant with no multicolinearity issues.
# # --> selecting via multicolinearity tells us to remove mat_m, which gives a samller RSquared and keeps only 2 variables.
# 
# 
# #################################science bit ########################
# dat <- dat_no.na[,c("sr_trans", "sub_trop_mbf", "pre_m")]
# names(dat) <- c("sr", "trf", "pre")
# cor(dat)
# lm1 <- lm(sr~trf, dat)
# lm2 <- lm(sr~pre, dat)
# lm1$coefficients
# lm2$coefficients
# lm3 <- lm(sr~trf+pre, dat)
# summary(lm3)
# mod_trf_pre <- "
#   sr ~ a*trf + b*pre
#   trf ~ c*pre
#   pre_ind := c*a
#   pre_tot := c*a + b
# "
# sem_trf_pre_bs500 <- sem(mod_trf_pre, data = dat, test = "bootstrap", bootstrap = 500)
# #sem_trf_pre <- sem(mod_trf_pre, data = dat)
# parameterestimates(sem_trf_pre_bs500, boot.ci.type = "bca.simp", standardized = TRUE)
# fitmeasures(sem_trf_pre_bs500, c("chisq", "df", "pvalue","cfi", "rmsea", "aic"))
# summary(sem_trf_pre_bs500, standardized = TRUE, rsq=TRUE)
# semPaths(object = sem_trf_pre_bs500, layout = "circle", rotation = 1,
#          whatLabels="std", label.cex = 1.3, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE
#          , edgeLabels = sem_sig_labels(sem_trf_pre_bs500))
# #####################################################################
# 
# 
# mod_clim <- "
#       sr_trans ~ a*pre_m + b*sub_trop_mbf + c*pet_m
#       sub_trop_mbf ~ d*pre_m + e*mat_m + f*pet_m
#       ind_pre := d*b
#       ind_mat := e*b
#       ind_pet := f*b
#       total_pre := d*b + a
#       total_pet := f*b + c
# "
# mod_clim_full <- "
#       sr_trans ~ a*pre_m + b*sub_trop_mbf + c*mat_m + d*pet_m + e*tra_m
#       sub_trop_mbf ~ f*pre_m + g*mat_m  + h*pet_m + i*tra_m
#       ind_pre := f*b
#       ind_mat := g*b
#       ind_pet := h*b
#       ind_tra := i*b
#       total_pre := c*b + a
#       total_mat := g*b + c
#       total_pet := h*b + d
#       total_tra := i*b + e
# "
# fsem_clim_full <- sem(mod_clim_full, data = dat_no.na)
# fsem_clim_boot <- sem(mod_clim_full, data = dat_no.na,  se = "bootstrap", bootstrap = 1000)
# fitmeasures(fsem_clim_full, c("chisq", "df", "pvalue", "cfi", "rmsea", "aic"))
# summary(fsem_clim_full, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
# summary(fsem_clim_boot, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
# fsem_clim <- sem(mod_clim, data = dat_no.na)
# fitmeasures(fsem_clim, c("chisq", "df", "pvalue","cfi", "rmsea", "aic"))
# summary(fsem_clim, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
# semPaths(object = fsem_clim, layout = "circle2", rotation = 1,
#          whatLabels="std", label.cex = 2, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          , edgeLabels = sem_sig_labels(fsem_clim))
# semPaths(object = fsem_clim_full, layout = "circle2", rotation = 1,
#          whatLabels="std", label.cex = 2, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE
#          , edgeLabels = sem_sig_labels(fsem_clim_full))
# resid(fsem_clim_full, "cor")
# # fitted(fsem_clim) # the model-implied covariance matrix
# # inspect(fsem_clim, what="cor.all") # the model-implied correlation among variables
# # lavCor(fsem_clim) # the observed correlations
# modificationindices(fsem_clim,sort=TRUE)[1:10,] # pet_m has to stay as direct effect
# 
# 
# # # removing direct effects of mat, pet, tra
# # mod_clim_sub <- "
# #       sr_trans ~ a*pre_m + b*sub_trop_mbf
# #       sub_trop_mbf ~ f*pre_m + g*mat_m + h*pet_m + i*tra_m
# #       ind_pre := f*b
# #       ind_mat := g*b
# #       ind_pet := h*b
# #       ind_tra := i*b
# #       total_pre := f*b + a
# # "
# # fsem_clim_sub <- sem(mod_clim_sub, data = dat_no.na)
# # fitmeasures(fsem_clim_sub, c("chisq", "df", "pvalue", "cfi", "rmsea", "aic"))
# # summary(fsem_clim_sub, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
# # anova(fsem_clim_sub, fsem_clim)
# # # full model is better with partial mediation, even though the infleunces are not significant
# # 
# # fitted(fsem_clim_sub) # the model-implied covariance matrix
# # inspect(fsem_clim_sub, what="cor.all") # the model-implied correlation among variables
# # lavCor(fsem_clim) # the observed correlations
# # resid(fsem_clim_sub, "cor") # positive values: model unterpredicts correlation, negative vice versa. values > |.1| are woth a look
# # # the model underpredicts the correlation of SR with mat_m and pet_m
# # modificationindices(fsem_clim_sub,sort=TRUE)[1:10,] # pet_m has to stay as direct effect
# 
# # SUBSET MODELS
# # indirect effect of pre_m via sub_trop_mbf ####
# mod_pre_pm <- "  # partial mediation model (direct and indirect effects)
#       sr_trans ~ a*pre_m + b*sub_trop_mbf
#       sub_trop_mbf ~ c*pre_m
#       ind_pre := c*b
#       total_pre := c*b + a
# "
# fsem_pre_pm <- sem(mod_pre_pm, data = dat_no.na)
# mod_pre_fm <- "  # full mediation model (only indirect effects)
#       sr_trans ~ a*sub_trop_mbf
#       sub_trop_mbf ~ b*pre_m
#       ind_pre := a*b
# "
# fsem_pre_fm <- sem(mod_pre_fm, data = dat_no.na)
# summary(fsem_pre_fm, rsq=TRUE)
# # note: fsem_pre_pm has no DF left, the model is saturated and we cannot get a chisq statistic! but we can use anova to compare both models
# anova(fsem_pre_fm, fsem_pre_pm)
# # The models signif differ and the partial mediation model is preferred (lower AIC). If the difference was not significant, one would prefer the more parsimonious model and conclude that there is no direct effect of pre_m on SR.
# 
# 
# # indirect effect of tra_m via sub_trop_mbf ####
# mod_tra_pm <- "  # partial mediation model (direct and indirect effects)
#       sr ~ a*tra_m + b*sub_trop_mbf
#       sub_trop_mbf ~ c*tra_m
#       ind_pre := c*b
#       total_pre := c*b + a
# "
# fsem_tra_pm <- sem(mod_tra_pm, data = dat_no.na)
# mod_tra_fm <- "  # full mediation model (only indirect effects)
#       sr ~ a*sub_trop_mbf
#       sub_trop_mbf ~ b*tra_m
#       ind_pre := a*b
# "
# fsem_tra_fm <- sem(mod_tra_fm, data = dat_no.na)
# summary(fsem_tra_pm, rsq=TRUE, fit.measures=TRUE)
# anova(fsem_tra_fm, fsem_tra_pm)
# # models differ, pm is preferred because smaller AIC
# 
# 
# # indirect effect of mat_m via sub_trop_mbf ####
# mod_mat_pm <- "  # partial mediation model (direct and indirect effects)
#       sr ~ a*mat_m + b*sub_trop_mbf
#       sub_trop_mbf ~ c*mat_m
#       ind_mat := c*b
#       total_mat := c*b + a
# "
# fsem_mat_pm <- sem(mod_mat_pm, data = dat_no.na)
# mod_mat_fm <- "  # full mediation model (only indirect effects)
#       sr ~ a*sub_trop_mbf
#       sub_trop_mbf ~ b*mat_m
#       ind_mat := a*b
# "
# fsem_mat_fm <- sem(mod_mat_fm, data = dat_no.na)
# summary(fsem_mat_pm, rsq=TRUE, fit.measures=TRUE)
# anova(fsem_mat_fm, fsem_mat_pm)
# # models differ, pm is preferred because smaller AIC
# 
# 
# ###########################################################################################3
# 
# # make the good model better: add indirect area effects
# mod_base5_plus_ia <- "
#           # regressions
#           sr ~ a*soil + sub_trop_mbf + sub_trop_dbf + c*pre_sd + mont_gs + pre_m  + b*area 
#           + mrd + d*tri + e*tra_sd + elev_m + f*sea_m + g*tra_m
#           
#           mrd ~ pre_m + sea_m + tri + sea_sd + soil + tra_sd + pre_sd + elev_m 
#           + tra_m + sub_trop_mbf + pet_m + mont_gs
# 
#           #soil ~ a1*area          
#           #pre_sd ~ a2*area
#           #tri ~ a3*area
#           #tra_sd ~ a4*area
#           #sea_m ~ a5*area
#           #tra_m ~ a6*area
#           
#           #area_ind := a1*a   #+a4*e +a5*f #+a6*g #+ a2*c +a3*d
#           #area_tot := b +a1*a #  +a4*e +a5*f #+a6*g  #+ a2*c  +a3*d
# 
# "
# fsem_base5_plus_ia <- sem(mod_base5_plus_ia, data = dat_no.na)
# summary(fsem_base5_plus_ia, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
# fitmeasures(fsem_base5_plus_ia, c("chisq", "df", "pvalue", "cfi", "rmsea", "aic"))
# semPaths(object = fsem_base5_plus_ia, layout = "circle2", rotation =1,
#          whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE,
#          edgeLabels = sem_sig_labels(fsem_base5_plus_ia))
# anova(fsem_base5_plus, fsem_base5_plus_ia)
# # adding just one indirect area effect makes the model fit much worse - maybe it`s because the model to start with was overfitted??`
# 
# 
# 
# 
# 
# # CHECK AREA MEDIATION (via soil etc) #############################################
# cor(dat_no.na$sr_trans, dat_no.na$area)
# cov(dat_no.na$sr_trans, dat_no.na$area) # thats because we scaled the data
# summary(lm(sr_trans~area, data=dat_no.na))
# 
# summary(lm(sr_trans~area + soil, data=dat_no.na))
# # this says area has no effect on species richness, although we know its correlated (and it makes sense)
# mod_sasr_simp <- "
#   sr_trans ~ a*area + b*soil
# "
# mod_sasr_pm <- "
#   sr_trans ~ a*area + b*soil
#   soil ~ c*area
#   indirect := b*c
#   total := b*c + a
# "
# mod_sasr_fm <- "
#   sr_trans ~ b*soil
#   soil ~ c*area
#   indirect := b*c
# "
# sasr_simp <- sem(mod_sasr_simp, data = dat_no.na)
# sasr_pm <- sem(mod_sasr_pm, data = dat_no.na)
# sasr_fm <- sem(mod_sasr_fm, data = dat_no.na)
# 
# # MODEL DIAGNOSTICS
# fitted(sasr_pm) # the model-implied covariance matrix
# inspect(sasr_pm, what="cor.all") # the model-implied correlation among variables
# lavCor(sasr_pm) # the observed correlations
# resid(sasr_pm, "cor") # positive values: model unterpredicts correlation, negative vice versa. values > |.1| are woth a look
# 
# anova(sasr_pm, sasr_fm)
# 
# semPaths(object = sasr_pm, layout = "circle2", rotation = 1,
#          whatLabels="std", label.cex = 2, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE
#          , edgeLabels = sem_sig_labels(sasr_pm))
# semPaths(object = sasr_fm, layout = "circle2", rotation = 1,
#          whatLabels="std", label.cex = 2, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE
#          , edgeLabels = sem_sig_labels(sasr_fm))
# summary(sasr_pm)
# # models do not differ, prefer full mediation due to parsimony and better AIC
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # LATENT VARIABLE STUFF modifications to allow for measuring heterogeneity as one factor ####
# ## using hg =~ tri + tra_sd + sea_m + tra_m does not work, model does not converge. maybe because all variables are not significant?
# mod_base5_plush <- "
#           # latent variable
# #          hg1 =~ tri + tra_sd + sea_m + tra_m    # all negative influence
#            hg2 =~ soil + area + pre_sd   #  have positive influence on SR, all rather unsuitable for mrd: spatial heterogeneity
#           # regressions
#           sr_trans ~ hg2 
#           #+ soil
#           + sub_trop_mbf 
#           + sub_trop_dbf
#           #+ area
#           #+ pre_sd 
#           + mont_gs 
#           + pre_m 
#           + mrd 
#           + tri
#           + tra_sd
#           + elev_m
#           + sea_m
#           + tra_m
#           
#           mrd ~ pre_m + 
#           #+ hg
#           + sea_m 
#           + tri
#           + sea_sd 
#           + soil 
#           + tra_sd 
#           + pre_sd
#           + elev_m 
#           + tra_m 
#           + sub_trop_mbf
#           + pet_m 
#           + mont_gs
# "
# fsem_base5_plush <- sem(mod_base5_plush, data = dat_no.na) 
# #summary(fsem_base5_plush, standardized = TRUE, fit.measures=TRUE)
# table2<-parameterEstimates(fsem_base5_plush,standardized=TRUE)[!is.na(parameterEstimates(fsem_base5_plush)$pvalue),]
# sig <- rep(" ",nrow(table2))
# sig[table2$pvalue<=0.05] <- "*"
# sig[table2$pvalue<=0.01] <- "**"
# sig[table2$pvalue<=0.001] <- "***"
# b2 <- paste0(round(table2$std.all,2), sig)
# semPaths(object = fsem_base5_plush, layout = "tree", rotation = 1,
#          whatLabels="std", label.cex = 1.5, edge.label.cex = 1, what = "std", 
#          exoCov = FALSE, exoVar = FALSE,
#          nCharNodes = 0, fade=FALSE,
#          edgeLabels = b2)
# fitmeasures(fsem_base5_plush, c("aic", "cfi", "rmsea", "srmr"))
# modificationindices(fsem_base5, sort = TRUE)
# # worse model fit
# 
# 
# # adding heterogeneity as a latent variable
# mod3 <- "# latent variable
#          heterogeneity =~ sea_m  + tri + area + soil
#         # tri ~~ area
#         # regressions
#           sr_trans ~ heterogeneity + mrd + area + sea_m
#           mrd ~ heterogeneity
#         # mods
# #        sea_m  ~~ sr_trans
# "
# 
# set.seed(1234)
# fsem3 <- sem(mod3, data = dat_no.na)#, se = "bootstrap", bootstrap = 100) 
# #summary(fsem3, standardized = TRUE, fit.measures=TRUE)
# fitmeasures(fsem3, c("aic", "cfi", "rmsea", "srmr"))
# # this model is really good!
# semPaths(object = fsem3, layout = "tree", rotation = 1,
#          whatLabels="std", label.cex = 1.5, edge.label.cex = 1.5, what="std", exoCov = FALSE, exoVar = FALSE)
# # negative connection between sr + mrd
# # heterogeneity is mostly defined soil, area, seasonality
# 
# 
# mod4 <- "# latent variable
#          heterogeneity =~ sea_m  + tri + area + soil
#         # tri ~~ area
#         # regressions
#           sr_trans ~ heterogeneity + area + mrd + sea_m
#           mrd ~ heterogeneity
#         # mods
#         #sea_m  ~~ sr_trans
# "
# 
# set.seed(1234)
# fsem3 <- sem(mod3, data = dat_no.na)#, se = "bootstrap", bootstrap = 100) 
# summary(fsem3, standardized = TRUE, fit.measures=TRUE, rsq=TRUE)
# 
# 
# 
# 
# 
# 
# 
# 
# ### MERGE variable importances from subset models with global model
# 
# varImp(gbmFit4)
# glob <- cbind(rownames(summary(gbmFit4)), rank(summary(gbmFit4)$rel.inf, ))
# 
# 
# 
# 
# 
# 
# 
# 
# 
