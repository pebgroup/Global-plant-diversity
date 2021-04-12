# Analysis script for WCSP paper
# loaded data:
# comm: community matrix. rows represents regions, columns WCSP species IDs
# phytoregions shapefile
# mean root distance

rm(list=ls())
library(ggplot2)
library(tidyverse)
library(plyr)
library(ape)
library(rgdal)
library(rgeos)
library(phytools)
library(VIM)
library(ggsci)

theme_set(theme_bw())


# DATA PREP #############################################




## trees
phylo <- read.tree("trees/allmb_matched_added_species_Nov20.tre") # allmb_matched_added_species_3.tre

#plot.phylo(phylo, show.tip.label = FALSE, type="fan")
## occurrence matrix
#load("data/comm.RData")
sr <- readRDS("data/comm_Feb2021.rds")
## MRD
#mrd <- readRDS("data/polytomie_RD_results.rds")
mrd <- readRDS("data/mrd_feb2021.rds")
mdr <- readRDS("data/mdr_Mar2021.rds")

## shapefile
shape <- readOGR("shapefile/level3.shp")
trueCentroids = as.data.frame(gCentroid(shape,byid=TRUE))



length(phylo$tip.label)
dim(sr)
dim(mrd)


# check how many botanical countries will be removed due to incomplete climate data:
dat_no.na <- readRDS("data/sem_input_data.rds")
table(rownames(sr) %in% rownames(dat_no.na))
# ex_countries <- rownames(sr)[which(!rownames(sr) %in% rownames(dat_no.na))]
# sr2 <- sr[rownames(sr) %in% rownames(dat_no.na),]
# sr2 <- na.omit(sr2)
# dim(sr)
# dim(sr2)



# Calc species richness
sr_df <- data.frame(species_richness = rowSums(sr), region =row.names(sr))
# to compare to old SR estimates that accidentally included ferns+mosses, read old community matrix comm_Nov2020.rds
# sr_df$species_richness_phylo_data_only <- rowSums(sr_sub)
# plot(sr_df$species_richness, sr_df$species_richness_phylo_data_only)
# abline(1,1)
# cor(sr_df$species_richness, sr_df$species_richness_phylo_data_only)
# sr_df$diff <- sr_df$species_richness - sr_df$species_richness_phylo_data_only
# hist(sr_df$diff, breaks=40)

rm(phylo)
#saveRDS(sr_df, "number_of_species_")


# Get standardized MRD z score
## subtract observed and expected devide by standard error (x-mu/sigma)
##(Hawkins, B. A., Diniz-Filho, J. A. F., & Soeller, S. A. (2005). Water links the historical and contemporary components of the Australian bird diversity gradient. Journal of Biogeography, 32(6), 1035–1042. https://doi.org/10.1111/j.1365-2699.2004.01238.x)
## --> Z-scores higher than 1.96 indicate that there is a 95% chance that the MRD in the cell is higher than would be expected if species were random global sample. This is supposed to account for differences in species richness between the plots. 

# 
# 
# # plot random distribution with observed value
# mrd.p <- mrd %>% 
#   as.data.frame() %>%
#   add_column(region = row.names(mrd)) %>%
#   pivot_longer(
#     cols = starts_with("rnd")
#   )
# 
# # plot in chunks of 10x10
# ggplot(mrd.p[mrd.p$region %in% unique(mrd.p$region)[1:50],], aes(x=value)) + 
#   geom_histogram() + 
#   geom_vline(aes(xintercept = obs), col="red") + 
#   facet_wrap(~region, scales = "free")
# 
# #  z <- (mrd[,1] - apply(mrd[,c(2:100)], 1, mean)) / (apply(mrd[,c(2:100)], 1, sd)/sqrt(99))
   z <- (mrd[,1] - apply(mrd[,c(2:100)], 1, mean)) / (apply(mrd[,c(2:100)], 1, sd))

#mrd <- cbind(mrd, as.numeric(z))


# Add SR and MRD and centroids to shapefile
shape@data$sr <- sr_df$species_richness
#shape@data$sr_new <- sr_df$species_richness_phylo_data_only
shape@data$mrd <- mrd[,1]
shape@data$mdr <- mdr[,1]
shape@data$mrd_z_score <- z
shape@data$mrd_z_score_binary <- shape@data$mrd_z_score
shape@data$mrd_z_score_binary[which(shape@data$mrd_z_score< -1.97)] <- "-"
shape@data$mrd_z_score_binary[which(shape@data$mrd_z_score>1.97)] <- "+"
shape@data$mrd_z_score_binary[which(shape@data$mrd_z_score>-1.97 & shape@data$mrd_z_score<1.97)] <- "zero"
shape@data$lng <- trueCentroids[,1]
shape@data$lat <- trueCentroids[,2]

MRD.sd <- readRDS("results/MRD_standard_deviation.rds")
MDR.sd  <- readRDS("results/MDR_standard_deviation.rds")
shape@data$MRD.sd <- MRD.sd
shape@data$MDR.sd <- MDR.sd

range(shape@data$lat)
breaks <- seq(-85,85,10)
# specify interval/bin labels
tags <- seq(-80, 80, 10)
# bucketing values into bins
lat_bins <- cut(shape@data$lat, 
                breaks=breaks, 
                include.lowest=TRUE, 
                right=FALSE, 
                labels=tags)
shape@data$lat_bin <- lat_bins

cor.test(shape@data$sr, shape@data$mrd, method="p")
cor.test(shape@data$sr, shape@data$mrd_z_score, method="p")
cor.test(shape@data$sr, shape@data$mdr, method="p")

rm(sr)


#### Add environmental variables to shapefile #####################################


soil <- readRDS("data/soil.rds")
shape@data$soil <- soil$number_soils
shape@data$soil_simp <- soil$simp_soil
shape@data$soil_even <- soil$even_soil

climate <- readRDS("data/climate.rds")
shape@data <- cbind(shape@data, climate)

topography <- readRDS("data/topography.rds")
shape@data$elev_mean <- topography$elev_mean
shape@data$elev_sd <- topography$elev_sd
shape@data$elev_n <- topography$elev_n
shape@data$tri <- topography$ruggedness_mean

# add area
library(geosphere)
shape$area <- areaPolygon(shape)


# add biomes ##########################################################
biomes <- readRDS("data/biomes_olson.rds")
biome_names <- c("(sub)tropical moist broadleaf forest",
                 "(sub)tropical dry broadleaf forest",
                 "(sub)tropical coniferious forest",
                 "temperate broadleaf/mixed forest",
                 "temperate coniferous forest",
                 "boreal forests/taiga",
                 "(sub)tropical grasslands, savannas, shrublands",
                 "temperate grasslands, savannas, shrublands",
                 "flooded grasslands and savannas",
                 "montane grasslands and shrublands",
                 "tundra",
                 "mediterranean forests, woodlands, scrub",
                 "deserts and xeric shrublands",
                 "mangroves")
names(biomes)[grepl("^X", names(biomes))] <- biome_names
shape@data <- cbind(shape@data, biomes[,-1])


# add paleocim ##########################################################
paleoclim <-readRDS("data/paleoclim.rds")
head(paleoclim)
all(shape@data$LEVEL_3_CO==rownames(paleoclim)) # identical?
shape@data <- cbind(shape@data, paleoclim)
names(shape@data)



# # add human footprint index
# hfp <- readRDS("data/hfp.rds")
# shape@data$hfp90 <- hfp$hfp90
# shape@data$hfp_mean <- hfp$hfp_mean

# Transform to data frame for ggplot
library(sf)
shp <- st_as_sf(shape)



# mean root distance between countries
hist(shp$mrd, breaks=20, xlab="mean root distance per country")





# missing data ##########################################################################


dat <- st_drop_geometry(shp)
names(dat)
dat <- dat[,c(12:ncol(dat))]

aggr_plot <- aggr(dat[,c(12:ncol(dat))], col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE,
                  labels=names(data), cex.axis=.9, gap=3, ylab=c("Histogram of missing data","Pattern"))




# 0.84 complete cases for all variables (0.86 without hfp)

# check without the SDs
aggr_plot <- aggr(dat[,-grep("_sd", names(dat))], col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE,
                  labels=names(data), cex.axis=.9, gap=3, ylab=c("Histogram of missing data","Pattern"))

# 0.91 complete cases for all variables except SDs
# 0.86 complete cases for all variables except SDs AND HFP


plot(shp$mrd, shp$mdr)
cor.test(shp$mrd, shp$mdr, method="s")
cor.test(shp$sr, shp$mrd)
cor.test(shp$sr, shp$mdr)

plot(shp$sr, shp$mrd)
plot(shp$sr, shp$mdr)

# SR MAP ##########################################################################
ggplot(shp) + 
  geom_sf(aes(fill = sr), lwd=0.1) + 
  scale_fill_viridis_c(option = "plasma", trans = "sqrt")+ #, 
#  scale_fill_gradient("Species richness", guide="colorbar", low="white", high="darkred",
#                      limits=c(min(sr_df$species_richness),max(sr_df$species_richness)))+
  guides(fill = guide_colourbar(barwidth = 20, direction="horizontal"))+ # stretch that colorbar
  theme_void()+
  theme(legend.position = "bottom")
ggsave(filename=paste0("results/sr_map_", gsub("-", "_", Sys.Date()), ".png"), dpi=600, width=10, height=7)

temp <- readRDS("data/data_for_SEM.rds")
temp2 <- shp[shp$LEVEL_3_CO %in% rownames(temp),]

#get buffer around tiny countries
# centroids <- st_centroid(temp2)
# pt <- centroids[which(temp2$area<1200000000),] 
# pt_buffer <- st_buffer(pt,2) 
# # thicker lines
thicc_lines <-temp2[which(temp2$area<1200000000),]
# temp2$thicc_lines <- 0.01
# temp2$thicc_lines[which(temp2$area<1200000000)] <- 3
lcol <- min(thicc_lines$sr)/max(temp2$sr)
ucol <- max(thicc_lines$sr)/max(temp2$sr)

ggplot(temp2) + 
  geom_sf(aes(fill=sr),lwd=0.1) + 
  geom_sf(data=thicc_lines, lwd=1.5, aes(col=sr), show.legend=F)+
  scale_colour_viridis_c("SR", option = "plasma", trans = "sqrt", 
                         begin = lcol, end = sqrt(ucol))+
  scale_fill_viridis_c("SR", option = "plasma", trans = "sqrt")+ #, 
  #  scale_fill_gradient("Species richness", guide="colorbar", low="white", high="darkred",
  #                      limits=c(min(sr_df$species_richness),max(sr_df$species_richness)))+
  #guides(fill = guide_colourbar(barwidth = 20, direction="vertical"))+ # stretch that colorbar
  theme_void()+
  theme(legend.position = c(0.21, 0.3), 
        legend.background = element_blank(),
        legend.key = element_blank())
ggsave(filename=paste0("results/sr_map_red_", gsub("-", "_", Sys.Date()), ".png"), dpi=600, width=10, height=7)

# ggplot(shp) + 
#   geom_sf(aes(fill = sr_new), lwd=0.1) + 
#   scale_fill_viridis_c(option = "plasma", trans = "sqrt")+ #, 
#   #  scale_fill_gradient("Species richness", guide="colorbar", low="white", high="darkred",
#   #                      limits=c(min(sr_df$species_richness),max(sr_df$species_richness)))+
#   guides(fill = guide_colourbar(barwidth = 20, direction="horizontal"))+ # stretch that colorbar
#   theme_void()+
#   theme(legend.position = "bottom")
# ggsave(filename=paste0("results/sr_new_map_", gsub("-", "_", Sys.Date()), ".png"), dpi=600, width=10, height=7)


# MRD MAP ##########################################################################
lcol <- (min(thicc_lines$mrd)-min(temp2$mrd))/(max(temp2$mrd)-min(temp2$mrd))
ucol <- (max(thicc_lines$mrd)-min(temp2$mrd))/(max(temp2$mrd)-min(temp2$mrd))
ggplot(temp2) + 
  geom_sf(aes(fill = mrd), lwd=0.1) + 
  geom_sf(data=thicc_lines, lwd=1.5, aes(col=mrd), show.legend=F)+
  scale_colour_viridis_c("SR", option = "plasma", 
                         begin = lcol, end = ucol)+
  scale_fill_viridis_c(option = "plasma")+ #, 
  theme_void()+
  theme(legend.position = c(0.21, 0.3), 
        legend.background = element_blank(),
        legend.key = element_blank())
ggsave(filename=paste0("results/mrd_map_red_", gsub("-", "_", Sys.Date()), ".png"), dpi=600, width=10, height=7)

# ggplot(shp) + 
#   geom_sf(aes(fill = mrd), lwd=0.1) + 
#   scale_fill_viridis_c(option = "plasma")+ #, 
#   guides(fill = guide_colourbar(barwidth = 20, direction="horizontal"))+ # stretch that colorbar
#   theme_void()+
#   theme(legend.position = "bottom")
# ggsave(filename=paste0("results/mrd_map_", gsub("-", "_", Sys.Date()), ".png"), dpi=600, width=10, height=7)
# 
# 
# ggplot(shp) + 
#   geom_sf(aes(fill = mrd_z_score), lwd=0.1) + 
#   scale_fill_viridis_c(option = "plasma")+ #, 
#   guides(fill = guide_colourbar(barwidth = 20, direction="horizontal"))+ # stretch that colorbar
#   theme_void()+
#   theme(legend.position = "bottom")
# ggsave(filename=paste0("results/mrd_z_score_map", gsub("-", "_", Sys.Date()), ".png"), dpi=600, width=10, height=7)
# 

# EQUAL SPLITS MAP ##########################################################################
ggplot(shp) + 
  geom_sf(aes(fill = mdr), lwd=0.1) + 
  scale_fill_viridis_c(option = "plasma")+ #, 
  guides(fill = guide_colourbar(barwidth = 20, direction="horizontal"))+ # stretch that colorbar
  theme_void()+
  theme(legend.position = "bottom")
ggsave(filename=paste0("results/mdr_map_", gsub("-", "_", Sys.Date()), ".png"), dpi=600, width=10, height=7)

library(gridExtra)
grid.arrange(ncol=1,
ggplot(shp) + 
  geom_sf(aes(fill = mrd), lwd=0.1) + 
  scale_fill_viridis_c(option = "plasma")+ #, 
#  guides(fill = guide_colourbar(barwidth = 20, direction="horizontal"))+ # stretch that colorbar
  theme_void()
,
ggplot(shp) + 
  geom_sf(aes(fill = mdr), lwd=0.1) + 
  scale_fill_viridis_c("1/ES", option = "plasma", trans = "log")+ #, 
 # guides(fill = guide_colourbar(barwidth = 20, direction="horizontal"))+ # stretch that colorbar
  theme_void()
)



## uncertainty maps ###########################################################################3
grid.arrange(ncol=1,
             ggplot(shp) + 
               geom_sf(aes(fill = MRD.sd), lwd=0.1) + 
               scale_fill_viridis_c("MRD sd", option = "plasma")+ #, 
               #  guides(fill = guide_colourbar(barwidth = 20, direction="horizontal"))+ # stretch that colorbar
               theme_void()
             ,
             ggplot(shp) + 
               geom_sf(aes(fill = MDR.sd), lwd=0.1) + 
               scale_fill_viridis_c("Equal splits SD", option = "plasma")+ #, 
               # guides(fill = guide_colourbar(barwidth = 20, direction="horizontal"))+ # stretch that colorbar
               theme_void()
)
cor.test(shp$MRD.sd, shp$MDR.sd)

# SCALED COMPARISON MAPS ########################################################
# a map that shows....? deviations from the average SR? what about diffs between SR+mrd...
# range01 <- function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))}
# #ztrans <- function(x){(x - mean(x, na.rm = T)) / sd(x, na.rm = T)}
# shp$sr_scaled <- range01(shp$sr)
# #shp$sr_ztrans <- ztrans(shp$sr)
# #shp$mrd_ztrans <- ztrans(shp$mrd)
# shp$mrd_scaled <- range01(shp$mrd)
# ggplot(shp) + 
#   geom_sf(aes(fill = mrd_scaled-sr_scaled), lwd=0.1) + 
#   scale_fill_viridis_c(option = "plasma")+ #, 
#   #  scale_fill_gradient("Species richness", guide="colorbar", low="white", high="darkred",
#   #                      limits=c(min(sr_df$species_richness),max(sr_df$species_richness)))+
#   guides(fill = guide_colourbar(barwidth = 20, direction="horizontal"))+ # stretch that colorbar
#   theme_void()+
#   theme(legend.position = "bottom")
# this shows the largest discrepancies between mrd and sr (smallest values), but it does not tell if its due to low MRD or big SR.


# LDG / longitudinal gradients #########################################################################
## global lat gradient: scatterplot SR + MRD vs latitude
library(scales)
(global_scatterplot <- ggplot(shp, aes(lat, sr))+
  geom_point(alpha=0.5)+
  scale_y_continuous("species richness", trans = "sqrt")+
  scale_x_continuous("latitude centroid botanical country")+
  geom_smooth(data=shp[shp$lat<0,], method="lm", col="red")+
  geom_smooth(data=shp[shp$lat>0,], method="lm", col="blue")
)
(global_scatterplot_mrd <- ggplot(shp, aes(lat, mrd))+
    geom_point(alpha=0.5)+
    scale_y_continuous("mean root distance")+
    scale_x_continuous("latitude centroid botanical country")+    
    geom_smooth(data=shp[shp$lat<0,], method="lm", col="red")+
    geom_smooth(data=shp[shp$lat>0,], method="lm", col="blue")
)
grid.arrange(nrow=1, global_scatterplot, global_scatterplot_mrd)
#, label = unit_format(unit = "K"), limits=c(0,max(shp$sr)/1000)
#  geom_smooth(data=shp[shp$lat<0,], method="lm", col="red")+
#  geom_smooth(data=shp[shp$lat>0,], method="lm", col="blue")
#geom_text(aes(label=LEVEL_3_CO, col=factor(mrd_z_score_binary)),hjust=0, vjust=0)+
cor.test(abs(shp$lat), shp$sr, method = "p") # pearson correlation = -0.27
cor.test(shp$lat[shp$lat<0], shp$sr[shp$lat<0], method = "p")
cor.test(shp$lat[shp$lat>0], shp$sr[shp$lat>0], method = "p")
summary(lm(sr~lat, data=shp[shp$lat<0,]))
summary(lm(sr~lat, data=shp[shp$lat>0,]))
summary(lm(sr~abs(lat), data=shp))

cor.test(abs(shp$lat), shp$mrd, method = "p")
cor.test(shp$lat[shp$lat<0], shp$mrd[shp$lat<0], method = "p")
cor.test(shp$lat[shp$lat>0], shp$mrd[shp$lat>0], method = "p")
summary(lm(mrd~lat, data=shp[shp$lat<0,]))
summary(lm(mrd~lat, data=shp[shp$lat>0,]))
summary(lm(mrd~abs(lat), data=shp))
# Both hemispheres have significant LDG away from the equator

## continental gradients: scatterplot SR vs latitude for each continent (just use the ones in shapefile)
ztrans <- function(x){(x - mean(x, na.rm = T)) / sd(x, na.rm = T)}
range01 <- function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))}
shpbins <- shp
shpbins$sr_scaled <- range01(shp$sr)
shpbins$mrd_scaled <- range01(shp$mrd)
range(shpbins$mrd_scaled, na.rm=T)
range(shpbins$sr_scaled, na.rm=T)
shpbins <- st_drop_geometry(shpbins)
temp <- pivot_longer(shpbins[,c("lat","mrd_scaled", "sr_scaled", "CONTINENT")], cols = c("mrd_scaled", "sr_scaled"))

(continent_scatterplot <- ggplot(temp, aes(lat, value, col=name))+
  geom_point(alpha=0.5)+
  scale_y_continuous("scaled MRD and SR")+#,
                     #sec.axis = sec_axis(~.*1, name="mean root distance"))+ 
  scale_x_continuous("latitude centroid botanical country")+
  facet_wrap(~CONTINENT, nrow = 3, scale="free_x")+ #, scale="free"
  geom_smooth(method="loess", se=F)+
  scale_colour_startrek(labels=c("MRD", "SR"), name="")
)
ggsave("scatterplot_SR_MRD_latitude_continents.png", dpi=600, width=6, height=4)


## same plots for longitude: boxplots for x degree increments for longitude
# (global_scatterplot_lng <- ggplot(shp, aes(lng, sr))+
#     geom_point(alpha=0.5)+
#     scale_y_continuous("species richness", trans = "sqrt")+
#     scale_x_continuous("longitude centroid botanical country")#+
# )
# (global_scatterplot_mrd_lng <- ggplot(shp, aes(lng, mrd))+
#     geom_point(alpha=0.5)+
#     scale_y_continuous("mean root distance")+
#     scale_x_continuous("longitude centroid botanical country")#+
# )

## LDG stuff reduced data that has been used in SEMs ####
shpbins <- temp2
shpbins$sr_scaled <- range01(shpbins$sr)
shpbins$mrd_scaled <- range01(shpbins$mrd)
range(shpbins$mrd_scaled, na.rm=T)
range(shpbins$sr_scaled, na.rm=T)
shpbins <- st_drop_geometry(shpbins)
temp <- pivot_longer(shpbins[,c("lat","mrd_scaled", "sr_scaled", "CONTINENT")], cols = c("mrd_scaled", "sr_scaled"))

(continent_scatterplot <- ggplot(temp, aes(lat, value, col=name))+
    geom_point(alpha=0.5)+
    scale_y_continuous("scaled MRD and SR")+#,
    #sec.axis = sec_axis(~.*1, name="mean root distance"))+ 
    scale_x_continuous("latitude centroid botanical country")+
    facet_wrap(~CONTINENT, nrow = 3, scale="free_x")+ #, scale="free"
    geom_smooth(method="loess", se=F)+
    scale_colour_startrek(labels=c("MRD", "SR"), name="")
)
ggsave("scatterplot_SR_MRD_latitude_continents_red.png", dpi=600, width=6, height=4)

(global_scatterplot_red <- ggplot(temp2, aes(lat, sr))+
    geom_point(alpha=0.5)+
    scale_y_continuous("species richness", trans = "sqrt")+
    scale_x_continuous("latitude centroid botanical country")+
    geom_smooth(data=shp[shp$lat<0,], method="lm", col="red")+
    geom_smooth(data=shp[shp$lat>0,], method="lm", col="blue")
)
(global_scatterplot_mrd_red <- ggplot(temp2, aes(lat, mrd))+
    geom_point(alpha=0.5)+
    scale_y_continuous("mean root distance")+
    scale_x_continuous("latitude centroid botanical country")+    
    geom_smooth(data=shp[shp$lat<0,], method="lm", col="red")+
    geom_smooth(data=shp[shp$lat>0,], method="lm", col="blue")
)
grid.arrange(nrow=1, global_scatterplot_red, global_scatterplot_mrd_red)
#, label = unit_format(unit = "K"), limits=c(0,max(shp$sr)/1000)
#  geom_smooth(data=shp[shp$lat<0,], method="lm", col="red")+
#  geom_smooth(data=shp[shp$lat>0,], method="lm", col="blue")
#geom_text(aes(label=LEVEL_3_CO, col=factor(mrd_z_score_binary)),hjust=0, vjust=0)+
cor.test(abs(temp2$lat), temp2$sr, method = "p") # pearson correlation = -0.38
cor.test(temp2$lat[temp2$lat<0], temp2$sr[temp2$lat<0], method = "p")
cor.test(temp2$lat[temp2$lat>0], temp2$sr[temp2$lat>0], method = "p")
# summary(lm(sr~lat, data=shp[shp$lat<0,]))
# summary(lm(sr~lat, data=shp[shp$lat>0,]))
# summary(lm(sr~abs(lat), data=shp))

cor.test(abs(temp2$lat), temp2$mrd, method = "p")
cor.test(temp2$lat[temp2$lat<0], temp2$mrd[temp2$lat<0], method = "p")
cor.test(temp2$lat[temp2$lat>0], temp2$mrd[temp2$lat>0], method = "p")
summary(lm(mrd~lat, data=temp2[temp2$lat<0,]))
summary(lm(mrd~lat, data=temp2[temp2$lat>0,]))
summary(lm(mrd~abs(lat), data=temp2))

## continental gradients: scatterplot SR vs latitude for each continent (just use the ones in shapefile)
ztrans <- function(x){(x - mean(x, na.rm = T)) / sd(x, na.rm = T)}
range01 <- function(x){(x-min(x, na.rm = TRUE))/(max(x, na.rm = TRUE)-min(x, na.rm = TRUE))}
shpbins <- shp
shpbins$sr_scaled <- range01(shp$sr)
shpbins$mrd_scaled <- range01(shp$mrd)
range(shpbins$mrd_scaled, na.rm=T)
range(shpbins$sr_scaled, na.rm=T)
shpbins <- st_drop_geometry(shpbins)
temp <- pivot_longer(shpbins[,c("lat","mrd_scaled", "sr_scaled", "CONTINENT")], cols = c("mrd_scaled", "sr_scaled"))

(continent_scatterplot <- ggplot(temp, aes(lat, value, col=name))+
    geom_point(alpha=0.5)+
    scale_y_continuous("scaled MRD and SR")+#,
    #sec.axis = sec_axis(~.*1, name="mean root distance"))+ 
    scale_x_continuous("latitude centroid botanical country")+
    facet_wrap(~CONTINENT, nrow = 3, scale="free_x")+ #, scale="free"
    geom_smooth(method="loess", se=F)+
    scale_colour_startrek(labels=c("MRD", "SR"), name="")
)
ggsave("scatterplot_SR_MRD_latitude_continents.png", dpi=600, width=6, height=4)



# 
# # bins
# breaks <- seq(-180,180,20)
# # specify interval/bin labels
# tags <- seq(-170, 170, 20)
# hist(shape$lng)
# lng_bins <- cut(shape@data$lng, 
#                 breaks=breaks, 
#                 include.lowest=TRUE, 
#                 right=FALSE, 
#                 labels=tags)
# shape@data$lng_bin <- lng_bins
# # plot 10 degree bins
# ggplot(shp, aes(lng_bin, sr)) + 
#   geom_boxplot(varwidth = TRUE, fill=c(rep("blue",5),rep("red",5),rep("blue",6)))
# 
# shape@data$lat_bin <- lat_bins
# # plot 10 degree bins
# ggplot(shp, aes(lat_bin, sr)) + 
#   geom_boxplot(varwidth = TRUE, fill=c(rep("blue",5),rep("red",5),rep("blue",6))) +
#   coord_flip()
# ggsave(filename=paste0("results/lat_bands_SR", gsub("-", "_", Sys.Date()), ".png"), dpi=600, width=6, height=4)
# 
# kruskal.test(shp$sr, shp$lat_bin)
# pairwise.wilcox.test(shp$sr, shp$lat_bin, p.adjust.method = "fdr")
# 
# 
# ggplot(shp, aes(lat_bin, mrd_z_score)) + 
#   geom_boxplot(varwidth = TRUE, fill=c(rep("blue",5),rep("red",5),rep("blue",6))) +
# #  scale_y_continuous(limits = c(-20,30))+
#   coord_flip()
# ggsave(filename=paste0("results/lat_bands_MRD_EcolLetFig3", gsub("-", "_", Sys.Date()), ".png"), dpi=600, width=6, height=4)
# 
# 
# table(shp$lat_bin)
# which(shp$lat < -25 | shp$lat > 25)
# shp$trops <- "tropical"
# shp$trops[which(shp$lat < -25 | shp$lat > 25)] <- "non-tropical"
# ggplot(shp, aes(trops, mrd)) + 
#   geom_boxplot(fill=c("blue", "red"))
# ggsave(filename=paste0("results/MRD_EcolLetFig1", gsub("-", "_", Sys.Date()), ".png"), dpi=600, width=3, height=3)
# kruskal.test(shp$mrd_z_score, shp$trops)

# # rank ordered distribution of mrd_z_score (figure1 ecol)
# qqnorm(y=sort(shp$mrd_z_score[shp$trops=="tropical"]),col="red", datax = TRUE)
# qqnorm(y=sort(shp$mrd_z_score[shp$trops=="non-tropical"]),col="blue", datax = TRUE, add=TRUE)
# #lines(sort(shp$mrd_z_score[shp$trops=="non-tropical"]),col="blue", type="l")
# #qqline()

# ggplot(shp, aes(sr, lat))+
#   geom_point(aes(col=mrd_z_score_binary))+
#   scale_x_log10()



# BIOME MAP ###############################################

ggplot(shp) + 
  geom_sf(aes(fill = factor(predominant_biome)), lwd=0.1) + 
  scale_fill_viridis_d()+
  theme_void()#+
#  theme(legend.position = "bottom")









saveRDS(shp, "data/shp_object_fin_analysis.RDS")



















#### Scaling MRD effects ###########################################


shp <- readRDS("data/shp_object_fin_analysis.RDS")
library(sf)
dat <- st_drop_geometry(shp)

# Remove NAs
# dat.cl <- na.omit(dat)

# # Is scaling MRD really necessary?
# # check residuals for area pattern
# lm.mrd <- lm(dat.cl$mrd~dat.cl$mrd_z_score)
# plot(lm.mrd)
# dat.cl$mrd_residuals <- lm.mrd$residuals
# 
# 
# library(gridExtra)
# grid.arrange(ncol=2,
#   ggplot(dat.cl, aes(sr, mrd, label=LEVEL_3_CO))+
#     geom_text(aes(label=LEVEL_3_CO, hjust=0, vjust=0), size=3)+
#     scale_x_log10()+
#     geom_smooth(method="lm"),
#   ggplot(dat.cl, aes(sr, mrd_z_score, label=LEVEL_3_CO))+
#     geom_text(aes(label=LEVEL_3_CO),hjust=0, vjust=0, size=3)+
#     scale_x_log10()+
#     geom_smooth(method="lm"),
#   ggplot(dat.cl, aes(mrd, mrd_z_score, label=LEVEL_3_CO))+
#     geom_text(aes(label=LEVEL_3_CO, col=log(area)),hjust=0, vjust=0, size=3)+
#     geom_smooth(method="lm"),
#   ggplot(dat.cl, aes(area, mrd_residuals, label=LEVEL_3_CO))+
#     geom_text(aes(label=LEVEL_3_CO),hjust=0, vjust=0, size=3)+
#     scale_x_log10()+
#     geom_smooth(method="lm")
# )
# ggplot(dat.cl, aes(mrd, mrd_z_score, label=LEVEL_3_CO))+
#   geom_text(aes(label=LEVEL_3_CO, col=log(sr)),hjust=0, vjust=0, size=3)+
#   geom_smooth(method="lm")
# ggplot(dat.cl, aes(sr, mrd_residuals, label=LEVEL_3_CO))+
#   geom_text(aes(label=LEVEL_3_CO),hjust=0, vjust=0, size=3)+
#   scale_x_log10()+
#   geom_smooth(method="lm")

# no
#plot(dat.cl$mrd, dat.cl$mrd_z_score, pch=factor(dat.cl$LEVEL_3_CO))
#abline(lm(dat.cl$mrd_z_score~dat.cl$mrd))
# countries affected by scaling:






# CORRELATIONS and MULTICOLINERITY ########################################################


rownames(dat) <- dat$LEVEL_3_CO
dat <- dat[,c(grep("sr", names(dat)):ncol(dat))]
dat_means <- dat[,-grep("lng|lat_bin|z_score|trops|binary|_n|_sd|_even|_simp", names(dat))]
dat_sds <- dat[,-grep("lng|lat_bin|z_score|trops|binary|_n|_mean|tri|predominant", names(dat))]


aggr_plot <- aggr(dat_means, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE,   # [,c(12:ncol(dat))]
                   labels=names(data), cex.axis=.9, gap=3, ylab=c("Histogram of missing data","Pattern"))
aggr_plot <- aggr(dat_sds, col=c('navyblue','red'), numbers=TRUE, sortVars=TRUE,   # [,c(12:ncol(dat))]
                  labels=names(data), cex.axis=.9, gap=3, ylab=c("Histogram of missing data","Pattern"))


dat_means <- na.omit(dat_means)
dat_sds <- na.omit(dat_sds)


library(ggcorrplot)
cor.dat.means<- cor(dat_means, method="s")
p.dat.means <- cor_pmat(dat_means, method="s")
ggcorrplot(cor.dat.means[,c(1:12)], lab=TRUE, p.mat = p.dat.means[,c(1:12)], insig = "blank")

cor.dat.sds<- cor(dat_sds, method="s")
p.dat.sds <- cor_pmat(dat_sds, method="s")
ggcorrplot(cor.dat.sds[,c(1:13)], lab=TRUE, p.mat = p.dat.sds[,c(1:13)], insig = "blank")


plot(log1p(dat$area), log1p(dat$sr))
plot(log1p(dat$area), log1p(dat$mrd))
plot(dat$lat, dat$olson_percent)

hist(log(dat$area[which(is.na(dat$pet_mean))]))

ggplot(dat, aes(area, sr, group=is.na(pet_mean)))+
  geom_point(aes(col=is.na(dat$pet_mean)))+
  scale_x_log10()+
  scale_y_log10()

ggplot(dat, aes(area, mrd_z_score, group=is.na(pet_mean)))+
  geom_point(aes(col=is.na(dat$pet_mean)))+
  scale_x_log10()+
  scale_y_continuous(trans = "log1p")


# biome correlations
dat.biome <- dat_means[grep("sr|mrd|forest|land|tundra|mangrove|number_biomes|predominant", names(dat_means))]
cor.biome <- cor(dat.biome, method="s")
p.biome <- cor_pmat(dat.biome, method="s")
ggcorrplot(cor.biome[,c(1:2)], lab=TRUE, p.mat = p.biome[,c(1:2)], insig = "blank")

# without excluding missing data countries:
biomes$sr <- shp$sr
biomes$mrd <- shp$mrd
cor.test(biomes$sr, biomes$`(sub)tropical moist broadleaf forest`, method="s")
cor.test(biomes$mrd, biomes$`(sub)tropical moist broadleaf forest`, method="s")
# only complete countries for means
cor.test(dat_means$sr, dat_means$`(sub)tropical moist broadleaf forest`, method="s")
cor.test(dat_means$mrd, dat_means$`(sub)tropical moist broadleaf forest`, method="s")
# only complete countries for SDs
cor.test(dat_sds$sr, dat_sds$`(sub)tropical moist broadleaf forest`, method="s")
cor.test(dat_sds$mrd, dat_sds$`(sub)tropical moist broadleaf forest`, method="s")

ggplot(dat.biome, aes(y=sr, x=predominant_biome, group=predominant_biome))+
  geom_boxplot(varwidth = TRUE)+
  scale_x_continuous(breaks=c(1:14), labels=biome_names)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
pairwise.wilcox.test(dat.biome$sr, dat.biome$predominant_biome, p.adjust.method = "fdr")

ggplot(dat.biome, aes(y=mrd, x=predominant_biome, group=predominant_biome))+
  geom_boxplot(varwidth = TRUE)+
  scale_x_continuous(breaks=c(1:14), labels=biome_names)+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
pairwise.wilcox.test(dat.biome$mrd, dat.biome$predominant_biome, p.adjust.method = "fdr")



## most common biome combinations
#get biome ranks per country
biomes <- dat.biome
#biomes <- readRDS("data/biomes_olson.rds")
library(tidyr)
long_biomes <- biomes %>% pivot_longer(cols=X1:X14)
order(s)
ggplot(long_biomes[1:100,], aes(x=country, y=value, fill=name))+
  geom_bar(stat="identity", position="stack")

boxplot(dat.biome$number_biomes~dat.biome$predominant_biome)

count <- biomes[,c(2:15)]!=0
count <- apply(count, 1, which)
count <- sapply(count, names)
combs <- lapply(count, paste, collapse=" ")
combs <- unlist(combs)
as.data.frame(sort(table(combs)))

# which countries have the most similar composition?
bin <- biomes[,c(2:15)]
rownames(bin) <- biomes$country
#bin[bin==TRUE] <- 1
dbin <- dist(bin)
klust <- hclust(dbin)
colors()

plot(klust, cex=0.6, label=colors(distinct=TRUE)[biomes$predominant_biome])
hcd <- as.dendrogram(klust)
library("phytools")
plot(as.phylo(klust), cex = 0.6, tip.color=colors()[biomes$predominant_biome])

# Cut tree into groups
sub_grp <- cutree(klust, k = 2)
plot(as.phylo(klust), cex = 0.6, tip.color=c("blue", "red")[sub_grp])
sub_grp <- cutree(klust, k = 3)
plot(as.phylo(klust), cex = 0.6, tip.color=c("blue", "red", "green")[sub_grp])
sub_grp <- cutree(klust, k = 4)
plot(as.phylo(klust), cex = 0.6, tip.color=c("blue", "red", "green", "yellow")[sub_grp])
sub_grp <- cutree(klust, k = 5)
plot(as.phylo(klust), cex = 0.6, tip.color=c("blue", "red", "green", "yellow", "grey")[sub_grp])

# try PCA for variable reduction
pca.bin <- prcomp(bin, center = TRUE, scale. = TRUE)
pca.bin2 <- prcomp(bin[,paste0("X", c(1,2,3,6,7,10,11,14))], center = TRUE, scale. = TRUE)
summary(pca.bin)
#library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)
ggbiplot(pca.bin)+ # , groups = biomes$sr
  geom_point(size = scale(biomes$sr)+2, alpha=0.5)

# plot(hcd, type = "rectangle", ylab = "Height")
# # Define nodePar
# nodePar <- list(lab.cex = 0.6, pch = c(NA, 19), 
#                 cex = 0.7, col = colors(distinct=TRUE)[biomes$predominant_biome])
# plot(hcd, ylab = "Height", nodePar = nodePar)


# whats the most common biome (the biome that takes up most percentages)
ggplot(long_biomes[long_biomes$value!=0,], aes(x=value))+
  geom_histogram()+
  geom_vline(data = ddply(long_biomes, "name", summarize, wavg = mean(value)), aes(xintercept=wavg))+
  facet_wrap(~name, ncol=1, strip.position = "right")+
  theme(strip.background = element_blank(),
        panel.spacing.y = unit(0,"cm"))

# is the percentage the habitat is represented key to SR?
plot(long_biomes$sr[long_biomes$name=="X1"]~long_biomes$value[long_biomes$name=="X1"])

long_biomes$predominant <- long_biomes$predominant_biome==as.numeric(gsub("X", "", long_biomes$name))

ggplot(long_biomes[long_biomes$value!=0,], aes(x=value, y=sr))+
  geom_point(aes(col=predominant), alpha=.5)+
  geom_smooth(method="lm", se = FALSE)+
  facet_wrap(~name)

# this does not work, long_biomes is not a propper representation of dat.biome WHY THE FUCK
# ok got it. becuase its based on complete data....:]
cor.test(long_biomes$value[long_biomes$name=="X1"],
        long_biomes$sr[long_biomes$name=="X1"], ,method = "s")
# cor.test(long_biomes$value[long_biomes$name=="X1"],
#          long_biomes$sr[long_biomes$name=="X1"], method="s")
# cor.test(long_biomes$value[long_biomes$name=="X10" &long_biomes$value!=0],
#          long_biomes$sr[long_biomes$name=="X10" &long_biomes$value!=0])
# cor.test(long_biomes$value[long_biomes$name=="X13" &long_biomes$value!=0],
#          long_biomes$sr[long_biomes$name=="X13" &long_biomes$value!=0]).
# SR + biome percent negative correlated in X6+11


# X1 SR versus percent
plot(biomes$X1,  biomes$sr)
cor.test(biomes$X1,  biomes$sr, method="s")



ggplot(long_biomes, aes(number_biomes))+
  geom_bar()+
#  geom_smooth(method="lm", se = FALSE)+
  facet_wrap(~name)



# # human footprint:
# plot(dat.cl$hfp90, dat.cl$sr)
# plot(dat.cl$hfp90, dat.cl$mrd)
# ggplot(dat.cl, aes(y=hfp90, x=predominant_biome, group=predominant_biome))+
#   geom_boxplot()+
#   scale_x_continuous(breaks=c(1:14), labels=paste(biome_names, c(1:14)))+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
# 
# # remove hfp
# rem <-  which(grepl("hfp", names(dat.cl)))
# dat.cl <- dat.cl[,-rem]

