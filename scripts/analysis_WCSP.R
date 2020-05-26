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

theme_set(theme_bw())

# load data
## trees
#phylo <- read.tree("trees/allmb_matched_added_species.tre")
## occurrence matrix
#load("data/comm.RData")
sr <- readRDS("data/comm.rds")
## MRD
mrd <- readRDS("data/mrd.rds")
## shapefile
shape <- readOGR("shapefile/level3.shp")



# Calc species richness
sr_df <- data.frame(species_richness = rowSums(sr), region =row.names(sr))


# Get standardized MRD z score
## subtract observed and expected devide by standard error (x-mu/sigma)
##(Hawkins, B. A., Diniz-Filho, J. A. F., & Soeller, S. A. (2005). Water links the historical and contemporary components of the Australian bird diversity gradient. Journal of Biogeography, 32(6), 1035–1042. https://doi.org/10.1111/j.1365-2699.2004.01238.x)
## --> Z-scores higher than 1.96 indicate that there is a 95% chance that the MRD in the cell is higher than would be expected if species were random global sample. This is supposed to account for differences in species richness between the plots. 

# plot random distribution with observed value
mrd.p <- mrd %>% 
  as.data.frame() %>%
  add_column(region = row.names(mrd)) %>%
  pivot_longer(
    cols = starts_with("rnd")
  )

# plot in chunks of 10x10
ggplot(mrd.p[mrd.p$region %in% unique(mrd.p$region)[1:50],], aes(x=value)) + 
  geom_histogram() + 
  geom_vline(aes(xintercept = obs), col="red") + 
  facet_wrap(~region, scales = "free")

#  z <- (mrd[,1] - apply(mrd[,c(2:100)], 1, mean)) / (apply(mrd[,c(2:100)], 1, sd)/sqrt(99))
  z <- (mrd[,1] - apply(mrd[,c(2:100)], 1, mean)) / (apply(mrd[,c(2:100)], 1, sd))

#mrd <- cbind(mrd, as.numeric(z))


# Add SR and MRD to shapefile
shape@data$sr <- sr_df$species_richness
shape@data$mrd <- mrd[,1]
shape@data$mrd_z_score <- z
shape@data$mrd_z_score_binary <- shape@data$mrd_z_score
shape@data$mrd_z_score_binary[which(shape@data$mrd_z_score< -1.97)] <- -1
shape@data$mrd_z_score_binary[which(shape@data$mrd_z_score>1.97)] <- 1
shape@data$mrd_z_score_binary[which(shape@data$mrd_z_score>-1.97 & shape@data$mrd_z_score<1.97)] <- 0

cor.test(shape@data$sr, shape@data$mrd, method="s")
cor.test(shape@data$sr, shape@data$mrd_z_score, method="s")


# Transform to data frame for ggplot
library(sf)
shp <- st_as_sf(shape)


# SR MAP ##########################################################################
ggplot(shp) + 
  geom_sf(aes(fill = sr), lwd=0.1) + 
  scale_fill_viridis_c(option = "plasma", trans = "sqrt")+ #, 
#  scale_fill_gradient("Species richness", guide="colorbar", low="white", high="darkred",
#                      limits=c(min(sr_df$species_richness),max(sr_df$species_richness)))+
  guides(fill = guide_colourbar(barwidth = 20, direction="horizontal"))+ # stretch that colorbar
  theme_void()+
  theme(legend.position = "bottom")
ggsave(filename=paste0("results/sr_map", gsub(" |:", "_", date()), ".png"), dpi=600, width=10, height=7)



# MRD MAP ##########################################################################
ggplot(shp) + 
  geom_sf(aes(fill = mrd), lwd=0.1) + 
  scale_fill_viridis_c(option = "plasma")+ #, 
  guides(fill = guide_colourbar(barwidth = 20, direction="horizontal"))+ # stretch that colorbar
  theme_void()+
  theme(legend.position = "bottom")
ggsave(filename=paste0("results/mrd_map", gsub(" |:", "_", date()), ".png"), dpi=600, width=10, height=7)


ggplot(shp) + 
  geom_sf(aes(fill = mrd_z_score), lwd=0.1) + 
  scale_fill_viridis_c(option = "plasma")+ #, 
  guides(fill = guide_colourbar(barwidth = 20, direction="horizontal"))+ # stretch that colorbar
  theme_void()+
  theme(legend.position = "bottom")
ggsave(filename=paste0("results/mrd_z_score_map", gsub(" |:", "_", date()), ".png"), dpi=600, width=10, height=7)



# check how many obs are inside 2xSD = < or > 1.96 z-score
shp <- na.omit(shp)
ggplot(shp) + 
  geom_sf(aes(fill = factor(mrd_z_score_binary)), lwd=0.1, na.rm = TRUE) + 
#  scale_fill_viridis_d(option = "plasma")+ 
  scale_fill_manual(values = c("lightblue", "grey", "darkorange")) + 
  theme_void()+
  theme(legend.position = "bottom")
ggsave(filename=paste0("results/mrd_z_score_binary_map", gsub(" |:", "_", date()), ".png"), dpi=600, width=10, height=7)

table(shp$mrd_z_score_binary)/nrow(shp) # 30% lower, 40% higher, 30% null


cols <- c("blue", "grey", "red")
continent <- unique(shp$CONTINENT)
ggplot(shp, aes(sr, mrd_z_score, col=factor(mrd_z_score_binary)))+
  geom_point()+
  scale_x_log10()+
  geom_smooth(method="lm")

ggplot(shp, aes(sr, mrd_z_score, group=factor(CONTINENT), label=LEVEL_3_CO))+
#  geom_point()+
  geom_text(aes(label=LEVEL_3_CO, col=factor(mrd_z_score_binary)),hjust=0, vjust=0)+
  scale_x_log10()+
  geom_smooth(method="lm")+
  facet_wrap(~CONTINENT, scales = "free")


#Get percentage higher lower for each continent
perc_tab <- function(x){table(factor(x))/length(x)}
tapply(shp$mrd_z_score_binary, shp$CONTINENT, perc_tab)



# TRASH

# workspace size etc
# system('grep MemTotal /proc/meminfo') # 16GB
# size = 0
# for (x in ls() ){
#   thisSize = object.size(get(x))
#   size = size + thisSize
#   message(x, " = ", appendLF = F); print(thisSize, units='auto')
# }
# message("total workspace is ",appendLF = F); print(size, units='auto')
