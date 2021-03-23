# Script to get human footprint index by Venter et al 2016


## RUNTIME: 3 hours on server

#### Data sources #######################################################
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.052q5


library(raster)
library(vegan)
library(rgdal)

shape <- readOGR("shapefile/level3.shp")
hfp <- raster("data/HFP2009.tif")
proj4string(hfp)
proj4string(shape)

shape <- spTransform(shape, CRS(proj4string(hfp)))

#plot(shape)
#plot(hfp, add=TRUE)



# subset example
# shape_sub <- subset(shape, shape$LEVEL_3_CO=="RUW")
# plot(shape_sub)
# plot(soil, add=TRUE) # this does not plot all pixels, memory issues
# rest <- extract(soil, shape_sub)
# length(unique(rest[[1]]))
# # checks out
# table(rest[[1]])
# diversity(table(rest[[1]]), index="simp")
# H <- diversity(table(rest[[1]])) # shannon diversity
# # evenness: shannon/shannonmax==number of species
# H/log(length(unique(rest_clean)))
# 



Sys.time()
hfp90 <- c()
hfp_mean <- c()
for(i in 1:nrow(shape@data)){ #nrow(shape@data)
  shape_sub <- subset(shape, shape$LEVEL_3_CO==shape$LEVEL_3_CO[[i]])
  rest <- extract(hfp, shape_sub)
  rest <- na.omit(rest[[1]])
  hfp90 <- c(hfp90, quantile(rest, probs=c(0.95))[[1]])
  hfp_mean <- c(hfp_mean, mean(rest))
  if(!i%%1)cat(i,"\r")
}
Sys.time()

res <- data.frame(country=shape$LEVEL_3_CO, hfp90, hfp_mean)

saveRDS(res, file="hfp.rds")

# run on server 5h
res <- readRDS("data/hfp.rds")


# check
# Transform to data frame for ggplot
library(sf)
library(ggplot2)
shape$hfp90 <- res$hfp90
shape$hfp_mean <- res$hfp_mean
shp <- st_as_sf(shape)


# HFP MAPS ##########################################################################
plot(shp$hfp90, shp$hfp_mean)
cor(na.omit(shp$hfp90), na.omit(shp$hfp_mean), method="s")

ggplot(shp) +
  geom_sf(aes(fill = hfp_mean), lwd=0.1)+
  theme_void()

ggplot(shp) +
  geom_sf(aes(fill = hfp90), lwd=0.1)+
  theme_void()





