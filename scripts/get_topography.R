# topography

# antonelli 2018 uses this:
# Topographic relief measured as mean of elevation range values (max-min) within 2.5 km radii for each point

# topographic heterogeneity
# get mean elevation, elevation range, ruggedness a la Riley et al. 1999?

# world clim has 30 seconds resolution


rm(list=ls())

library(raster)
library(vegan)
library(rgdal)
library(spatialEco)

shape <- readOGR("shapefile/level3_mod.shp")
elev <- raster("data/worldclim/wc2.1_30s_elev.tif")

#library(rayshader)
#elmat = raster_to_matrix(elev)

#We use another one of rayshader's built-in textures:
# elmat %>%
#   sphere_shade(texture = "desert") %>%
#   plot_map()


Sys.time()


elev_mean <- c()
elev_sd <- c()
elev_n <- c()
for(i in 1:nrow(shape@data)){ #nrow(shape@data)
  shape_sub <- subset(shape, shape$LEVEL_3_CO==shape$LEVEL_3_CO[[i]]) # 
  rest <- raster::extract(elev, shape_sub)
  rest <- na.omit(rest[[1]])
  elev_mean <- c(elev_mean, mean(rest))
  elev_sd <- c(elev_sd, sd(rest))
  elev_n <- c(elev_n, length(rest))
  if(!i%%1)cat(i,"\r")
}
Sys.time()


res <- data.frame(country=shape$LEVEL_3_CO, elev_mean, elev_sd, elev_n)


# next part takes a while
Sys.time()
ruggedness_mean <- c()
for(i in 1:nrow(shape@data)){
  shape_sub <- subset(shape, shape$LEVEL_3_CO=="JNF") #shape$LEVEL_3_CO[[i]]
  elev_sub <- crop(elev, extent(shape_sub))
  elev_sub2 <- mask(elev_sub, shape_sub)
  ruggedness <- tri(elev_sub2, exact=FALSE)
  rug <- getValues(ruggedness)
  rest <- na.omit(rug)
  ruggedness_mean <- c(ruggedness_mean, mean(rest))
  if(!i%%1)cat(i,"\r")
}
Sys.time()

res$ruggedness_mean <- ruggedness_mean

saveRDS(res, file="topography.rds")



##### RESULTS ############################################################

# run on server 5h
res <- readRDS("data/topography.rds")


# check
# Transform to data frame for ggplot
library(sf)
library(ggplot2)
shape$elev_mean <- res$elev_mean
shape$elev_sd <- res$elev_sd
shape$elev_n <- res$elev_n
shape$tri <- res$ruggedness_mean
shp <- st_as_sf(shape)


# TOPO MAPS ##########################################################################
library(ggplotify)
gg.elev <- as.ggplot(~plot(elev), scale = 1, hjust = 0, vjust = 0)


plot.mean <- ggplot(shp) +
  geom_sf(aes(fill = elev_mean), lwd=0.1)+
  theme_void()+
  scale_fill_continuous(trans="sqrt")

plot.sd <- ggplot(shp) +
  geom_sf(aes(fill = elev_sd), lwd=0.1)+
  theme_void()+
  scale_fill_continuous(trans="sqrt")

plot.tri <- ggplot(shp) +
  geom_sf(aes(fill = tri), lwd=0.1)+
  theme_void()+
  scale_fill_continuous(trans="sqrt")

grid.arrange(gg.elev, plot.mean, plot.sd, plot.tri)



