# topography

# antonelli 2018 uses this:
# Topographic relief measured as mean of elevation range values (max-min) within 2.5 km radii for each point

# topographic heterogeneity
# get mean elevation, elevation range, range evenness?

# world clim has 30 seconds resolution


library(raster)
library(vegan)
library(rgdal)

shape <- readOGR("shapefile/level3.shp")
elev <- raster("data/wc2.1_30s_elev.tif")


Sys.time()
elev_mean <- c()
elev_sd <- c()
elev_n <- c()
for(i in 1:nrow(shape@data)){ #nrow(shape@data)
  shape_sub <- subset(shape, shape$LEVEL_3_CO==shape$LEVEL_3_CO[[i]])
  rest <- extract(elev, shape_sub)
  rest <- na.omit(rest[[1]])
  elev_mean <- c(elev_mean, mean(rest))
  elev_sd <- c(elev_sd, sd(rest))
  elev_n <- c(elev_n, length(rest))
  if(!i%%1)cat(i,"\r")
}
Sys.time()

res <- data.frame(country=shape$LEVEL_3_CO, elev_mean, elev_sd, elev_n)

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
shp <- st_as_sf(shape)


# SOIL MAPS ##########################################################################
ggplot(shp) +
  geom_sf(aes(fill = elev_n), lwd=0.1)+
  theme_void()+
  scale_fill_continuous(trans="sqrt")


