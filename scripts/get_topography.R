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
elev <- raster("data/worldclim/wc2.1_30s_elev.tif")


Sys.time()
elev_mean <- c()
elev_sd <- c()
even_elev <- c()
for(i in 1:nrow(shape@data)){ #nrow(shape@data)
  shape_sub <- subset(shape, shape$LEVEL_3_CO==shape$LEVEL_3_CO[[i]])
  rest <- extract(soil, shape_sub)
  rest <- na.omit(rest[[1]])
  number_soils <- c(number_soils, length(unique(rest)))
  simp_soil <- c(simp_soil, diversity(table(rest), index="simp"))
  even_soil <- c(even_soil, diversity(table(rest))/log(length(unique(rest))))
  if(!i%%1)cat(i,"\r")
}
Sys.time()

res <- data.frame(country=shape$LEVEL_3_CO, number_soils, simp_soil, even_soil)

saveRDS(res, file="soil.rds")