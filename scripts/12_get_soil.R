# Script to get number of soil types for botanical countries. 
# Saves a "soil.rds" file
# Input raster layer has been reduced in resolution using QGIS (0.00832 = approx 1km) for better computation time


## RUNTIME: approx 5h

#### Data sources #######################################################
# # https://www.isric.org/web-coverage-services-wcs
## https://files.isric.org/soilgrids/latest/data/wrb/MostProbable.vrt
## vrt for most probable downloaded and translated to geoTiff in QGIS with a lower resolution.

wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str())) 
library(raster)
library(rgdal)

shape <- readOGR("../data/shapefile_bot_countries/level3_mod.shp")
soil <- raster("../data/soil_raster_layer_000832.tif")

Sys.time()
number_soils <- c()
for(i in 1:nrow(shape@data)){ 
  shape_sub <- subset(shape, shape$LEVEL_3_CO==shape$LEVEL_3_CO[[i]])
  rest <- extract(soil, shape_sub)
  rest <- na.omit(rest[[1]])
  number_soils <- c(number_soils, length(unique(rest)))
  if(!i%%1)cat(i,"\r")
}
Sys.time()

res <- data.frame(country=shape$LEVEL_3_CO, number_soils)

saveRDS(res, file="../processed_data/soil.rds")
