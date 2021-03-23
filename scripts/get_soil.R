# Script to get soil types for botanical countries. 
# Getting number of soil types, simpsons diversity and evenness.
# Saves a "data/soil.rds" file: dataframe containing number of soil types, simpsons diversity and evenness.
# Raster layer has been reduced in resolution using QGIS (0.00832 = approx 1km)


## RUNTIME: approx 5h on 128 gb RAM

#### Data sources #######################################################
# # https://www.isric.org/web-coverage-services-wcs
## https://files.isric.org/soilgrids/latest/data/wrb/MostProbable.vrt
## vrt for most probable downloaded and translated to geoTiff in QGIS with a lower resolution.


library(raster)
library(vegan)
library(rgdal)

shape <- readOGR("shapefile/level3_mod.shp")
soil <- raster("soil_raster_layer_000832.tif")
#soil <- raster("data/soilgrids/soil_raster_layer_000832.tif")

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
number_soils <- c()
simp_soil <- c()
even_soil <- c()
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

# run on server 5h
res <- readRDS("data/soil.rds")


# check
# Transform to data frame for ggplot
library(sf)
library(ggplot2)
shape$soil <- res$number_soils
shape$soil_simp <- res$simp_soil
shape$soil_even <- res$even_soil
shp <- st_as_sf(shape)


# SOIL MAPS ##########################################################################
ggplot(shp) +
  geom_sf(aes(fill = soil), lwd=0.1)+
  theme_void()

ggplot(shp) +
  geom_sf(aes(fill = soil_simp), lwd=0.1)+
  theme_void()

ggplot(shp) +
  geom_sf(aes(fill = soil_even), lwd=0.1)+
  theme_void()


  
  
  
  
  
  
  
  
  
  


# WASTE BIN ##################################################################
# 
# library(XML)
# library(rgdal)
# library(gdalUtils)
# 
# voi = "nitrogen" # variable of interest
# depth = "5-15cm"
# quantile = "Q0.5"
# 
# voi_layer = paste(voi,depth,quantile, sep="_") # layer of interest 
# wcs_path = paste0("https://maps.isric.org/mapserv?map=/map/",voi,".map") # Path to the WCS. See maps.isric.org
# wcs_service = "SERVICE=WCS"
# wcs_version = "VERSION=2.0.1"
# 
# 
# 
# 
# wcs_request = "DescribeCoverage" 
# 
# wcs = paste(wcs_path, wcs_service, wcs_version, wcs_request, sep="&")
# 
# l1 <- newXMLNode("WCS_GDAL")
# l1.s <- newXMLNode("ServiceURL", wcs, parent=l1)
# l1.l <- newXMLNode("CoverageName", voi_layer, parent=l1)
# 
# # Save to local disk
# xml.out = "./sg.xml"
# saveXML(l1, file = xml.out)
# gdalinfo("./sg.xml")
# 
# 
# 
# 
# #Manual download from https://files.isric.org/soilgrids/latest/data/wrb/
# library(raster)
# test <- raster("data/soilgrids/MostProbable.vrt")
# test
# 
# test2 <- readGDAL("data/soilgrids/MostProbable.vrt")
# # this file is way too big
# 
# 
# 
# 
# library(gdalUtils)
# gdalwarp(t_srs="EPSG:4326", multi=TRUE, wm=200, 
#          co=c("BIGTIFF=YES", "COMPRESS=DEFLATE", "TILED=TRUE"),
#          tr=c(0.25,0.25), # Desired output resolution
#          verbose=T,
#          "/vsicurl?max_retry=3&retry_delay=1&list_dir=no&url=https://files.isric.org/soilgrids/latest/data/wrb/Acrisols.vrt", # Input VRT
#          "Acrisols.tif") # Output file



############

# Digital Soil Map of the World. Based on Harmonized Soil Map, using the FAO/UNESCO soil types
# Map download: http://www.fao.org/geonetwork/srv/en/main.home?uuid=446ed430-8383-11db-b9b2-000d939bc5d8
# Legend for soiltypes: http://www.fao.org/fileadmin/user_upload/soils/docs/Soil_map_FAOUNESCO/images/Legend_I.jpg
# library(rgdal)
# library(sf)
# soil <- readOGR("data/DSMW/DSMW.shp")
# names(soil@data)
# plot(soil, col=soil@data$DOMSOI)
# 
# shp <- st_as_sf(soil)
# library(ggplot2)
# ggplot(shp) + 
#   geom_sf(aes(fill = DOMSOI), lwd=0.1)+
#   theme_void()
# 
# table(shp$DOMSOI)
# 
# # Get the broader soil types (capital letters)
# split()
# soil.list <- strsplit(as.character(shp$DOMSOI), "")
# broad_soils <-c()
# for(i in 1:length(soil.list)){
#   broad_soils <- c(broad_soils, toString(soil.list[[i]][grep("[A-Z]", soil.list[[i]])]))
# }
# 
# shp$broad_soil <- broad_soils
# ggplot(shp) + 
#   geom_sf(aes(fill = broad_soil), lwd=0.1)+
#   theme_void()
# 
# shape <- readOGR("shapefile/level3.shp")
# proj4string(soil) <- proj4string(shape)
# # the shapefile seems corrupt :]
# 
# library(rgeos)
# 
# for(i in 1:length(nrow(shape@data))){
#   shape.sub <- subset(shape, shape@data$LEVEL_3_CO==shape@data$LEVEL_3_CO[1])
#   gIntersection(shape.sub, soil)
# }
# 
# 
# 
# #library(Hmisc)
# #test <- mdb.get("data/HWSD.mdb")
# 
# library(raster)
# test <- raster("data/HWSD_RASTER/hwsd.bil")
# 
# 
# 
# 
# 
# 
# # SOILGRIDS - does not work...
# library(rgdal)
# library(gdalUtils)
# 
# bb=c(-337500.000,1242500.000,152500.000,527500.000) # Example bounding box (homolosine)
# igh='+proj=igh +lat_0=0 +lon_0=0 +datum=WGS84 +units=m +no_defs' # proj string for Homolosine projection
# gdal_translate(
#   '/vsicurl?max_retry=3&retry_delay=1&list_dir=no&url=https://files.isric.org/soilgrids/latest/data/wrb/MostProbable.vrt', 
#   "./crop_roi_igh_r.tif",
#   tr=c(250,250),
#   projwin=bb,
#   projwin_srs =igh,
#   verbose=TRUE
# )
# 
