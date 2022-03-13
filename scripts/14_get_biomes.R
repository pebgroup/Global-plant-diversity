# Script to get biome types for botanical countries. 
# Getting percentage of country classified as tropical rainforest / tropics
# Biome maps: tropical rainforest following Corlett & Primack (2011), Olsen 2001


wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str())) 
library(raster)
library(rgdal)
library(rgeos)

shape <- readOGR("../data/shapefile_bot_countries/level3.shp")
olson <- readOGR("../data/shapefile_biomes/wwf_terr_ecos.shp")


# check projections
proj4string(olson)
proj4string(shape)


# remove rock and ice biomes
olson <- olson[!olson$BIOME %in% c(98,99),]


Sys.time()
res <- matrix(NA, nrow=nrow(shape), ncol=length(unique(olson@data$BIOME)))
for(i in 1:nrow(shape@data)){
  shape_sub <- subset(shape, shape$LEVEL_3_CO==shape$LEVEL_3_CO[[i]]) # try for CLM shape$LEVEL_3_CO[[i]]
  shape_sub <- gBuffer(shape_sub, byid=TRUE, width=0)
  # intersect from raster package
  for(j in 1:ncol(res)){
    olson_sub <- gBuffer(olson[olson$BIOME==j,], byid=TRUE, width=0) # throwing on a zero buffer to avoid ring self-intersection issues
    pi <- intersect(shape_sub, olson_sub)
    # Extract areas from polygon objects then attach as attribute
    if(class(pi)!="NULL"){
      trop_area <- sum(area(pi))
      shape_sub$area <- area(shape_sub)
      res[i,j] <- trop_area / shape_sub$area
    }else{res[i,j] <- 0}
  }
  if(!i%%1)cat(i,"\r")
}
Sys.time()

res.df <- data.frame(country=shape$LEVEL_3_CO, res)

saveRDS(res.df, file="../processed_data/biomes_olson.rds")



# Canopy ------------------------------------------------------------------

#dATA: P. Potapov, X. Li, A. Hernandez-Serna, A. Tyukavina, M.C. Hansen, A.
#Kommareddy, A. Pickens, S. Turubanova, H. Tang, C. E. Silva, J. Armston, R.
#Dubayah, J. B. Blair, M. Hofton (2020).
#https://doi.org/10.1016/j.rse.2020.112165
# 0-60 Forest canopy height, meters
# 101 Water
# 102 Snow/ice
# 103 No data

can <- raster("../data/canopy_height/combined_tif.tif")
range(getValues(can))
val <- 103/255*(getValues(can))
can <- setValues(can, val)
val[val>60] <- NA # remove water+ice
can <- setValues(can, val)

can_height <- c()
can_range <- c()
for(i in 1:nrow(shape@data)){
  shape_sub <- subset(shape, shape$LEVEL_3_CO==shape$LEVEL_3_CO[[i]])
  can_sub <- crop(can, extent(shape_sub)) # crops to the extent
  can_sub2 <- mask(can_sub, shape_sub) # masks non rectangular patterns
  can_height_tmp <- cellStats(can_sub2, "mean", na.rm=TRUE)
  can_height <- c(can_height, can_height_tmp)
  can_range_tmp <- cellStats(can_sub2, "range")
  can_range <- c(can_range, diff(can_range_tmp))
  if(!i%%1)cat(i,"\r")
}

hist(can_height)
hist(can_range)

res <- data.frame(level3=shape$LEVEL_3_CO, can_height, can_range)
saveRDS(res, file="../processed_data/canopy.rds")

