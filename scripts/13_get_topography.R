# Script to get topography for botanical countries. 
# Using Terrain ruggedness index (TRI) (Riley et al. 1999)

# TRI
Sys.time()
ruggedness_mean <- c()
for(i in 1:nrow(shape@data)){
  shape_sub <- subset(shape, shape$LEVEL_3_CO==shape$LEVEL_3_CO[[i]]) 
  elev_sub <- crop(elev, extent(shape_sub)) # crops to the extent
  elev_sub2 <- mask(elev_sub, shape_sub) # masks non rectangular patterns
  ruggedness <- tri(elev_sub2, exact=FALSE, s=5)
  rug <- getValues(ruggedness)
  rest <- na.omit(rug)
  ruggedness_mean <- c(ruggedness_mean, mean(rest))
  if(!i%%1)cat(i,"\r")
}
Sys.time()

# Elevational range
# extract_metric_per_country <- function(shape, env_data, metric){
#   # shape = a polygon from a shapefile
#   # env_data = world raster layer
#   # metric = character value of metric you need to extract (e.g. mean, sd, range)
#   # empty vector to take in results c()
#   elev_sub <- crop(env_data, extent(shape)) # crops to the extent
#   elev_sub2 <- mask(elev_sub, shape) # masks non rectangular patterns
#   elev_range <- cellStats(elev_sub2, metric)
#   if(metric=="range"){
#   elev_range <- diff(elev_range)
#   }
#   elevation_range <- c(elevation_range, elev_range)
#   return(elevation_range)
# }
# sapply(shape, extract_metric_per_country, shape=shape, env_data=elev, metric="range")
# apply(lapply(list.files(pattern="*.shp"), readShapePoly), gArea)


elevation_range <- c()
for(i in 1:nrow(shape@data)){
  shape_sub <- subset(shape, shape$LEVEL_3_CO==shape$LEVEL_3_CO[[i]])
  elev_sub <- crop(elev, extent(shape_sub)) # crops to the extent
  elev_sub2 <- mask(elev_sub, shape_sub) # masks non rectangular patterns
  elev_range <- cellStats(elev_sub2, "range")
  #rug <- getValues(ruggedness)
  #rest <- na.omit(rug)
  elevation_range <- c(elevation_range, diff(elev_range))
  if(!i%%1)cat(i,"\r")
}


res <- data.frame(country=shape$LEVEL_3_CO, ruggedness_mean, elevation_range)
saveRDS(res, file="processed_data/topography.rds")

# 
# # check some results
# plot(res$ruggedness_mean, res$elevation_range)
# cor.test(res$ruggedness_mean, res$elevation_range)
# 
# par(mfrow=c(2,2))
# ## countries with low elevation range but high ruggedness
# plot(res$ruggedness_mean, res$elevation_range, pch=NA)
# text(res$ruggedness_mean, res$elevation_range, labels=res$country, 
#       col=c("black", "red")[as.numeric(res$country=="CLM")+1])
# 
# ## JNF
# shape[shape$LEVEL_3_CO=="JNF",]
# shape_sub <- subset(shape, shape$LEVEL_3_CO=="JNF")
# elev_sub <- crop(elev, extent(shape_sub)) # crops to the extent
# elev_sub2 <- mask(elev_sub, shape_sub) # masks non rectangular patterns
# plot(elev_sub2, main="JNF, high ruggedness, low elevational range")
# 
# ## KAZ
# shape_sub <- subset(shape, shape$LEVEL_3_CO=="KAZ")
# elev_sub <- crop(elev, extent(shape_sub)) # crops to the extent
# elev_sub2 <- mask(elev_sub, shape_sub) # masks non rectangular patterns
# plot(elev_sub2, main="KAZ, low ruggedness, high elevational range")
# 
# ## CLM
# shape_sub <- subset(shape, shape$LEVEL_3_CO=="CLM")
# elev_sub <- crop(elev, extent(shape_sub)) # crops to the extent
# elev_sub2 <- mask(elev_sub, shape_sub) # masks non rectangular patterns
# plot(elev_sub2, main="CLM")

# 
# library(rayshader)
# #And convert it to a matrix:
# elmat = raster_to_matrix(elev_sub2)
# #We use another one of rayshader's built-in textures:
# elmat %>%
#   sphere_shade(texture = "desert") %>%
#   plot_3d()
# 
# elmat %>%
#   sphere_shade(texture = "desert") %>%
#   add_water(detect_water(elmat), color = "desert") %>%
#   add_shadow(ray_shade(elmat, zscale = 3), 0.5) %>%
#   add_shadow(ambient_shade(elmat), 0) %>%
#   plot_3d(elmat, zscale = 10, fov = 0, theta = 135, zoom = 0.75, phi = 45, windowsize = c(1000, 800))
# Sys.sleep(0.2)
# render_snapshot()
