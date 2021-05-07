# Script to get topography for botanical countries. 
# Using Terrain ruggedness index (TRI) (Riley et al. 1999)

# TRI
Sys.time()
ruggedness_mean <- c()
for(i in 1:nrow(shape@data)){
  shape_sub <- subset(shape, shape$LEVEL_3_CO==shape$LEVEL_3_CO[[i]]) 
  elev_sub <- crop(elev, extent(shape_sub)) # crops to the extent
  elev_sub2 <- mask(elev_sub, shape_sub) # masks non rectangular patterns
  ruggedness <- tri(elev_sub2, exact=FALSE)
  rug <- getValues(ruggedness)
  rest <- na.omit(rug)
  ruggedness_mean <- c(ruggedness_mean, mean(rest))
  if(!i%%1)cat(i,"\r")
}
Sys.time()

res <- data.frame(country=shape$LEVEL_3_CO, ruggedness_mean)



