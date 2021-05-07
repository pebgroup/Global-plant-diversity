# Script to get biome types for botanical countries. 
# Getting percentage of country classified as tropical rainforest / tropics
# Biome maps: tropical rainforest following Corlett & Primack (2011), Olsen 2001


# check projections
proj4string(olson) == 
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

