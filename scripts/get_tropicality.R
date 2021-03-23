# Script to get biome types for botanical countries. 
# Getting percentage of country classified as tropical rainforest / tropics
# Biome maps: tropical rainforest following Corlett & Primack (2011), Olsen 2001


library(raster)
library(vegan)
library(rgdal)


#### Data sources #######################################################


shape <- readOGR("shapefile/level3.shp")
olson <- readOGR("shapefile/olson/wwf_terr_ecos.shp")
corlett <- readOGR("shapefile/corlett/corlett_TRF_ext.shp")

#plot(olson[olson$BIOME==1,], color="black") # biome 1 = tropical and subtropical moist broadleaf forest (TRF def)
#plot(corlett) # corlett = binary, just tropics


# check projections
proj4string(olson)
proj4string(corlett)
proj4string(shape)

# # simplyfiy olson
# olson <- olson[olson$BIOME==1,]
# olson <- aggregate(olson) # yes there are warnings, but looks ok
# proj4string(shape) <- proj4string(olson)
# 
# plot(shape)
# plot(olson, add=TRUE, col="pink", border="pink") # this does not plot all pixels, memory issues
# plot(corlett, add=TRUE, col="yellow", border="yellow")
# 
# 
# 
# Sys.time()
# percent_olson <- c()
# percent_corlett <- c()
# for(i in 1:nrow(shape@data)){
#   shape_sub <- subset(shape, shape$LEVEL_3_CO==shape$LEVEL_3_CO[[i]]) # try for CLM 
#   # intersect from raster package
#   pi <- intersect(shape_sub, olson)
#   #plot(shape_sub, axes=T); plot(olson, add=T, border="green"); plot(pi, add=T, col='red')
#   
#   # Extract areas from polygon objects then attach as attribute
#   if(class(pi)!="NULL"){
#     pi$area <- area(pi)
#     shape_sub$area <- area(shape_sub)
#     percent_olson <- c(percent_olson, pi$area / shape_sub$area)
#   }else{percent_olson <- c(percent_olson, 0)}
#   if(!i%%1)cat(i,"\r")
# }
# Sys.time()
# 
# res <- data.frame(country=shape$LEVEL_3_CO, percent_olson)
# 
# saveRDS(res, file="biomes.rds")








# remove rock and ice
olson <- olson[!olson$BIOME %in% c(98,99),]
#olson <- aggregate(olson) # yes there are warnings, but looks ok
#proj4string(shape) <- proj4string(olson)

plot(shape)
plot(olson, add=TRUE) # this does not plot all pixels, memory issues
#plot(corlett, add=TRUE, col="yellow", border="yellow")



Sys.time()
res <- matrix(NA, nrow=nrow(shape), ncol=length(unique(olson@data$BIOME)))
for(i in 1:nrow(shape@data)){
  shape_sub <- subset(shape, shape$LEVEL_3_CO==shape$LEVEL_3_CO[[i]]) # try for CLM shape$LEVEL_3_CO[[i]]
  shape_sub <- gBuffer(shape_sub, byid=TRUE, width=0)
  # intersect from raster package
  for(j in 1:ncol(res)){
    olson_sub <- gBuffer(olson[olson$BIOME==j,], byid=TRUE, width=0) # throwing on a zero buffer to avoid ring self-intersection issues
    pi <- intersect(shape_sub, olson_sub)
    #plot(shape_sub, axes=T); plot(olson, add=T, border="green"); plot(pi, add=T, col='red')
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

# calc the predominant biome
res.df$predominant_biome <- apply(res, 1, which.max)
res.df$number_biomes <- apply(res, 1, function(x){length(which(x!=0))})


# could also get the "predominant" biome, but are those comparable? could be 90% or 51%

saveRDS(res.df, file="biomes_olson.rds")
