## MODIFY the shapefile
## the goal is to make the tiny islands slightly bigger to compensate for raster + shapefile being slightly off

library(geosphere)
library(sf)


shape <- readOGR("shapefile/level3.shp")
shape$area <- areaPolygon(shape)
table(shape$area>1e+9)
shape$size <- as.numeric(shape$area>1e+9)+1
plot(shape, border=c("red","blue")[shape$size])
# 1e+9 is chosen based on subjective observation: covers exclusively small islands (<1000km^2)

shape2 <- st_as_sf(shape)
buff <- st_buffer(shape2[shape2$size ==1,], dist = 0.2)
plot(shape[shape$size==1,])
plot(st_geometry(buff), add=TRUE, border="pink")

# rewrite shapefile
shape2[shape2$size ==1,] <- buff
shape3 <- as_Spatial(shape2)
shape3$area <- areaPolygon(shape3)
shape$area==shape3$area
plot(log(shape$area), log(shape3$area))

shape3 <- shape3[,-(6:7)]
writeOGR(shape3, "shapefile/level3_mod.shp", driver="ESRI Shapefile", layer = 'level3_buffered')

test <- readOGR("shapefile/level3_mod.shp")
