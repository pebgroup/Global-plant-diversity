# rayshader experiments

rm(list=ls())
library(ggplot2)
library(rayshader)
library(raster)
library(rgdal)

shp <- readRDS("data/shp_object_fin_analysis.RDS")

# map <- ggplot(shp) + 
#   geom_sf(aes(fill = sr,), lwd=0.1) + 
#   scale_fill_viridis_c(option = "plasma", trans = "sqrt")+ #, 
#   #  scale_fill_gradient("Species richness", guide="colorbar", low="white", high="darkred",
#   #                      limits=c(min(sr_df$species_richness),max(sr_df$species_richness)))+
#   guides(fill = guide_colourbar(barwidth = 20, direction="horizontal"))+ # stretch that colorbar
#   theme_void()+
#   theme(legend.position = "bottom")

# this takes forever....
#plot_gg(map, width =12, height = 12, zoom = 0.5, multicore = TRUE)
#render_snapshot(clear = TRUE)



# Try with elevation data to show topography index differences
# use BRY + CON

elevation <- raster("data/worldclim/wc2.1_30s_elev.tif")
shape <- readOGR("shapefile/level3.shp")
# crop + mask
bry_ele <- crop(elevation, shape[shape$LEVEL_3_CO=="BRY",])
bry_ele <- mask(bry_ele, shape[shape$LEVEL_3_CO=="BRY",])
#plot(bry_ele)





# convert it to a matrix:
elmat = raster_to_matrix(bry_ele)

elmat %>%
  sphere_shade() %>%
  add_water(detect_water(elmat), color = "desert") %>%
  add_shadow(ray_shade(elmat, zscale = 3), 0.5) %>%
  add_shadow(ambient_shade(elmat), 0) %>%
  plot_3d(elmat, zscale = 30, fov = 0, theta = 0, zoom = 0.6, phi = 45, windowsize = c(1000, 1000))
render_snapshot("bry_elevation.png")
#render_highquality("bry_elevation.png")

# CON

con_ele <- crop(elevation, shape[shape$LEVEL_3_CO=="CON",])
con_ele <- mask(con_ele, shape[shape$LEVEL_3_CO=="CON",])
plot(con_ele)

# convert it to a matrix:
con = raster_to_matrix(con_ele)

con %>%
  sphere_shade() %>% #texture = "desert"
  add_water(detect_water(con), color = "desert") %>%
  add_shadow(ray_shade(con, zscale = 3), 0.5) %>%
  add_shadow(ambient_shade(con), 0) %>%
  plot_3d(con, zscale = 30, fov = 0, theta = 0, zoom = 0.6, phi = 45, windowsize = c(1000, 800))
render_snapshot("con_elevation.png")




## GGPLOT #################################################################

#convert the raster to points for plotting
map.p <- rasterToPoints(bry_ele)
df <- data.frame(map.p)
colnames(df) <- c("Longitude", "Latitude", "MAP")
bry_gg <- ggplot(data=df, aes(y=Latitude, x=Longitude)) +
  geom_raster(aes(fill=MAP)) +
  scale_fill_viridis_c("altitude", limits=c(0,3000))+
  theme_void() +
  coord_equal() +
  theme(legend.position="none")

map.con <- rasterToPoints(con_ele)
df.con <- data.frame(map.con)
colnames(df.con) <- c("Longitude", "Latitude", "MAP")
con_gg <- ggplot(data=df.con, aes(y=Latitude, x=Longitude)) +
  geom_raster(aes(fill=MAP)) +
  scale_fill_viridis_c("altitude", limits=c(0,3000))+
  theme_void() +
  coord_equal() 

library(gridExtra)
png("ruggedness_extremes.png", width=20, height=12, units = "cm", res=300)
grid.arrange(bry_gg, con_gg, ncol=2)
dev.off()

df$country <- "BRY"
df.con$country <- "CON"
comb <- rbind(df, df.con)

ggplot(data=comb, aes(y=Latitude, x=Longitude)) +
  geom_raster(aes(fill=MAP)) +
  scale_fill_viridis_c("altitude")+
  theme_void() +
  coord_equal() #+ 
#  facet_wrap(~country, scales = "free")


plot_gg(con_gg, width = 3.5, multicore = TRUE, windowsize = c(800, 800), 
        zoom = 0.6, phi = 35, theta = 0, sunangle = 225, soliddepth = -100)

