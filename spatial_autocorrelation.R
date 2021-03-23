#load(url("http://github.com/mgimond/Spatial/raw/master/Data/moransI.RData"))
class(s1)
library(tmap)
tm_shape(s1) + tm_polygons(style="quantile", col = "Income") +
  tm_legend(outside = TRUE, text.size = .8) 

### NEIGHBOURS BASED ##############################################

# define neighbouring polygons
library(spdep)
nb <- poly2nb(s1, queen=TRUE)
nb[[1]]
s1$NAME[c(2,3,4,5)] # nighbouring counties

# you can assign weights to each neighbor. use here style="w" for equal weights. this can cause spatial edge effects, consider using "B" instead

lw <- nb2listw(nb, style="W", zero.policy=TRUE)
lw$weights[1]

# get spatially lagged values: average neighbour values of interest:
Inc.lag <- lag.listw(lw, s1$Income)
Inc.lag

# Create a regression model
M <- lm(Inc.lag ~ s1$Income)

# Plot the data
plot( Inc.lag ~ s1$Income, pch=20, asp=1, las=1)
abline(lm(Inc.lag ~ s1$Income), col="red")

# check if the slope is sign different from zero
# ...

# Morans I function with analytical p-value
moran.test(s1$Income,lw)

# MOrans I with MC simulation method:
MC<- moran.mc(s1$Income, lw, nsim=599)
MC
plot(MC, main="", las=1)


### DISTANCE BASED ##############################################


coo <- coordinates(s1) # get center of each polygon
S.dist  <-  dnearneigh(coo, 0, 50000) # define search radius as 50000 meters

lw <- nb2listw(S.dist, style="W",zero.policy=T) 
MI  <-  moran.mc(s1$Income, lw, nsim=599,zero.policy=T) 
plot(MI, las=1)
MI




### SPECIES RICHNESS #################################################

library(rgdal)
s1 <- readOGR("shapefile/level3.shp")
shp <- readRDS("data/shp_object_fin_analysis.RDS")
#library(sf)
#s1 <- as_Spatial(shp)
s1$SR <- shp$sr

coo <- coordinates(s1) # get center of each polygon

S.dist  <-  dnearneigh(coo, 0, 1000, longlat = TRUE) # define search radius. we got decimal degrees, so the number refers to kilometers

#nb <- poly2nb(s1, queen=TRUE)
s1$LEVEL_3_CO[S.dist[[1]]] # neighboring counties


lw <- nb2listw(S.dist, style="W",zero.policy=T) 
MI  <-  moran.mc(s1$SR, lw, nsim=599,zero.policy=T) 
plot(MI, las=1)
MI

# run for different distance classes
# Kreft & Jetz use for their global analysis distance classes between 0 and 5000km, in ca 100km steps

moran <- data.frame(dist.class = seq(100, 10000, 100),
                    moransI = NA,
                    moransp = NA)
for(i in 1:length(moran$dist.class)){
  S.dist  <-  dnearneigh(coo, 0, moran$dist.class[i], longlat = TRUE)
  lw <- nb2listw(S.dist, style="W",zero.policy=T) 
  MI <- moran.mc(s1$SR, lw, nsim=599,zero.policy=T) 
  moran$moransI[i] <- MI$statistic
  moran$moransp[i] <- MI$p.value
  if(!i%%1)cat(i,"\r")
}

plot(moran$moransI~moran$dist.class, 
     type="p", pch=20,
     col=c("red", "black")[as.numeric(moran$moransp<0.05)+1],
     xlab="Distance class (km)",
     ylab="Moran's I")





