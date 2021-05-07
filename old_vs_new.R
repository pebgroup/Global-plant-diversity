# compare old and new environmental variables

 clim <- readRDS(file="processed_data/climate.rds")
 clim_o <- readRDS(file="data/climate.rds")
 
 # soil <- readRDS(file="processed_data/soil.rds")
 # soil_o <- readRDS(file="data/soil.rds")
 
 topo <- readRDS(file="processed_data/topography.rds")
 topo_o <- readRDS(file="data/topography.rds")
 
 biom <- readRDS(file="processed_data/biomes_olson.rds")
 biom_o <- readRDS(file="data/biomes_olson.rds")

 # climate
par(mfrow=c(4,4))
for(i in 1:ncol(clim)){
  plot(clim[,i], clim_o[,i])
  print(cor.test(clim[,i], clim_o[,i]))
} 
 
 # topo
dev.off()
plot(topo$ruggedness_mean, topo_o$ruggedness_mean)
cor.test(topo$ruggedness_mean, topo_o$ruggedness_mean)

# biomes
par(mfrow=c(4,4))
for(i in 2:ncol(biom)){
  plot(biom[,i], biom_o[,i])
  print(cor.test(biom[,i], biom_o[,i]))
} 

