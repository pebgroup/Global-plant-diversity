# Analysis script for WCSP paper
# loaded data:
# comm: community matrix. rows represents regions, columns WCSP species IDs
# phytoregions shapefile
# mean root distance


wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str()))  
library(rgdal)
library(rgeos)
library(geosphere)
library(sf)

sr <- readRDS("../processed_data/comm_September2021.rds")
mrd <- readRDS("../processed_data/mrd_Sep2021.rds")
divr <- readRDS("../processed_data/DR_rates_Feb2022.rds")
shape <- readOGR("../data/shapefile_bot_countries/level3.shp")
soil <- readRDS("../processed_data/soil.rds")
climate <- readRDS("../processed_data/climate.rds")
pclimate <- readRDS("../processed_data/paleo_climate.rds")
topography <- readRDS("../processed_data/topography.rds")
biomes <- readRDS("../processed_data/biomes_olson.rds")

# get country centroids
trueCentroids = as.data.frame(gCentroid(shape,byid=TRUE))

# get species richness
sr_df <- data.frame(species_richness = rowSums(sr), region =row.names(sr))


# Add SR, MRD and centroids to shapefile ----------------------------------
shape@data$sr <- sr_df$species_richness
shape@data$mrd <- mrd$mrd
shape@data$mdr <- divr$MDR
#shape@data$lng <- trueCentroids[,1]
#shape@data$lat <- trueCentroids[,2]
rm(sr, mrd, divr)


# Add environmental variables to shape file -------------------------------
shape@data$soil <- soil$number_soils
shape@data <- cbind(shape@data, climate)
shape@data <- cbind(shape@data, pclimate[,grep("ano_mean", names(pclimate))])
  # try also absolute Miocene anomalies 
  shape@data$mio_mat_ano_abs <- abs(shape@data$mio_mat_ano_mean)
  shape@data$mio_pre_ano_abs <- abs(shape@data$mio_pre_ano_mean)
shape@data$elev_range <- topography$elevation_range
shape@data$tri <- topography$ruggedness_mean
rm(soil, climate, pclimate, topography)

# add area
shape$area <- areaPolygon(shape)


# Add biomes --------------------------------------------------------------
biome_names <- c("sub_trop_mbf", "sub_trop_dbf", "sub_trop_cf", "temp_bmf", "temp_cf",
  "boreal_f_taiga", "sub_trop_gss", "temp_gss", "flooded_gs", "mont_gs", "tundra",
  "medit_fws", "deserts_x_shrub", "mangroves")

names(biomes)[grepl("^X", names(biomes))] <- biome_names
shape@data <- cbind(shape@data, biomes[,-grep("country",names(biomes))])


# Transform to data frame for ggplot
shp <- st_as_sf(shape)
saveRDS(shp, "../processed_data/shp_object_fin_analysis.RDS")


