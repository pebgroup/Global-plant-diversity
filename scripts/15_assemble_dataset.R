# Analysis script for WCSP paper
# loaded data:
# comm: community matrix. rows represents regions, columns WCSP species IDs
# phytoregions shapefile
# mean root distance

# get country centroids
trueCentroids = as.data.frame(gCentroid(shape,byid=TRUE))

# get species richness
sr_df <- data.frame(species_richness = rowSums(sr), region =row.names(sr))
#sr_2020_df <- data.frame(species_richness = rowSums(sr_2020), region =row.names(sr_2020))


# Add SR, MRD and centroids to shapefile #########################################
shape@data$sr <- sr_df$species_richness
#shape@data$sr_2020 <- sr_2020_df$species_richness
shape@data$mrd <- mrd$mrd
#shape@data$mrd_2020 <- mrd_2020$mrd
shape@data$MRD.sd <- mrd$mrd_sd
shape@data$lng <- trueCentroids[,1]
shape@data$lat <- trueCentroids[,2]

rm(sr)


#### Add environmental variables to shape file #####################################

shape@data$soil <- soil$number_soils
shape@data <- cbind(shape@data, climate)
# shape@data$elev_mean <- topography$elev_mean
# shape@data$elev_sd <- topography$elev_sd
# shape@data$elev_n <- topography$elev_n
shape@data$tri <- topography$ruggedness_mean

# add area
shape$area <- areaPolygon(shape)


# add biomes #######################################################################3#
biome_names <- c("sub_trop_mbf", "sub_trop_dbf", "sub_trop_cf", "temp_bmf", "temp_cf",
  "boreal_f_taiga", "sub_trop_gss", "temp_gss", "flooded_gs", "mont_gs", "tundra",
  "medit_fws", "deserts_x_shrub", "mangroves")

names(biomes)[grepl("^X", names(biomes))] <- biome_names
shape@data <- cbind(shape@data, biomes[,-1])


# Transform to data frame for ggplot
shp <- st_as_sf(shape)



