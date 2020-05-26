# This script collects data for analysis, and brings them into a useful standard form
# To be sourced from analysis scripts; assumes correct working directory is set. 

library(rgdal)
library(rgeos)
library(raster) 

# data location
data_location <- "plant_sr_data/"
# shapefile location
shapefile_location <- "shapefile/"

# Load data for endemics - phylogenetic
# load("MRD_ps_endemics.RData")
# ps_bl_a_end <- ps_bl.a_a_a
# ps_bl_b_end <- ps_bl.b_a_a
# ps_no_a_end <- ps_no.a_a_a
# ps_no_b_end <- ps_no.b_a_a
# rm(MRD.a_a_a, MRD.b_a_a, ps_bl.a_a_a, ps_bl.b_a_a, ps_no.a_a_a, ps_no.b_a_a, comm)

# Load data for all species - phylogenetic    
load(paste0(data_location, "MRD_ps.RData"))
ps_bl_a_all <- ps_bl.a_a_a
ps_bl_b_all <- ps_bl.b_a_a
ps_no_a_all <- ps_no.a_a_a
ps_no_b_all <- ps_no.b_a_a
rm(MRD.a_a_a, MRD.b_a_a, ps_bl.a_a_a, ps_bl.b_a_a, ps_no.a_a_a, ps_no.b_a_a)

# keep community matrix (comm) and phylosim matrices: ps_bl.a_a_a = phylosim branch length, no= without branch length 

# Load data for all species - non-phylogenetic (non-phylogenetic Simpson dissimilarity based on a community matrix)
load(paste0(data_location, "bs_spp.RData"))

# Clean data for clustering

# excluding BOU (which has no species) to avoid NAs in dissimilarity matrix
for(i in 1:100){
  ps_bl_a_all[[i]] <- ps_bl_a_all[[i]][rownames(ps_bl_a_all[[i]])!="BOU",colnames(ps_bl_a_all[[i]])!="BOU"]
  ps_bl_b_all[[i]] <- ps_bl_b_all[[i]][rownames(ps_bl_b_all[[i]])!="BOU",colnames(ps_bl_b_all[[i]])!="BOU"]
  ps_no_a_all[[i]] <- ps_no_a_all[[i]][rownames(ps_no_a_all[[i]])!="BOU",colnames(ps_no_a_all[[i]])!="BOU"]
  ps_no_b_all[[i]] <- ps_no_b_all[[i]][rownames(ps_no_b_all[[i]])!="BOU",colnames(ps_no_b_all[[i]])!="BOU"]
  bs_spp <- bs_spp[rownames(bs_spp)!="BOU",colnames(bs_spp)!="BOU"]
}

# endemics require more extensive exclusion
exclude <- colnames(ps_no_a_end[[1]])[is.na(ps_no_a_end[[1]][2,])]
exclude <- exclude[exclude != "AFG"]

for(i in 1:100){
  ps_bl_a_end[[i]] <- ps_bl_a_end[[i]][!rownames(ps_bl_a_end[[i]]) %in% exclude, !colnames(ps_bl_a_end[[i]]) %in% exclude]
  ps_bl_b_end[[i]] <- ps_bl_b_end[[i]][!rownames(ps_bl_b_end[[i]]) %in% exclude, !colnames(ps_bl_b_end[[i]]) %in% exclude]
  ps_no_a_end[[i]] <- ps_no_a_end[[i]][!rownames(ps_no_a_end[[i]]) %in% exclude, !colnames(ps_no_a_end[[i]]) %in% exclude]
  ps_no_b_end[[i]] <- ps_no_b_end[[i]][!rownames(ps_no_b_end[[i]]) %in% exclude, !colnames(ps_no_b_end[[i]]) %in% exclude]
}

rm(exclude, i)

# Import shapefile for mapping

# load TDWG shapefile
shape <- readOGR(dsn = shapefile_location, layer = "level3")
shape <- shape[!shape@data$LEVEL_3_CO=="BOU",]

# project Behrmann
shape_behr <- spTransform(shape, CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs"))
  shape_behr@data$area=area(shape_behr)

# islands smaller than Sri Lanka - for visualisation
area_threshold <- shape_behr@data[shape_behr@data[,"LEVEL_3_CO"]=="SRL","area"]

shapeCentroids <- gCentroid(shape_behr, byid=TRUE)

shapeBuffer <- gBuffer(shapeCentroids, byid = TRUE, width = 250000)

save.image(paste0(data_location, "data_to_analysis.RData"))
