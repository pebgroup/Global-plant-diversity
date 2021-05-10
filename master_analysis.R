# May 4, 2021
# Written by Melanie Tietje, contributions by Wolf L. Eiserhardt
#  
# 
# This script sources all analysis scripts. 
# 
# Please note that the most time consuming calculations are muted by the following 
# and instead the results are loaded where needed. This affects scripts to
# taxonomically add species to the phylogeny, resolve polytomies, extract climate
# variables and run automated SEM selection procedure. Several analyses are
# set up for parallel processing (up to 16 cores) to speed up performance.
# 
# Processing partly large data files, tested on a linux system (Ubuntu 20.04) with 
# Intel© Core™ i7-8565U CPU @ 1.80GHz × 4 , 16GB RAM, RStudio Version 1.3.1073
# All used R packages are cited in the supplement


# setup  ----------------------------------------------------------------------------------------------
run_p = FALSE
run_s=FALSE
run_e=FALSE
run_ms=FALSE

dir.create("publish_figures")
dir.create("processed_data")

source("publish_scripts/functions.R")

library(ape)
library(caret)
library(castor)
library(cowplot) # for multiplot arrangements
library(data.table)
library(dplyr)
library(foreach)
library(gbm)
library(geosphere)
library(ggcorrplot)
library(ggplot2)
theme_set(theme_bw())
library(ggpubr)
library(ggsci)
library(lavaan)
library(parallel)
library(phytools)
library(plyr)
library(png)
library(raster)
library(rgdal)
library(rgeos)
library(semPlot)
library(sf)
library(Simpsons)
library(spatialEco)
library(spdep)
library(stringr)
library(tidyr)
library(writexl) # to export model results as table 

# library(tidyverse)
# library(scales)
# library(vegan)
# library(VIM)
# library(doParallel)
# library(ggthemes)



# Process WCVP taxonomic reference data  -----------------------------------------------------------
# clear workspace, keep functions and settings
rm(list = setdiff(ls(), c(lsf.str(), "run_p", "run_s", "run_e", "run_ms"))) 


wcp <- fread("data/wcs_jun_20/world_checklist_names_and_distribution_jun_2020/checklist_names.txt", quote="")
saveRDS(wcp, "processed_data/wcp_jun_20.rds")



# APG family lookup ------------------------------------------------------------------------------
db <- "WCP" 
data_folder_path <- "./processed_data/"
WCP_input_filename <- "wcp_jun_20.rds"
source("publish_scripts/APG_family_lookup.R")



# Clean WCP --------------------------------------------------------------------------------------
wcp <- readRDS("processed_data/apg_wcp_jun_20.rds")
source("publish_scripts/WCP_cleanup.R")
saveRDS(wcp2, "processed_data/apg_wcp_jun_20_clean.rds")





# PHYLOGENY  #####################################################################################


## 1) match tip labels Smith & Brown phylogeny with WCVP names
## 2) add missing species taxonomically
## 3) calculate root distances for all tip labels, then average per TDWG level 3 unit

# Get tip label name sources from OTT, match NCBI names, export GBIF names ------------------------
rm(list = setdiff(ls(), c(lsf.str(), "run_p", "run_s", "run_e", "run_ms"))) 
tax_level <- "elevated_to_species_id" 

phylo <- read.tree("trees/ALLMB.tre") # phylogeny
ott <- readRDS("data/ott.rds") # OTL taxonomy to identify tip labels
matches <- readRDS("data/fin_species_match_NCBI.rds") # NCBI-WCSP matches from the taxonomy matcher
source("publish_scripts/match_data.R")
# export tip labels that have not been matched with NCBI
saveRDS(fin[!is.na(fin$gbif_id),], file="processed_data/SB_tip_labels.rds")
saveRDS(fin, file="processed_data/fin.rds")



# Create common format to source into taxonomy matcher --------------------------------------------
rm(list = setdiff(ls(), c(lsf.str(), "run_p", "run_s", "run_e", "run_ms"))) 
# Load Smith & Brown tip labels not covered by NCBI 
sb <- readRDS("processed_data/SB_tip_labels.rds")
source("publish_scripts/common_format_creator_SB.R")
saveRDS(input, "processed_data/input_tip_labels.rds")



# Summon the APG lookup script for tip labels -----------------------------------------------------
db <- "GBIF" 
data_folder_path <- "./processed_data/"
gbif_input_filename <- "input_tip_labels.rds"
source("publish_scripts/APG_family_lookup.R")
# writes apg_input_tip_labels.rds



# Run taxonomy matcher for GBIF tip labels --------------------------------------------------------
rm(list = setdiff(ls(), c(lsf.str(), "run_p", "run_s", "run_e", "run_ms"))) 

## setup
DB.name <- "GBIF"
data_folder_path <- "./processed_data/" 
results_folder_path <- "./processed_data/"
wcvp_input_filename <- "apg_wcp_jun_20_clean.rds"
gbif_input_filename <- "apg_input_tip_labels.rds"
fin_file_species <- paste0("fin_species_match_", DB.name, ".rds")
bry <- read.csv("data/bryophyta.csv")

source("publish_scripts/taxonomic_matcher.v.1.4.R")
# output: fin_species_match_GBIF.rds


# Load results and match remaining tip labels, remove duplicated tip labels ----------------------
rm(list = setdiff(ls(), c(lsf.str(), "run_p", "run_s", "run_e", "run_ms"))) 
matches <- readRDS("processed_data/fin_species_match_GBIF.rds")
fin <- readRDS("processed_data/fin.rds")
phylo <- read.tree("trees/ALLMB.tre") # phylogeny
tax_level <- "elevated_to_species_id" 
phylogb_a <- read.tree("trees/GBMB.tre") # genbank tree to identify tips with molecular data behind them
source("publish_scripts/match_data_part2.R")
write.tree(clean_tree, "processed_data/allmb_matched_clean.tre")



# Add missing species taxonomically --------------------------------------------------------------
rm(list = setdiff(ls(), c(lsf.str(), "run_p", "run_s", "run_e", "run_ms"))) 
phylo <- read.tree("processed_data/allmb_matched_clean.tre")
wcp <- readRDS("processed_data/apg_wcp_jun_20_clean.rds")
ferns <- as.vector(read.csv("data/fern_list.txt")[,1])

if(run_s){
  source("publish_scripts/add_species_server_version.R", print.eval=T)
}else{
  fin_tree <- read.tree("processed_data/allmb_matched_added_species_clean.tre") # phylogeny
}



# Resolve polytomies and get root distance for each tip ------------------------------------------
## we recommend running this in parallel on a machine that can handle at least 16 cores
## 16 cores takes about 1h computation time

rm(list = setdiff(ls(), c(lsf.str(), "run_p", "run_s", "run_e", "run_ms"))) 

if(run_p){
  n_cores =16
  reps <- 1000
  source("publish_scripts/resolve_polytomies.R", print.eval = TRUE)
  saveRDS(res, "processed_data/polytomie_RD_results.rds")
}




# GEOGRAPHY #####################################################################################



# Process geography -------------------------------------------------------------------------------
rm(list = setdiff(ls(), c(lsf.str(), "run_p", "run_s", "run_e", "run_ms"))) 

phylo <- read.tree("processed_data/allmb_matched_added_species_clean.tre") 
wcp <- readRDS("processed_data/apg_wcp_jun_20_clean.rds")
shape = readOGR(dsn = "shapefile", layer = "level3")
dist <- fread("data/wcs_jun_20/world_checklist_names_and_distribution_jun_2020/checklist_distribution.txt")
source(process_geography.R, print.eval = TRUE)
saveRDS(comm, file="processed_data/comm_April2021.rds")



# Mean root distance for each TDGW level 3 unit --------------------------------------------------
rm(list = setdiff(ls(), c(lsf.str(), "run_p", "run_s", "run_e", "run_ms"))) 

res <- readRDS("processed_data/polytomie_RD_results_Apr21.rds")
phylo <- read.tree("processed_data/allmb_matched_added_species_clean.tre") 
comm <- readRDS("processed_data/comm_April2021.rds")

# compare phylogenetic with geographic data again
table(colnames(comm) %in% phylo$tip.label) # should be all represented

n_cores=4
source("publish_scripts/phylostruct.R", print.eval = TRUE)
saveRDS(mrd.res, "processed_data/mrd_Apr2021.rds")






# Environmental variables ##########################################################################
# running all scripts takes about 6 hours
if(run_e){
  rm(list = setdiff(ls(), c(lsf.str(), "run_p", "run_s", "run_e", "run_ms"))) 
  #library(cruts)

  # raster files are loaded inside the scripts
  dir.create("processed_data/cru")
  
  source("publish_scripts/get_climate.R")
  saveRDS(res, "processed_data/climate.rds")
  
  rm(list = setdiff(ls(), c(lsf.str(), "run_p", "run_s", "run_e", "run_ms"))) 
  shape <- readOGR("shapefile/level3_mod.shp")
  soil <- raster("data/soilgrids/soil_raster_layer_000832.tif")
  source("publish_scripts/get_soil.R")
  saveRDS(res, file="processed_data/soil.rds")
  
  rm(list = setdiff(ls(), c(lsf.str(), "run_p", "run_s", "run_e", "run_ms"))) 
  shape <- readOGR("shapefile/level3_mod.shp")
  elev <- raster("data/worldclim/wc2.1_30s_elev.tif")
  source("publish_scripts/get_topography.R")
  saveRDS(res, file="processed_data/topography.rds")
  
  rm(list = setdiff(ls(), c(lsf.str(), "run_p", "run_s", "run_e", "run_ms"))) 

  shape <- readOGR("shapefile/level3.shp")
  olson <- readOGR("shapefile/olson/wwf_terr_ecos.shp")
  source("publish_scripts/get_biomes.R")
  saveRDS(res.df, file="processed_data/biomes_olson.rds")
  
}



# Assemble final data set -------------------------------------------------------------------------
rm(list = setdiff(ls(), c(lsf.str(), "run_p", "run_s", "run_e", "run_ms"))) 

sr <- readRDS("processed_data/comm_April2021.rds")
mrd <- readRDS("processed_data/mrd_Apr2021.rds")
shape <- readOGR("shapefile/level3.shp")
soil <- readRDS("processed_data/soil.rds")
climate <- readRDS("processed_data/climate.rds")
topography <- readRDS("processed_data/topography.rds")
biomes <- readRDS("processed_data/biomes_olson.rds")

source("publish_scripts/assemble_dataset.R")
saveRDS(shp, "processed_data/shp_object_fin_analysis.RDS")



# Data prep and checks ---------------------------------------------------------------------------
# variable importance, correlations, scaling and variable distributions
rm(list = setdiff(ls(), c(lsf.str(), "run_p", "run_s", "run_e", "run_ms"))) 

shp <- readRDS("processed_data/shp_object_fin_analysis.RDS")
source("publish_scripts/data_prep_and_checks.R")
# saves internally the final dataset with only complete cases: sem_input_data.rds





# STRUCTURAL EQUATION MODELING ##################################################################

# SEM model selection ----------------------------------------------------------------------------
if(run_ms){
  rm(list = setdiff(ls(), c(lsf.str(), "run_p", "run_s", "run_e", "run_ms"))) 

  dat_no.na <- readRDS("data/sem_input_data.rds")
  
  # set max number of variables (to add to the fixed ones) used and register cores
  var.max <- 10
  registerDoParallel(16)
  
  source("publish_scripts/mod_selection_SEM.R")
  save(temp, file="mod_selection_Apr2021.RData")
}


# Evaluate SEMs ---------------------------------------------------------------------------------------
# get stats from model selection runs and manually adjust models for better fit
rm(list = setdiff(ls(), c(lsf.str(), "run_p", "run_s", "run_e", "run_ms"))) 

load("processed_data/mod_selection_Apr2021.RData")
dat_no.na <- readRDS("processed_data/sem_input_data.rds")
source("publish_scripts/mod_selection_SEM_analysis.R", print.eval = TRUE)

save.image("processed_data/all_models.RData")
save(max.var.mod_mod7, max.var.mod.fit.mod7, file="processed_data/best_model.RData")




# Anaylsis using final model - Maps, Figures, Output --------------------------------------------------
## closer model inspection, local effects, spatial autocorrelation assessment and correction,
## latitudinal SR and MRD patterns, robustness
rm(list = setdiff(ls(), c(lsf.str(), "run_p", "run_s", "run_e", "run_ms"))) 

load("processed_data/best_model.RData")
shp <- readRDS("processed_data/shp_object_fin_analysis.RDS") # load spatial object
dat_no.na <- readRDS("processed_data/sem_input_data.rds")
#shape <- readOGR("shapefile/level3.shp")

source("publish_scripts/sem_results_further_analysis.R")


# FIN #####


