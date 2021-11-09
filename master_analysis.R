# September 28, 2021
# Written by Melanie Tietje, contributions by Wolf L. Eiserhardt
#  
# 
# This script sources all analysis scripts in the right order. 
# 
# Please note that some scripts may have a long runtime (hours) and require more
# computational resources than the average computer (e.g. RAM and multicore
# parallel processing). To save time, core results of steps requiring a long
# runtime are provided to allow to repeat following steps without further delay.
# Scripts with long runtime and/or requiring multiple cores:
# - add_species_server_version.R
# - resolve_polytomies.R
# - get_climate.R
# - get_soil.R
# - mod_selection_SEM.R
# 
# Scripts 1-18 are data compilation, prep and model selection, for working with
# the final data see 19_sem_results_further_analysis.R
#
#
# Tested with R version 3.6.3, in RStudio Version 1.3.1073
# 
# All used R packages are cited in the supplement




# Setup -------------------------------------------------------------------

ncores = 16

dir.create("figures")
dir.create("processed_data")

source("scripts/0_functions.R")
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
library(rgbif)
library(writexl) # to export model results as table 



# WCVP PREP ###############################################################

# Process WCVP data -------------------------------------------------------
# clear workspace, keep functions
rm(list = setdiff(ls(), lsf.str())) 

wcp <- fread("data/WCVP/world_checklist_names_and_distribution_JUL_21/checklist_names.txt",
             quote="", na.strings = "")
saveRDS(wcp, "processed_data/wcp_jul_21.rds")



# Fix family names following APG IV ---------------------------------------
db <- "WCP" 
data_folder_path <- "./processed_data/"
WCP_input_filename <- "wcp_jul_21.rds"
source("1_APG_family_lookup.R")



# Clean WCVP --------------------------------------------------------------
# clear workspace, keep functions
rm(list = setdiff(ls(), lsf.str())) 

wcp <- readRDS("processed_data/apg_wcp_jul_21.rds")
source("2_WCP_cleanup.R")
saveRDS(wcp2, "processed_data/apg_wcp_jul_21_clean.rds")





# PHYLOGENY  ##############################################################


## 1) match tip labels Smith & Brown phylogeny with WCVP names
## 2) add missing species taxonomically
## 3) calculate root distances for all tip labels, then average per TDWG level 3 unit


# Match phyologeny tip labels to WCVP, part 1 -----------------------------
# clear workspace, keep functions
rm(list = setdiff(ls(), lsf.str())) 

tax_level <- "elevated_to_species_id" 

phylo <- read.tree("data/ALLMB.tre") # phylogeny from Smith&Brown 2018
ott <- readRDS("data/ott.rds") # OTL taxonomy to identify tip labels
matches <- readRDS("data/fin_species_match_NCBI_wcvp21.rds") # NCBI-WCSP matches from the taxonomy matcher

source("3_match_data.R")

# export tip labels that have not been matched with NCBI
saveRDS(fin[!is.na(fin$gbif_id),], file="processed_data/SB_tip_labels_21.rds")
# export all tip labels with source
saveRDS(fin, file="processed_data/fin.rds")



# Prep unmatched tip labels for taxonomy matching ------------------------
# clear workspace, keep functions
rm(list = setdiff(ls(), lsf.str())) 

# Load Smith & Brown tip labels not covered by NCBI (-->GBIF)
sb <- readRDS("processed_data/SB_tip_labels_21.rds")


# GBIF download to complete tax info --------------------------------------
rm(list = setdiff(ls(), lsf.str())) 
sb <- readRDS("SB_tip_labels_21.rds")
ids <- sb$gbif_id
res <- list()
for(i in 1:length(ids)){
  tryCatch({q
    res[[i]] <- name_usage(ids[i], return="data")
    if(!i%%10)cat(i,"\r")
    # failsafe
    if(i %in% c(50000, 100000, 150000, 200000)){saveRDS(res, "gbif.rds")}
    if(i==length(ids)){saveRDS(res, "gbif_SB_tip_labels_21.rds")}
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



# Prep tip labels for taxonomy matching -----------------------------------
rm(list = setdiff(ls(), lsf.str())) 
gbif <- readRDS("data/gbif_SB_tip_labels_21.rds")
source("4_common_format_creator_SB.R")
saveRDS(input, "processed_data/input_tip_labels.rds")



# Fix family names following APG IV ---------------------------------------
db <- "GBIF" 
data_folder_path <- "./processed_data/"
gbif_input_filename <- "input_tip_labels.rds"
source("1_APG_family_lookup.R")
# writes apg_input_tip_labels.rds



# Run taxonomy matchig for GBIF tip labels --------------------------------
# clear workspace, keep functions
rm(list = setdiff(ls(), lsf.str())) 

DB.name <- "GBIF"
data_folder_path <- "./processed_data/" 
results_folder_path <- "./processed_data/"
wcvp_input_filename <- "apg_wcp_jul_21_clean.rds"
gbif_input_filename <- "apg_input_tip_labels.rds"
fin_file_species <- paste0("fin_species_match_", DB.name, ".rds")
bry <- read.csv("data/bryophyta.csv")

source("5_taxonomic_matcher.v.1.4.R")
# output: fin_species_match_GBIF.rds


# Match phyologeny tip labels to WCVP, part 2 -----------------------------
rm(list = setdiff(ls(), lsf.str())) 

matches <- readRDS("processed_data/fin_species_match_GBIF.rds")
fin <- readRDS("processed_data/fin.rds")
phylo <- read.tree("data/ALLMB.tre") # phylogeny
tax_level <- "elevated_to_species_id" 
phylogb_a <- read.tree("data/GBMB.tre") # genbank tree to identify tips with molecular data
source("6_match_data_part2.R")
write.tree(clean_tree, "processed_data/allmb_matched_clean.tre")



# Adding species ----------------------------------------------------------
## LONG RUNTIME ! ca 11 hours

rm(list = setdiff(ls(), lsf.str())) 

phylo <- read.tree("processed_data/allmb_matched_clean.tre")
wcp <- readRDS("processed_data/apg_wcp_jul_21_clean.rds")
ferns <- as.vector(read.csv("data/fern_list.txt")[,1])

source("7_add_species_server_version.R", print.eval=T)
# writes allmb_matched_added_species_Sep21_clean.tre



# Resolve polytomies, get root distances -----------------------------------
## LONG RUNTIME ! (ca 1h on 16 cores)
## PARALLEL PROCESSING

rm(list = setdiff(ls(), lsf.str())) 


n_cores =ncores
reps <- 1000
source("8_resolve_polytomies.R", print.eval = TRUE)
saveRDS(res, "processed_data/polytomie_RD_results.rds")
# writes polytomie_RD_results_Sep21.rds



# GEOGRAPHY ##################################################################



# Build presence matrix ----------------------------------------------------
# clear workspace, keep functions
rm(list = setdiff(ls(), lsf.str())) 


phylo <- read.tree("processed_data/allmb_matched_added_species_Sep21_clean.tre") 
wcp <- readRDS("processed_data/apg_wcp_jul_21_clean.rds")
shape <- readOGR("data/shapefile_bot_countries/level3.shp")
dist <- fread("data/WCVP/world_checklist_names_and_distribution_JUL_21/checklist_distribution.txt")
source("9_process_geography.R", print.eval = TRUE)
saveRDS(comm, file="processed_data/comm_September2021.rds")



# Mean root distance for each TDGW level 3 unit ----------------------------
rm(list = setdiff(ls(), lsf.str())) 


res <- readRDS("processed_data/polytomie_RD_results_Sep21_2.rds")
phylo <- read.tree("processed_data/allmb_matched_added_species_Sep21_clean.tre") 
comm <- readRDS("processed_data/comm_September2021.rds")

# compare phylogenetic with geographic data again
table(colnames(comm) %in% phylo$tip.label) # should be all represented

n_cores=4
source("10_phylostruct.R", print.eval = TRUE)
saveRDS(mrd.res, "processed_data/mrd_Sep2021.rds")






# Environmental variables ####################################################

## LONG RUNTIME ! ca 6 hours

rm(list = setdiff(ls(), lsf.str())) 

# raster files are loaded inside the scripts
dir.create("processed_data/cru")

source("11_get_climate.R")
saveRDS(res, "processed_data/climate.rds")


rm(list = setdiff(ls(), lsf.str())) 

shape <- readOGR("data/shapefile_bot_countries/level3_mod.shp")
soil <- raster("data/soil_raster_layer_000832.tif")
source("12_get_soil.R")
saveRDS(res, file="processed_data/soil.rds")

rm(list = setdiff(ls(), lsf.str())) 
shape <- readOGR("data/shapefile_bot_countries/level3_mod.shp")
elev <- raster("data/wc2.1_30s_elev.tif")
source("13_get_topography.R")
saveRDS(res, file="processed_data/topography.rds")

rm(list = setdiff(ls(), lsf.str()))  

shape <- readOGR("data/shapefile_bot_countries/level3.shp")
olson <- readOGR("data/shapefile_biomes/wwf_terr_ecos.shp")
source("14_get_biomes.R")
saveRDS(res.df, file="processed_data/biomes_olson.rds")



# Assemble final dataset --------------------------------------------------
rm(list = setdiff(ls(), lsf.str()))  

sr <- readRDS("processed_data/comm_September2021.rds")
#sr_2020 <- readRDS("processed_data/comm_April2021.rds")
mrd <- readRDS("processed_data/mrd_Sep2021.rds")
#mrd_2020 <- readRDS("processed_data/mrd_Apr2021.rds")
shape <- readOGR("data/shapefile_bot_countries/level3.shp")
soil <- readRDS("processed_data/soil.rds")
climate <- readRDS("processed_data/climate.rds")
topography <- readRDS("processed_data/topography.rds")
biomes <- readRDS("processed_data/biomes_olson.rds")

source("15_assemble_dataset.R")
saveRDS(shp, "processed_data/shp_object_fin_analysis.RDS")



# Data prep and checks -----------------------------------------------------
# variable importance, correlations, scaling and variable distributions
rm(list = setdiff(ls(), lsf.str()))  

shp <- readRDS("processed_data/shp_object_fin_analysis.RDS")
source("16_data_prep_and_checks.R")
# saves the final dataset with only complete cases: sem_input_data.rds





# STRUCTURAL EQUATION MODELING ###########################################

# SEM model selection -----------------------------------------------------
## LONG RUNTIME ! ca 5 hours on 16 cores
## PARALLEL PROCESSING

rm(list = setdiff(ls(), lsf.str()))  

dat_no.na <- readRDS("data/sem_input_data.rds")

# set max number of variables (to add to the fixed ones) used and register cores
var.max <- 10
registerDoParallel(ncores)

source("17_mod_selection_SEM.R")
#save(temp, file="processed_data/mod_selection_Sep2021.RData")



# Process SEMs ------------------------------------------------------------
# get stats from model selection runs and manually adjust models for better fit
rm(list = setdiff(ls(), lsf.str()))  

load("processed_data/mod_selection_Sep2021.RData")
dat_no.na <- readRDS("processed_data/sem_input_data.rds")
source("18_mod_selection_SEM_analysis.R", print.eval = TRUE)

save.image("processed_data/all_models.RData")
save(max.var.mod_mod7, max.var.mod.fit.mod7, file="processed_data/best_model.RData")




# Maps, Figures, Output --------------------------------------------------
## model parameters, local effects, spatial autocorrelation assessment and
## correction, latitudinal SR and MRD patterns, robustness
rm(list = setdiff(ls(), lsf.str()))  

load("processed_data/best_model.RData")
shp <- readRDS("processed_data/shp_object_fin_analysis.RDS")
dat_no.na <- readRDS("processed_data/sem_input_data.rds")

source("19_sem_results_further_analysis.R")






# FIN #####


