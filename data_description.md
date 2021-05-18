# Data description

## Primary data

File                        | Description
---------------------------- | --------------------------------------------------------------------------
checklist_names.txt         | Text file containing World Checklist of Vascular Plants taxonomy, taxon status, unique ID and synonyms
checklist_distribution.txt  | Text file containing World Checklist of Vascular Plants taxon presence data in botanical countries
ALLMB.tre		                | Phylogney with GenBank and Open Tree of Life taxa with a backbone provided by Magallón et al. (2015). From Smith & Brown 2018
GBMB.tre			              | Phylogney with GenBank taxa with a backbone provided by Magallón et al. (2015). From Smith & Brown 2018. Used to identify taxa with molecular information
ott.rds                     | Open Tree of Life taxonomy file (.rds version of taxonomy.tsv, Version: 3.0. Available on https://tree.opentreeoflife.org/about/taxonomy-version/ott3.0). Used to identify source for tip label name (NCBI or GBIF)
fin_species_match_NCBI.rds  | List with taxonomy matching results for NCBI phylogeny tip labels with WCVP IDs
gbif_all.rds                | Download of taxonomic information from GBIF for GBIF phylogeny tip labels 
bryophyta.csv               | List of mosses and liverworts families (Bryophyes), available on theplantlist.org/browse
fern_list.txt               | List of fern and fern allies families (Pteridophytes), available on theplantlist.org/browse
soil_raster_layer_000832.tif| Soil raster layer with most probable soil types, available on https://files.isric.org/soilgrids/latest/data/wrb/MostProbable.vrt and saved as geotiff with lower resolution for reasonable computation time
shapefile_biomes/wwf_terr_ecos.*	| Shapefile with terrestrial ecoregions of the world (biomes). Available on https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world
shapefile_bot_countries/level3.*   | Shapefile with polygons for all 369 TDWG level 3 units

 
## Primary data not included but online:
File                            | Description
---------------------------- | --------------------------------------------------------------------------
wc2.1_30s_elev.tif		          | WorldClim elevation data, available on https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_elev.zip
cru_ts4.04.1901.2019.pet.dat.nc | CRU TS v. 4.04, available on https://catalogue.ceda.ac.uk/uuid/89e1e34ec3554dc98594a5732622bce9.  Free account required for download.
cru_ts4.04.1901.2019.tmp.dat.nc | see above 
cru_ts4.04.1901.2019.pre.dat.nc | see above





## Processed data
Includes intermediate results built during analysis, allowing processing of later steps without starting from scretch.

File                        | Description
---------------------------- | --------------------------------------------------------------------------
apg_wcp_jun_20.rds          | .rds file version of checklist_names.txt file for faster loading 
apg_wcp_jun_20_clean.rds	  | World Checklist of Vascular Plants taxonomy working file. Cleaned and family names adjusted for the latest accepted version based on the APG IV (2016) system
SB_tip_labels.rds		        | Non-NCBI based tip labels in ALLMB.tre
fin_species_match_GBIF.rds	|  List with taxonomy matching results for GBIF phylogeny tip labels with WCVP IDs
fin.rds			                | Tip labels in ALLMB.tre
allmb_matched_clean.tre	    | Phylogeny with original tip labels replaced with WCVP unique IDs
level3_mod.shp              | Shapefile with polygons for all 369 TDWG level 3 units and buffer around small islands
allmb_matched_added_species_clean.tre | Phylogeny used in analysis with WCVP tip labels and taxonomically added species
comm_April2021.rds		      | Community matrix containing World Checklist of Vascular Plants taxon presence data in botanical countries, using WCVP unique IDs
polytomie_RD_results_Apr21.rds | Mean root distances for each taxon
mrd_Apr2021.rds		          | Mean root distances for each TDWG level 3 unit
soil.rds 			              | Number of soil types per TDWG level 3 unit
topography.rds			        | Terrain ruggedness index for each TDWG level 3 unit
biomes_olson.rds		        | Biome types for each TDWG level 3 unit
climate.rds			            | Mean and standard deviations of climatic variables for each TDWG level 3 unit
cru/.*                      | Processed CRU climate data
shp_object_fin_analysis.RDS	| Assembled dataset incl. species richness, mean root distance, environmental variables
sem_input_data.rds		      | Final dataset without incomplete cases
mod_selection_Apr2021.RData	| A list called "temp" with model selection results and stats from each 1046429 model runs 
all_models.RData		        | Processed model selection results and best models
best_model.RData		        | final structural equation model lavaan object
