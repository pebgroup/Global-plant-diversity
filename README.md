# Drivers of global plant diversity

## About

We integrate a recent compilation of global species richness on the botanical
country level with novel estimates of diversification rate and environmental
variables, to analyse how their interactions shape today’s global species
richness distribution. All scripts to run this analysis are stored in this
repository.

### Reproducibility and data
We provide here a final dataset (_processed_data/shp_object_fin_analysis.RDS_)
which contains measures on the botanical country level for species richness,
diversification rates and environmental variables, as well as scripts to process
figures, the SEM and to build this data set.

While all files used and produced throughout each step are described here, the
source data needs to be downloaded (e.g. CRU TS, wordclim, soil data, WCVP) from
the respective sources. Version number and information how to obtain data is
listed in the description part below.

## Scripts

Scripts are accessed through *master_analysis.R*, which sources all other
scripts. The header in the master script contains system requirements and
execution time information.

Script                          | Job
--------------------------------|--------------------------------------------------------------------------
0_functions.R                     | Functions called mostly in phylogeny processing and SEM analysis and spatial autocorrelation correction, functions for cleaner structural equation model plots
1_APG_family_lookup.R             | Adjust family names for APG 4 System [1]
2_WCP_cleanup.R                   | Data cleaning
3_match_data.R                    | Get tip label name sources, match NCBI-based tip labels with WCVP accepted taxa IDs
4_common_format_creator_SB.R      | Builds common format for taxonomic names from GBIF data to feed into the taxonomy matcher
5_taxonomic_matcher.v.1.4.R       | Match taxonomic names with WCVP accepted IDs
6_match_data_part2.R              | Match GBIF-based tip labels with WCVP accepted taxa IDs. Replace tip labels with WCVP accepted IDs
7_add_species_server_version.R    | Adding species taxonomically to a phylogeny
8_resolve_polytomies.R            | Repeatedly resolve polytomies in phylogeny and record average root distance for each tip
9_process_geography.R             | Building presence matrix for WCVP accepted taxa in TDWG level 3 units
10_phylostruct.R                   | Computing mean root distance for TDWG level 3 units
11_get_climate.R                   | Process climate data
12_get_soil.R                      | Process soil data
13_get_topography.R                | Process elevation data for terrain ruggedness index
14_get_biomes.R                    | Process biome data
15_assemble_dataset.R              | Assemble final data set (environmental, evolutionary and species richness data)
16_data_prep_and_checks.R          | Get correlations, variable importance analysis, distribution + transformation where necessary, multicolinearity tests, scaling data
17_mod_selection_SEM.R             | Building Structural Equation Models with possible variable combinations, collecting models stats
18_mod_selection_SEM_analysis.R    | Get stats from model selection, manual model adjustments
19_sem_results_further_analysis.R  | Loads the final SEM, spatial autocorrelation analysis, figures and maps 


[1] APG IV. An update of the Angiosperm Phylogeny Group classification for the
orders and families of flowering plants: APG IV. Botanical Journal of the
Linnean Society, 2016, 181, 1–20.


## Data
Primary data includes files obtained online / from collaborators. Processed data
are files which are produced during the analysis steps and stored to shorten
execution time for following analysis steps.

Data set in the process of publication and therefore not provided as raw data:
_checklist_distribution.txt_


### Primary data

File                        | Description
---------------------------- | --------------------------------------------------------------------------
checklist_names.txt         | Text file containing World Checklist of Vascular Plants taxonomy, taxon status, unique ID and synonyms
checklist_distribution.txt  | Text file containing World Checklist of Vascular Plants taxon presence data in botanical countries, provided by collaborators (publication in prep)
ALLMB.tre		                | Phylogeny with GenBank and Open Tree of Life taxa with a backbone provided by Magallón et al. (2015). From Smith & Brown 2018
GBMB.tre			              | Phylogeny with GenBank taxa with a backbone provided by Magallón et al. (2015). From Smith & Brown 2018. Used to identify taxa with molecular information
ott.rds                     | Open Tree of Life taxonomy file (.rds version of taxonomy.tsv, Version: 3.0. Available on https://tree.opentreeoflife.org/about/taxonomy-version/ott3.0). Used to identify source for tip label name (NCBI or GBIF)
fin_species_match_NCBI.rds  | List with taxonomy matching results for NCBI phylogeny tip labels with WCVP IDs
gbif_all.rds                | Download of taxonomic information from GBIF for GBIF phylogeny tip labels 
bryophyta.csv               | List of mosses and liverworts families (Bryophyes), available on theplantlist.org/browse
fern_list.txt               | List of fern and fern allies families (Pteridophytes), available on theplantlist.org/browse
soil_raster_layer_000832.tif| Soil raster layer with most probable soil types, available on https://files.isric.org/soilgrids/latest/data/wrb/MostProbable.vrt and saved as geotiff with lower resolution for reasonable computation time
shapefile_biomes/wwf_terr_ecos.*	| Shapefile with terrestrial ecoregions of the world (biomes). Available on https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world
shapefile_bot_countries/level3.*   | Shapefile with polygons for all 369 TDWG level 3 units
wc2.1_30s_elev.tif		          | WorldClim elevation data, available on https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_elev.zip
cru_ts4.04.1901.2019.pet.dat.nc | CRU TS v. 4.04, available on https://catalogue.ceda.ac.uk/uuid/89e1e34ec3554dc98594a5732622bce9.  Free account required for download.
cru_ts4.04.1901.2019.tmp.dat.nc | see above 
cru_ts4.04.1901.2019.pre.dat.nc | see above





### Processed data
Files built during analysis.

File                        | Description
---------------------------- | --------------------------------------------------------------------------
apg_wcp_jun_21.rds          | .rds file version of checklist_names.txt file for faster loading 
apg_wcp_jun_21_clean.rds	  | World Checklist of Vascular Plants taxonomy working file. Cleaned and family names adjusted for the latest accepted version based on the APG IV (2016) system
SB_tip_labels.rds		        | Non-NCBI based tip labels in ALLMB.tre
fin_species_match_GBIF.rds	|  List with taxonomy matching results for GBIF phylogeny tip labels with WCVP IDs
fin.rds			                | Tip labels in ALLMB.tre
allmb_matched_clean.tre	    | Phylogeny with original tip labels replaced with WCVP unique IDs
level3_mod.shp              | Shapefile with polygons for all 369 TDWG level 3 units and buffer around small islands
allmb_matched_added_species_clean.tre | Phylogeny used in analysis with WCVP tip labels and taxonomically added species
comm_April2021.rds		      | Community matrix containing World Checklist of Vascular Plants taxon presence data in botanical countries, using WCVP unique IDs
polytomie_RD_results_Sep21.rds | Mean root distances for each taxon
mrd_Sep2021.rds		          | Mean root distances for each TDWG level 3 unit
soil.rds 			              | Number of soil types per TDWG level 3 unit
topography.rds			        | Terrain ruggedness index for each TDWG level 3 unit
biomes_olson.rds		        | Biome types for each TDWG level 3 unit
climate.rds			            | Mean and standard deviations of climatic variables for each TDWG level 3 unit
cru/.*                      | Processed CRU climate data
shp_object_fin_analysis.RDS	| Assembled dataset incl. species richness, mean root distance, environmental variables
sem_input_data.rds		      | Final dataset without incomplete cases
mod_selection_Sep2021.RData	| A list called "temp" with model selection results and stats from each 1046429 model runs 
all_models.RData		        | Processed model selection results and best models
best_model.RData		        | final structural equation model lavaan object


## Analysis workflow
<p align="left">
<img src="flow_chart_wcvp.png" width=1200/>  
</p>  