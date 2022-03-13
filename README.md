# Drivers of global plant diversity

## About
Code and data description to reproduce Tietje et al. (in prep) *Drivers of global plant diversity*

### Reproducibility and data
We provide here a final dataset (_processed_data/shp_object_fin_analysis.RDS_)
which contains measures on the botanical country level for species richness,
diversification rates and environmental variables, as well as scripts to process
figures, the SEM and to build this data set.

All files used and produced throughout each step are described here, however due
to the large size and data ownership requirements, the source data needs to be
downloaded (e.g. CRU TS, wordclim, soil data, WCVP) from the respective sources.
Version number and information how to obtain data is listed in the description
part below.

## Scripts
Scripts are described in *master_analysis.R*, and also source some scripts. More
information about system requirements and descriptions of each script can be
found there.

## Data
Primary data includes files obtained online / from collaborators. Processed data
is the final datafile used for analysis.


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
fern_list.txt               | List of fern and fern allies families (Pteridophytes), available on theplantlist.org/browse
soil_raster_layer_000832.tif| Soil raster layer with most probable soil types, available on https://files.isric.org/soilgrids/latest/data/wrb/MostProbable.vrt and saved as geotiff with lower resolution for reasonable computation time
shapefile_biomes/wwf_terr_ecos.*	| Shapefile with terrestrial ecoregions of the world (biomes). Available on https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world
shapefile_bot_countries/level3.*   | Shapefile with polygons for all 369 TDWG level 3 units
wc2.1_30s_elev.tif		          | WorldClim elevation data, available on https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_30s_elev.zip
cru_ts4.04.1901.2019.pet.dat.nc | CRU TS v. 4.04, available on https://catalogue.ceda.ac.uk/uuid/89e1e34ec3554dc98594a5732622bce9.  Free account required for download.
cru_ts4.04.1901.2019.tmp.dat.nc | See above 
cru_ts4.04.1901.2019.pre.dat.nc | See above
tczti_mat_ann_fsy1.nc           | Miocene mean annual temperature from https://www.paleo.bristol.ac.uk/ummodel/users/Bradshaw_et_al_2012/new2/
tczti_precip_ann_fsy1.nc        | Miocene annual precipitation, source as above
tcztk-tczti_mat_ann_fsy1.nc     | Miocene mean annual temperature anomaly
tcztk-tczti_precip_ann_fsy1.nc  | Miocene annual precipitation anomaly
CHELSA_cur_V1_2B_r10m/10min/bio_1.tif | CHELSA Last Glacial Maxiumum Bioclims, Version 1.2, 5/29/2018, obrtained via paleoclim.org
CHELSA_cur_V1_2B_r10m/10min/bio_12.tif | see above
chelsa_LGM_v1_2B_r10m/10min/bio_1.tif |see above
chelsa_LGM_v1_2B_r10m/10min/bio_12.tif |see above
Islands_TDWG_AllData.txt        | Island classifications from http://dx.doi.org/10.1111/j.1466-8238.2011.00728.x




### Final data
Files to make figures and inspect final SEM (scripts 17-19). 

File                        | Description
---------------------------- | --------------------------------------------------------------------------
sem_input_data.rds          | Data.frame to feed to the Structural Equation Model selection process
shp_object_fin_analysis.RDS | Spatial object including all data in sem_input_data.rds
best_model.RData            | Best SEM lavaan model object



## Analysis workflow
<p align="left">
<img src="flow_chart_wcvp.png" width=1200/>  
</p>  
