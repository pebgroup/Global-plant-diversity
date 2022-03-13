# process geography

wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str())) 
library(ape)
library(rgdal)
library(data.table)

phylo <- read.tree("../processed_data/allmb_matched_added_species_Sep21_clean.tre") 
wcp <- readRDS("../processed_data/apg_wcp_jul_21_clean.rds")
shape <- readOGR("../data/shapefile_bot_countries/level3.shp")
dist <- fread("../data/WCVP/world_checklist_names_and_distribution_JUL_21/checklist_distribution.txt")


dist = dist[dist$introduced==0,] #remove introduced and extinct records
dist$area_code_l3 = toupper(dist$area_code_l3)

# process only species that are included in the phylogeny
all(phylo$tip.label %in% wcp$accepted_plant_name_id)
sub <- wcp[wcp$accepted_plant_name_id %in% phylo$tip.label,]
all(phylo$tip.label %in% sub$accepted_plant_name_id)

# create empty community matrix
## rows = botanical countries
## cols = number of unique accepted plant IDs = number of taxa in the phylogeny

comm = matrix(nrow=nrow(shape@data), ncol = length(unique(sub$accepted_plant_name_id)))
rownames(comm) = as.vector(shape@data$LEVEL_3_CO)
rm(phylo, shape)
#if("genus_hybrid_marker" %in% names (wcp)){colnames(comm) = goodspp$checklist_id}else{colnames(comm) = goodspp$plant_name_id} 
colnames(comm) = unique(sub$accepted_plant_name_id)
comm[is.na(comm)] <- 0

# remove distribution records that cannot be assigned to proper TDWG units
dist = dist[dist$area_code_l3 %in% rownames(comm),]

# synonymize distributions table by adding accepted_plant_id column
dist <- merge(dist, wcp[,c("plant_name_id", "accepted_plant_name_id")], all.x=TRUE)
table(is.na(dist$accepted_plant_name_id)) # distribution occurrences with and without accepted plant name ID

# remove distribution records of all taxa that 
dist <- dist[dist$accepted_plant_name_id %in% colnames(comm),]
all(dist$accepted_plant_name_id %in% colnames(comm))

# indexing rownames for faster processing
idx = 1:ncol(comm)
names(idx) = colnames(comm)
dist$accepted_name_idx = idx[dist$accepted_plant_name_id]
rm(idx)


# populating community table
Sys.time()
for(i in 1:nrow(dist)){
  comm[dist[i,"area_code_l3"][[1]], dist[i,"accepted_name_idx"][[1]]] = 1
  if(!i%%1000)cat(i,"\r")
}
rm(i)
Sys.time()

saveRDS(comm, file="../processed_data/comm_September2021.rds")




