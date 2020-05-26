setwd("~/Documents/WOLF/PROJECTS/58 World Checklist paper/analyses 2019")

library(rgdal)

#read WCSP
published <- fread("data/published_names_19_10_2018.csv", header=TRUE, sep="|")
unpublished <- fread("data/unpublished_names_19_10_2018.csv", header=TRUE, sep="|")
wcsp <- rbind(published,unpublished)
rm(published, unpublished)
#wcsp <- readRDS("data/WCSP.apg.rds")
#wcsp <- readRDS("data/wcp_dec_19.rds")

# SWITCH

# 5558 IDs from old version not in new version anymore, 321 IDs were accepted names. Check compare_checklist_versions.R for details
# stelis atra var. atra --> accepted in wcsp.old, not in wcsp new
# Carex longi wcs-228268 acepted ssp in wcsp old, not in wcsp new

#read distributions 
fread("data/published_distribution_19_10_2018.csv", header=TRUE, sep="|") -> published
fread("data/unpublished_distribution_19_10_2018.csv", header=TRUE, sep="|") -> unpublished
dist <- rbind(published,unpublished)
rm(published, unpublished)
dist = dist[dist$introduced==0 & dist$extinct==0,] #remove introduced and extinct records
dist$area_code_l3 = toupper(dist$area_code_l3)

shape = readOGR(dsn = "shapefile", layer = "level3")

goodspp = wcsp[wcsp$genus_hybrid_marker == "" & wcsp$species_hybrid_marker == "" & wcsp$infraspecific_rank == "" & wcsp$species != "" & wcsp$taxon_status_description == "Accepted",]

# create empty community matrix
comm = matrix(nrow=nrow(shape@data), ncol = nrow(goodspp))
rownames(comm) = as.vector(shape@data$LEVEL_3_CO)
colnames(comm) = goodspp$checklist_id
comm[is.na(comm)] <- 0
rm(goodspp)

# remove distribution records that cannot be assigned to proper TDWG units
dist = dist[dist$area_code_l3 %in% rownames(comm),]

# synonymise distributions table by adding accepted_plant_id column
accids = wcsp$accepted_name_id
names(accids) = wcsp$checklist_id
dist$accepted_name_id = accids[as.vector(dist$checklist_id)]
dist$accepted_name_id = as.vector(dist$accepted_name_id)
rm(accids)

# remove distribution records that have no accepted id
dist = dist[dist$accepted_name_id != "",]

# remove distribution records of infraspecific taxa
dist = dist[dist$accepted_name_id %in% colnames(comm),]

# indexing rownames for faster processing
idx = 1:ncol(comm)
names(idx) = colnames(comm)
dist$accepted_name_idx = idx[dist$accepted_name_id]
rm(idx)

# populating community table
dist <- as.data.frame(dist)
for(i in 1:nrow(dist)){
  comm[dist[i,"area_code_l3"], dist[i,"accepted_name_idx"]] = 1
  if(i%%16400==0) print("*")
}
rm(i)



saveRDS(comm, "comm_old_version.rds")
save(comm, file="comm_old_version.RData")
