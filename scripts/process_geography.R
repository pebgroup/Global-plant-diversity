
library(rgdal)
library(data.table)


#read WCSP
#published <- fread("data/published_names_19_10_2018.csv", header=TRUE, sep="|")
#unpublished <- fread("data/unpublished_names_19_10_2018.csv", header=TRUE, sep="|")
#wcsp <- rbind(published,unpublished)
#rm(published, unpublished)
wcsp <- readRDS("data/WCSP_clean.apg.rds")
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
#dist <- readRDS("data/dist.rds")
dist = dist[dist$introduced==0 & dist$extinct==0,] #remove introduced and extinct records
dist$area_code_l3 = toupper(dist$area_code_l3)

shape = readOGR(dsn = "shapefile", layer = "level3")

if("genus_hybrid_marker" %in% names (wcsp)){
  goodspp = wcsp[wcsp$genus_hybrid_marker == "" & 
                   wcsp$species_hybrid_marker == "" & 
                   wcsp$infraspecific_rank == "" & 
                   wcsp$species != "" & 
                   wcsp$taxon_status_description == "Accepted",]
}else{
  goodspp = wcsp[wcsp$genus_hybrid == "" & 
                   wcsp$species_hybrid == "" & 
                   wcsp$infraspecific_rank == "" & 
                   wcsp$species != "" & 
                   wcsp$taxon_status == "Accepted",]
  
}






# compatibility old and new checklist ###################################################
# distribution data is older than the checklist version, check names

# remove non matching species from distribution data
## create old ID format
id_split <- strsplit(as.character(dist$checklist_id), "-")
source <- sapply(id_split, "[[", 1)
source[source=="atoz"] <- "az"
inverse_ID <- paste(sapply(id_split, "[[", 2), source, sep="-")
dist$checklist_id_inverse <- inverse_ID

# check excluded species
length(which(!dist$checklist_id_inverse %in% wcsp$plant_name_id))  # 2836
length(which(!dist$checklist_id %in% wcsp$checklist_id)) # 33
sort(table(dist$area_code_l3[which(!dist$checklist_id_inverse %in% wcsp$plant_name_id)]))

if("genus_hybrid_marker" %in% names (wcsp)){
  dist <- dist[-which(!dist$checklist_id %in% wcsp$checklist_id),]
}else{
  dist <- dist[-which(!dist$checklist_id_inverse %in% wcsp$plant_name_id),]
}

rm(list = c("id_split", "inverse_ID", "source"))


#######################################################################################











# create empty community matrix
comm = matrix(nrow=nrow(shape@data), ncol = nrow(goodspp))
rownames(comm) = as.vector(shape@data$LEVEL_3_CO)
if("genus_hybrid_marker" %in% names (wcsp)){colnames(comm) = goodspp$checklist_id}else{colnames(comm) = goodspp$plant_name_id} 
comm[is.na(comm)] <- 0
rm(goodspp)

# remove distribution records that cannot be assigned to proper TDWG units
dist = dist[dist$area_code_l3 %in% rownames(comm),]

# synonymise distributions table by adding accepted_plant_id column
if("genus_hybrid_marker" %in% names (wcsp)){
  accids = wcsp$accepted_name_id
  names(accids) = wcsp$checklist_id
  dist$accepted_name_id = accids[as.vector(dist$checklist_id)]
  dist$accepted_name_id = as.vector(dist$accepted_name_id)
  table(dist$accepted_name_id=="") # 10 586 occurrence records have no accepted plant id
}else{
  accids = wcsp$accepted_plant_name_id # accepted names
  names(accids) = wcsp$plant_name_id # synonyms etc are name of each entry
  dist$accepted_name_id = accids[as.vector(dist$checklist_id_inverse)]
  dist$accepted_name_id = as.vector(dist$accepted_name_id)
  table(dist$accepted_name_id=="") # 10 170 occurrence records have no accepted plant id
}

rm(accids)



# remove distribution records that have no accepted id
if(length(which(dist$accepted_name_id == "")) > 0){
  dist = dist[dist$accepted_name_id != "",]
}

# remove distribution records of infraspecific taxa
nrow(dist)
dist = dist[dist$accepted_name_id %in% colnames(comm),]
nrow(dist)

# indexing rownames for faster processing
idx = 1:ncol(comm)
names(idx) = colnames(comm)
dist$accepted_name_idx = idx[dist$accepted_name_id]
rm(idx)

# populating community table
Sys.time()
for(i in 1:nrow(dist)){
  comm[dist[i,"area_code_l3"][[1]], dist[i,"accepted_name_idx"][[1]]] = 1
  if(!i%%1000)cat(i,"\r")
}
rm(i)
Sys.time()


saveRDS(comm, file="data/comm.rds")
#saveRDS(comm, file="data/comm_old_wcsp.rds")
















