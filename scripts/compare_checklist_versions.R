# compare checklist versions
 ## I don`t correct familiy names and stuff, this is data as downloaded - tax errors are the same in both`

# #wcsp1 <- readRDS("data/wcsp_old.rds")
# #wcsp <- readRDS("data/WCSP.apg.rds")
# dec19 <- readRDS("data/wcp_dec_19.rds")
# jun20 <- readRDS("data/wcp_jun_20.rds")

library(data.table)
#dec19 <- fread("../data/wcs_dec_19/checklist_names.txt", quote = "")
jun20 <- fread("../data/wcs_jun_20/world_checklist_names_and_distribution_jun_2020/checklist_names.txt", quote = "")
jun20_2 <- read.csv("../data/wcs_jun_20/world_checklist_names_and_distribution_jun_2020/checklist_names.txt", quote = "")
wcspnames <- read.csv("../data/wcs_jun_20/world_checklist_names_and_distribution_jun_2020/checklist_names.txt", sep="|")

# accepted angiosperm species according to checklist
sub <- jun20[jun20$taxon_rank == "Species" & jun20$taxon_status == "Accepted",]
length(unique(sub$accepted_plant_name_id))




# alexanders code:
wcsp.names <- fread("../data/wcs_jun_20/world_checklist_names_and_distribution_jun_2020/checklist_names.txt", quote = "")
# Subsetting for only accepted species and removing entries at genus level
wcsp.names.subset <- wcsp.names[grepl("Accepted", wcsp.names$taxon_status),] # Including only accepted species
wcsp.names.subset <- wcsp.names.subset[!(wcsp.names.subset$species == ""), ] # Excluding entries at genus level
length(wcsp.names.subset$accepted_plant_name_id) # I now get just 204284 species - I get 391307

sub <- wcsp.names[wcsp.names$taxon_rank == "Species" & wcsp.names$taxon_status == "Accepted",]
length(unique(sub$accepted_plant_name_id))


###########################

# Ferns + mosses
ferns <- as.vector(read.csv("../data/fern_list.txt")[,1])
jun20.ferns <- jun20[jun20$family %in% ferns,]
length(unique(jun20.ferns$accepted_plant_name_id))
fern_IDs <- unique(jun20.ferns$accepted_plant_name_id)


goodssp2 = dec19[dec19$genus_hybrid == "" & 
                 dec19$species_hybrid == "" & 
                 dec19$infraspecific_rank == "" & 
                 dec19$species != "" & 
                 dec19$taxon_status == "Accepted",]
length(unique(goodssp2$accepted_plant_name_id))
goodssp2 <- goodspp2[!goodspp2$accepted_plant_name_id %in% fern_IDs,]
length(unique(goodssp2$accepted_plant_name_id))

goodssp3 = jun20[jun20$genus_hybrid == "" & 
                   jun20$species_hybrid == "" & 
                   jun20$infraspecific_rank == "" & 
                   jun20$species != "" & 
                   jun20$taxon_status == "Accepted",]
length(unique(goodssp3$accepted_plant_name_id))
goodssp3 <- goodssp3[!goodssp3$accepted_plant_name_id %in% fern_IDs,]
length(unique(goodssp3$accepted_plant_name_id))

dec19.dist <- fread("../data/wcs_dec_19/checklist_distribution.txt")
jun20.dist <- fread("../data/wcs_jun_20/checklist_distribution.txt")

# fixx names in dist file
## wcsp dec19
accids = dec19$accepted_plant_name_id # accepted names
names(accids) = dec19$plant_name_id # synonyms etc are name of each entry
dec19.dist$accepted_name_id = accids[as.vector(dec19.dist$db_id)]
dec19.dist$accepted_name_id = as.vector(dec19.dist$accepted_name_id)
table(dec19.dist$accepted_name_id=="") # 12 948 occurrence records have no accepted plant name ID

## wcsp jun20
accids = jun20$accepted_plant_name_id # accepted names
names(accids) = jun20$plant_name_id # synonyms etc are name of each entry
jun20.dist$accepted_name_id = accids[as.vector(jun20.dist$plant_name_id)]
jun20.dist$accepted_name_id = as.vector(jun20.dist$accepted_name_id)
table(jun20.dist$accepted_name_id=="") # 12 813 occurrence records have no accepted plant name ID


table(dec19.dist$db_id %in% dec19$plant_name_id)
table(jun20.dist$plant_name_id %in% jun20$plant_name_id)

length(unique(dec19.dist$accepted_name_id[dec19.dist$accepted_name_id %in% dec19$accepted_plant_name_id])) 
length(unique(dec19.dist$accepted_name_id[dec19.dist$accepted_name_id %in% goodssp2$accepted_plant_name_id])) 

length(unique(jun20.dist$accepted_name_id[jun20.dist$accepted_name_id %in% jun20$accepted_plant_name_id])) 
length(unique(jun20.dist$accepted_name_id[jun20.dist$accepted_name_id %in% goodssp3$accepted_plant_name_id])) 

# 
# 
# phylo <- read.tree("../trees/allmb_matched_added_species_2.tre")
# length(phylo$tip.label)

# 
# # no hybrids, no infra ranks, no genera, only species, accepted
# length(unique(goodspp2$plant_name_id))
# length(unique(goodspp3$plant_name_id))
# 
# # Compare checklists
# 
# ## create new ID format
# id_split <- strsplit(as.character(wcsp.old$checklist_id), "-")
# source <- sapply(id_split, "[[", 1)
# source[source=="atoz"] <- "az"
# new_ID <- paste(sapply(id_split, "[[", 2), source, sep="-")
# wcsp.old$checklist_id_new <- new_ID
# 
# 
# # Compare new checklist versions
# table(wcsp$plant_name_id %in% wcsp.no$plant_name_id) # adjusted families new to original new: all included
# 
# table(wcsp.old$checklist_id_new %in% wcsp$plant_name_id) # 5558 IDs not included anymore
# table(wcsp$plant_name_id %in% wcsp.old$checklist_id_new) # 17398 new IDs in new version
# 
# # Compare missing
# 
# old_missing <- wcsp.old[which(!wcsp.old$checklist_id_new %in% wcsp$plant_name_id),]
# table(old_missing$taxon_status_description)
# 
# View(old_missing[old_missing$taxon_status_description=="Accepted",])
# 
# 
# 
# 
# 
# ########### compare distribution data interlinked with wcsp #############333
# comm <- readRDS("data/comm.rds")
# comm_old <- readRDS("data/comm_old_wcsp.rds")
# comm_org <- readRDS("comm_old_version.rds")
# 
# comb <- data.frame(region =row.names(comm), 
#                    sr_new=rowSums(comm), 
#                    sr_old = rowSums(comm_old),
#                    sr_org = rowSums(comm_org))
# plot(log(comb$sr_new), log(comb$sr_org))
# 
# 
# 
# 
#      
#      