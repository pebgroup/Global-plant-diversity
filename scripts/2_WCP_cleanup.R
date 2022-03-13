# Script to check inconsistencies in WCVP #####


wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str())) 
library(dplyr)

wcp <- readRDS("../processed_data/apg_wcp_jul_21.rds")

# Remove entries that don't have accepted IDs
wcp <- wcp[-which(is.na(wcp$accepted_plant_name_id)),]

# Remove not-accepted taxa pointing to themselves
## This is mostly artificial hybrids - do we want to remove them or keep them?
wcp <- wcp[-which(as.character(wcp$plant_name_id)==as.character(wcp$accepted_plant_name_id) 
                    & wcp$taxon_status!="Accepted"),]

# Create index to identify faulty entries
ind <- which(wcp$taxon_status!="Accepted" & wcp$plant_name_id %in% wcp$accepted_plant_name_id)
wcp.tmp <- wcp[ind,]
table(wcp.tmp$taxon_status) # status of "supposed to be accepted"-taxa

# down the rabbit hole: get accepted IDs from falsely linked as accepted taxa: are those accepted?
lev1 <- wcp[wcp$plant_name_id %in% wcp.tmp$accepted_plant_name_id,]
table(lev1$taxon_status)
# the next level link actually points to accepted. change those entries

# subset wcp.tmp to those species that are accepted in lev 1 -> in those cases, the wcp.tmp accepted ID is the right accepted ID
wcp.tmp.sub <- wcp.tmp[wcp.tmp$accepted_plant_name_id %in% lev1$plant_name_id[lev1$taxon_status=="Accepted"],]

# merging first level: match original accepted with level 1 plant name id to get new accepted ID. The valid accepted is always the last one added
wcp2 <- merge(wcp, wcp.tmp[,c("plant_name_id", "accepted_plant_name_id")], 
              by.x="accepted_plant_name_id", by.y="plant_name_id",
              all.x=TRUE)

# remove already solved from next matching step
lev1.tmp <- lev1[-which(lev1$taxon_status=="Accepted"),]
nrow(lev1.tmp)
# lev1.tmp <- lev1.tmp %>%
#   rename(accepted_plant_name_id_lev2 = accepted_plant_name_id)

# level 2 - get accepted IDs from still falsely linked as accepted taxa
# lev2 <- wcp[wcp$plant_name_id %in% lev1.tmp$accepted_plant_name_id_lev2,]
# table(lev2$taxon_status)
# 2 more accepted species

# subset lev1.tmp to those species that are accepted in lev 2 -> in those cases, the lev1.tmp accepted ID is the right accepted ID
#lev1.tmp.sub <- lev1.tmp[lev1.tmp$accepted_plant_name_id_lev2 %in% lev2$plant_name_id[lev2$taxon_status=="Accepted"],]

# merging second level: the valid accepted is always the last one added
# wcp2 <- merge(wcp2, lev1.tmp.sub[,c("plant_name_id", "accepted_plant_name_id_lev2")], 
#               by.x="accepted_plant_name_id.y", by.y="plant_name_id",
#               all.x=TRUE)

#View(wcp2[,c("accepted_plant_name_id", "accepted_plant_name_id.y", "accepted_plant_name_id_lev2")])


# Remove taxa stuck in a loop of synonyms
#View(wcp2[loops,])

loops <- which(wcp2$plant_name_id %in% lev1$plant_name_id[lev1$taxon_status!="Accepted"])



# Cleanup ####
## list to update already existing trees and files etc.
update_list <- wcp2[!is.na(wcp2$accepted_plant_name_id.y),]
# accepted_plant_name_id.y must replace accepted plant_name_id

#update_list$accepted_plant_name_id.y[!is.na(update_list$accepted_plant_name_id_lev2)] <- 
#  update_list$accepted_plant_name_id_lev2[!is.na(update_list$accepted_plant_name_id_lev2)]
update_list <- update_list %>%
  dplyr::rename(old_accepted_ID = accepted_plant_name_id) %>%
  dplyr::rename(new_accepted_ID = accepted_plant_name_id.y) %>%  
  dplyr::select(c(plant_name_id, old_accepted_ID, new_accepted_ID))


## building a new clean wcp version
# first half: selecting entries with old accepted IDs (in update_list)
# second half: replace with new accepted ID
wcp2$accepted_plant_name_id[!is.na(wcp2$accepted_plant_name_id.y)] <- wcp2$accepted_plant_name_id.y[!is.na(wcp2$accepted_plant_name_id.y)]

# remove the second accepted column
wcp2 <- wcp2 %>%
  dplyr::select(-c(accepted_plant_name_id.y))

table(wcp2$taxon_status[wcp2$plant_name_id %in% wcp2$accepted_plant_name_id])
which(wcp2$taxon_status!="Accepted" & wcp2$plant_name_id %in% wcp2$accepted_plant_name_id)

table(wcp2$accepted_plant_name_id %in% wcp2$plant_name_id) # this is not zero because some of the accepted IDs point to taxon IDs that have no accepted IDs, and those have been removed earlier. This is alright, has been like that in the old version 


saveRDS(wcp2, "processed_data/apg_wcp_jul_21_clean.rds")
