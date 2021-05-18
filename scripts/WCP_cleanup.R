# Script to check inconsistencies in the wcp, using download from June 2020 #####

library(tidyverse)
library(data.table)


# Fix the wcp ######


wcp <- wcp[-which(wcp$accepted_plant_name_id==""),]

# Remove invalid species pointing to themselves
## This is mostly artificial hybrids - do we want to remove them or keep them?
wcp <- wcp[-which(as.character(wcp$plant_name_id)==as.character(wcp$accepted_plant_name_id) 
                    & wcp$taxon_status!="Accepted"),]

# Create index to identify faulty entries
ind <- which(wcp$taxon_status!="Accepted" & wcp$plant_name_id %in% wcp$accepted_plant_name_id)
wcp.tmp <- wcp[ind,]
table(wcp.tmp$taxon_status) # status of "supposed to be accepted"-taxa

# level 1 - get accepted IDs from falsely linked as accepted taxa
lev1 <- wcp[wcp$plant_name_id %in% wcp.tmp$accepted_plant_name_id,]
table(lev1$taxon_status)

# subset wcp.tmp to those species that are accepted in lev 1 -> in those cases, the wcp.tmp accepted ID is the right accepted ID
wcp.tmp.sub <- wcp.tmp[wcp.tmp$accepted_plant_name_id %in% lev1$plant_name_id[lev1$taxon_status=="Accepted"],]

# merging first level: match original accepted with level 1 plant name id to get new accepted ID. The valid accepted is always the last one added
wcp2 <- merge(wcp, wcp.tmp[,c("plant_name_id", "accepted_plant_name_id")], 
              by.x="accepted_plant_name_id", by.y="plant_name_id",
              all.x=TRUE)

# remove already solved from next matching step
lev1.tmp <- lev1[-which(lev1$taxon_status=="Accepted"),]
lev1.tmp <- lev1.tmp %>%
  rename(accepted_plant_name_id_lev2 = accepted_plant_name_id)

# level 2 - get accepted IDs from still falsely linked as accepted taxa
lev2 <- wcp[wcp$plant_name_id %in% lev1.tmp$accepted_plant_name_id_lev2,]
table(lev2$taxon_status)
# 2 more accepted species

# subset lev1.tmp to those species that are accepted in lev 2 -> in those cases, the lev1.tmp accepted ID is the right accepted ID
lev1.tmp.sub <- lev1.tmp[lev1.tmp$accepted_plant_name_id_lev2 %in% lev2$plant_name_id[lev2$taxon_status=="Accepted"],]

# merging second level: the valid accepted is always the last one added
wcp2 <- merge(wcp2, lev1.tmp.sub[,c("plant_name_id", "accepted_plant_name_id_lev2")], 
              by.x="accepted_plant_name_id.y", by.y="plant_name_id",
              all.x=TRUE)

#View(wcp2[,c("accepted_plant_name_id", "accepted_plant_name_id.y", "accepted_plant_name_id_lev2")])


# Remove taxa stuck in a loop of synonyms
#View(wcp2[loops,])

loops <- which(wcp2$plant_name_id %in% lev2$plant_name_id[lev2$taxon_status!="Accepted"])
wcp2 <- wcp2[-loops,]



# Cleanup ####
## list to update already existing trees and files etc.
update_list <- wcp2[!is.na(wcp2$accepted_plant_name_id.y),]

update_list$accepted_plant_name_id.y[!is.na(update_list$accepted_plant_name_id_lev2)] <- 
  update_list$accepted_plant_name_id_lev2[!is.na(update_list$accepted_plant_name_id_lev2)]
update_list <- update_list %>%
  rename(old_accepted_ID = accepted_plant_name_id) %>%
  rename(new_accepted_ID = accepted_plant_name_id.y) %>%  
  select(c(plant_name_id, old_accepted_ID, new_accepted_ID))


## building a new clean wcp version
wcp2$accepted_plant_name_id[!is.na(wcp2$accepted_plant_name_id.y)] <- wcp2$accepted_plant_name_id.y[!is.na(wcp2$accepted_plant_name_id.y)]
wcp2$accepted_plant_name_id[!is.na(wcp2$accepted_plant_name_id_lev2)] <- wcp2$accepted_plant_name_id_lev2[!is.na(wcp2$accepted_plant_name_id_lev2)]
wcp2 <- wcp2 %>%
  select(-c(accepted_plant_name_id.y, accepted_plant_name_id_lev2))

table(wcp2$taxon_status[wcp2$plant_name_id %in% wcp2$accepted_plant_name_id])


#saveRDS(update_list, "processed_data/wcp_update_list.rds")




