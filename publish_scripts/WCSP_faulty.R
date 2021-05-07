# Script to check inconsistencies in the WCSP, using download from June 2020 #####

library(tidyverse)
library(data.table)



#wc_all <- readRDS("./data/WCSP.apg.rds")
wc_all <- fread("data/wcs_jun_20/world_checklist_names_and_distribution_jun_2020/checklist_names.txt", quote = "")
dim(wc_all)

# remove all with no accepted ID
wc_all <- wc_all[-which(wc_all$accepted_plant_name_id==""),]

table(wc_all$taxon_status[wc_all$plant_name_id %in% wc_all$accepted_plant_name_id])
View(wc_all[,c(1,2,3,4,5,24)])
# --> Some accepted IDs do not link to accepted species according to taxon status

# Extracting taxa with accepted plant name ID but non accepted status
wc_weird <- wc_all[wc_all$plant_name_id %in% wc_all$accepted_plant_name_id & wc_all$taxon_status!="Accepted",]
table(wc_weird$taxon_status)
# 4348 entries, most of those are hybrids, but there is others


#  Checking what their accepted plant name ID is
# Down the rabbit hole
wc_weird1 <- wc_all[wc_all$plant_name_id %in% wc_weird$accepted_plant_name_id,]
table(wc_weird1$taxon_status)
View(wc_weird1[,c(1,2,3,4,5,24)])
## going one level deeper solves 441 cases (mostly synonyms)

wc_weird2 <- wc_all[wc_all$plant_name_id %in% wc_weird1$accepted_plant_name_id,]
table(wc_weird2$taxon_status)
## going two levels deeper solves a synonym loop (?), but does not lead to more accepted


wc_weird3 <- wc_all[wc_all$plant_name_id %in% wc_weird2$accepted_plant_name_id,]
table(wc_weird3$taxon_status)
# going past second level does not fix more 


# Check the remaining non-accepted
## remove accepted 
wc_weird3 <- wc_weird3[-which(wc_weird3$taxon_status=="Accepted"),]
View(wc_weird3[wc_weird3$taxon_status=="Unplaced",c(1,2,3,4,5,24)])

#529673-wcs points to itself while being marked as "Unplaced"
table(as.character(wc_weird3$plant_name_id) == as.character(wc_weird3$accepted_plant_name_id))
# Almost all of these entries point to themselves


wc_weird4 <- wc_weird3[as.character(wc_weird3$plant_name_id) != as.character(wc_weird3$accepted_plant_name_id),]

# the remaining create a loop:
wc_weird4[,c("plant_name_id", "accepted_plant_name_id")]

View(wc_weird4[,c(1,2,3,4,5,24)])




# Fix the WCSP ####################################################################################

#wcsp <- readRDS("data/WCSP.apg.rds") # old version
wcsp <- readRDS("data/apg_wcp_jun_20.rds")
wcsp <- wcsp[-which(wcsp$accepted_plant_name_id==""),]

# Remove invalid species pointing to themselves
## This is mostly artificial hybrids - do we want to remove them or keep them?
wcsp <- wcsp[-which(as.character(wcsp$plant_name_id)==as.character(wcsp$accepted_plant_name_id) 
                    & wcsp$taxon_status!="Accepted"),]

# Create index to identify faulty entries
ind <- which(wcsp$taxon_status!="Accepted" & wcsp$plant_name_id %in% wcsp$accepted_plant_name_id)
wcsp.tmp <- wcsp[ind,]
table(wcsp.tmp$taxon_status) # status of "supposed to be accepted"-taxa

# level 1 - get accepted IDs from falsely linked as accepted taxa
lev1 <- wcsp[wcsp$plant_name_id %in% wcsp.tmp$accepted_plant_name_id,]
table(lev1$taxon_status)

# subset wcsp.tmp to those species that are accepted in lev 1 -> in those cases, the wcsp.tmp accepted ID is the right accepted ID
wcsp.tmp.sub <- wcsp.tmp[wcsp.tmp$accepted_plant_name_id %in% lev1$plant_name_id[lev1$taxon_status=="Accepted"],]

# merging first level: match original accepted with level 1 plant name id to get new accepted ID. The valid accepted is always the last one added
wcsp2 <- merge(wcsp, wcsp.tmp[,c("plant_name_id", "accepted_plant_name_id")], 
              by.x="accepted_plant_name_id", by.y="plant_name_id",
              all.x=TRUE)

# remove already solved from next matching step
lev1.tmp <- lev1[-which(lev1$taxon_status=="Accepted"),]
lev1.tmp <- lev1.tmp %>%
  rename(accepted_plant_name_id_lev2 = accepted_plant_name_id)

# level 2 - get accepted IDs from still falsely linked as accepted taxa
lev2 <- wcsp[wcsp$plant_name_id %in% lev1.tmp$accepted_plant_name_id_lev2,]
table(lev2$taxon_status)
# 2 more accepted species

# subset lev1.tmp to those species that are accepted in lev 2 -> in those cases, the lev1.tmp accepted ID is the right accepted ID
lev1.tmp.sub <- lev1.tmp[lev1.tmp$accepted_plant_name_id_lev2 %in% lev2$plant_name_id[lev2$taxon_status=="Accepted"],]

# merging second level: the valid accepted is always the last one added
wcsp2 <- merge(wcsp2, lev1.tmp.sub[,c("plant_name_id", "accepted_plant_name_id_lev2")], 
              by.x="accepted_plant_name_id.y", by.y="plant_name_id",
              all.x=TRUE)

View(wcsp2[,c("accepted_plant_name_id", "accepted_plant_name_id.y", "accepted_plant_name_id_lev2")])


# Remove taxa stuck in a loop of synonyms
View(wcsp2[loops,])

loops <- which(wcsp2$plant_name_id %in% lev2$plant_name_id[lev2$taxon_status!="Accepted"])
wcsp2 <- wcsp2[-loops,]



# Cleanup ####
## list to update already existing trees and files etc.
update_list <- wcsp2[!is.na(wcsp2$accepted_plant_name_id.y),]

update_list$accepted_plant_name_id.y[!is.na(update_list$accepted_plant_name_id_lev2)] <- 
  update_list$accepted_plant_name_id_lev2[!is.na(update_list$accepted_plant_name_id_lev2)]
update_list <- update_list %>%
  rename(old_accepted_ID = accepted_plant_name_id) %>%
  rename(new_accepted_ID = accepted_plant_name_id.y) %>%  
  select(c(plant_name_id, old_accepted_ID, new_accepted_ID))


## building a new clean WCSP version
wcsp2$accepted_plant_name_id[!is.na(wcsp2$accepted_plant_name_id.y)] <- wcsp2$accepted_plant_name_id.y[!is.na(wcsp2$accepted_plant_name_id.y)]
wcsp2$accepted_plant_name_id[!is.na(wcsp2$accepted_plant_name_id_lev2)] <- wcsp2$accepted_plant_name_id_lev2[!is.na(wcsp2$accepted_plant_name_id_lev2)]
wcsp2 <- wcsp2 %>%
  select(-c(accepted_plant_name_id.y, accepted_plant_name_id_lev2))

table(wcsp2$taxon_status[wcsp2$plant_name_id %in% wcsp2$accepted_plant_name_id])


saveRDS(wcsp2, "data/apg_wcp_jun_20_clean.rds")
saveRDS(update_list, "data/wcsp_update_list.rds")




