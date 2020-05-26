# compare checklist versions

wcsp.old <- readRDS("data/wcsp_old.rds")
wcsp <- readRDS("data/WCSP.apg.rds")
wcsp.no <- readRDS("data/wcp_dec_19.rds")

# Compate checklists

## create new ID format
id_split <- strsplit(as.character(wcsp.old$checklist_id), "-")
source <- sapply(id_split, "[[", 1)
source[source=="atoz"] <- "az"
new_ID <- paste(sapply(id_split, "[[", 2), source, sep="-")
wcsp.old$checklist_id_new <- new_ID


# Compare new checklist versions
table(wcsp$plant_name_id %in% wcsp.no$plant_name_id) # adjusted families new to original new: all included

table(wcsp.old$checklist_id_new %in% wcsp$plant_name_id) # 5558 IDs not included anymore
table(wcsp$plant_name_id %in% wcsp.old$checklist_id_new) # 17398 new IDs in new version

# Compare missing

old_missing <- wcsp.old[which(!wcsp.old$checklist_id_new %in% wcsp$plant_name_id),]
table(old_missing$taxon_status_description)

View(old_missing[old_missing$taxon_status_description=="Accepted",])





########### compare distribution data interlinked with wcsp #############333
comm <- readRDS("data/comm.rds")
comm_old <- readRDS("data/comm_old_wcsp.rds")
comm_org <- readRDS("comm_old_version.rds")

comb <- data.frame(region =row.names(comm), 
                   sr_new=rowSums(comm), 
                   sr_old = rowSums(comm_old),
                   sr_org = rowSums(comm_org))
plot(log(comb$sr_new), log(comb$sr_org))




     
     