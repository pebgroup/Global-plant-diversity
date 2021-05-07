# script to transform initial WCSP to RDS format

library(data.table)
wcp19 <- fread("data/wcs_dec_19/checklist_names.txt")
wcp20 <- fread("data/wcs_jun_20/world_checklist_names_and_distribution_jun_2020/checklist_names.txt")

#saveRDS(wcsp, "data/wcsp/wcp_jun_20.rds")

# quick version check
wcsp <- readRDS("data/apg_wcp_jun_20_clean.rds") # file produced in BIEN project
