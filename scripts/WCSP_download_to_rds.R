# script to transform initial WCSP to RDS format

library(data.table)
wcsp <- fread("data/wcs_dec_19/checklist_names.txt")

saveRDS(wcsp, "data/wcsp/wcp_jun_20.rds")
