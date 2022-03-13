# matching tip labels with WCVP IDs, part I

wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str())) 
library(phytools)

tax_level <- "elevated_to_species_id" 
phylo <- read.tree("../data/ALLMB.tre") # phylogeny from Smith&Brown 2018
ott <- readRDS("../data/ott.rds") # OTL taxonomy to identify tip labels
matches <- readRDS("../data/fin_species_match_NCBI_wcvp21.rds") # NCBI-WCSP matches from the taxonomy matcher

fin <- data.frame(tips = phylo$tip.label)
fin$tips_mod <- gsub("_", " ", fin$tips) # replace underscore with space



### STATS AND SOURCES FOR TIP LABELS #################################################


#table(fin$tips_mod %in% ncbi$name) / length(fin$tips_mod) # 0.398
table(fin$tips_mod %in% ott$name) / length(fin$tips_mod) # 0.997
      
ott <- ott[ott$name %in% fin$tips_mod,]
ott_sub <- ott[-which(!ott$uniqname==""),] # remove bacteria and marine taxa
rm(ott)

# check sources for the tip labels
table(grepl("gbif", ott_sub$sourceinfo))/nrow(ott_sub)
table(grepl("ncbi", ott_sub$sourceinfo))/nrow(ott_sub)
table(grepl("ncbi|gbif", ott_sub$sourceinfo))/nrow(ott_sub) # covers almost every entry
table(grepl("ncbi.*gbif|gbif.*ncbi", ott_sub$sourceinfo))/nrow(ott_sub)


# mark tip labels according to source
ott_ncbi <- ott_sub[grepl("ncbi", ott_sub$sourceinfo),]
fin$source <- NA
fin$source[which(fin$tips_mod %in% ott_ncbi$name)] <- "ncbi"
fin$source[which(is.na(fin$source))] <- "gbif"






## get NCBI id from OTT #######################################################
one <- gsub("ncbi:", "", ott_ncbi$sourceinfo)
ott_ncbi$ncbi_id <- gsub(",gbif:[0-9].*", "", one)
ott_ncbi$ncbi_id <- gsub(",.*", "", ott_ncbi$ncbi_id)

# MATCH NCBI TIP LABELS 
# combine elevated to species column with the rest - what happens if species is not defined?
matches$elevated_to_species_id[which(is.na(matches$elevated_to_species_id))] <- matches$accepted_plant_name_id[which(is.na(matches$elevated_to_species_id))]
matches <- matches[,c(tax_level, "id")]

ott_ncbi <- merge(ott_ncbi, matches, 
             by.x="ncbi_id", by.y="id", all.x=TRUE)
table(is.na(ott_ncbi$elevated_to_species_id)) / nrow(ott_ncbi) # 85.5% matches 

## merge into fin
fin <- merge(fin, ott_ncbi[,c("name", "ncbi_id", tax_level)], 
              by.x="tips_mod", by.y="name", all.x=TRUE)
rm(matches)


## get GBIF id from OTT #######################################################
ott_gbif <- ott_sub[!ott_sub$uid %in% ott_ncbi$uid,]
ott_gbif$gbif_id <- gsub(".*gbif:", "", ott_gbif$sourceinfo)
ott_gbif$gbif_id <- gsub(",.*", "", ott_gbif$gbif_id)
fin <- merge(fin, ott_gbif[,c("name", "gbif_id")], 
             by.x="tips_mod", by.y="name", all.x=TRUE)


# export tip labels that have not been matched with NCBI
saveRDS(fin[!is.na(fin$gbif_id),], file="../processed_data/SB_tip_labels_21.rds")
# export all tip labels with source
saveRDS(fin, file="../processed_data/fin.rds")

