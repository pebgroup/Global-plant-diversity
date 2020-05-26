rm(list=ls())
library(ape)
library(phytools)
library(tidyverse)


# SET PARAMETERS ########################################################################################
## use all accepted taxa or just species level?
tax_level <- "elevated_to_species_id" # "accepted_plant_name_id" 

#source("plant_sr/functions.R")

# read trees from Smith & Brown 2018
## ALLTB = GenBank and Open Tree of Life taxa with a backbone provided by Open Tree of Life version 9.1
## ALLMB = GenBank and Open Tree of Life taxa with a backbone provided by Magallón et al. (2015)
### Smith & Brown 2018 say "Primarily, the Magallón et al. (2015) and Open Tree of Life back- bones were 
### similar but with the Open Tree of Life backbone providing more resolution toward the tips, that can 
### be useful when there are no molecular data."

## CHOOSE TREE ############################################################################################
phylo <- read.tree("trees/ALLMB.tre")
#phylo_b <- read.tree("trees/ALLOTB.tre")

fin <- data.frame(tips = phylo$tip.label)
fin$tips_mod <- gsub("_", " ", fin$tips) # replace underscore with space


# exclude rogue Asteraceae
#phylo <- drop.tip(phylo, c("Schmidtia_capensis", "Hypochaeris_arachnoides"))
#phylo_b <- drop.tip(phylo_b, c("Schmidtia_capensis", "Hypochaeris_arachnoides"))


# read OTL taxonomy to identify tip labels
# ott <- read.csv("data/ott/taxonomy.tsv",  sep="\t")
# ott <- ott[,!grepl("X", names(ott))]
# saveRDS(ott, "data/ott.rds")
ott <- readRDS("data/ott.rds")




### STATS AND SOURCES FOR TIP LABELS #################################################


#table(fin$tips_mod %in% ncbi$name) / length(fin$tips_mod) # 0.398
table(fin$tips_mod %in% ott$name) / length(fin$tips_mod) # 0.997
      
ott <- ott[ott$name %in% fin$tips_mod,]
ott_sub <- ott[-which(!ott$uniqname==""),] # remove bacteria and marine stuff
rm(ott)

# check sources for the tip labels
table(grepl("gbif", ott_sub$sourceinfo))/nrow(ott_sub)
table(grepl("ncbi", ott_sub$sourceinfo))/nrow(ott_sub)
table(grepl("ncbi|gbif", ott_sub$sourceinfo))/nrow(ott_sub)
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
# read NCBI-WCSP matches from the taxonomy matcher
matches <- readRDS("data/fin_species_match_NCBI.rds")
# combine elevated to species colum with the rest - what happens if species is not defined?
matches$elevated_to_species_id[which(is.na(matches$elevated_to_species_id))] <- matches$accepted_plant_name_id[which(is.na(matches$elevated_to_species_id))]
matches <- matches[,c(tax_level, "id")]

ott_ncbi <- merge(ott_ncbi, matches, 
             by.x="ncbi_id", by.y="id", all.x=TRUE)
table(is.na(ott_ncbi$elevated_to_species_id)) / nrow(ott_ncbi) # 85% matches for ncbi tip labels

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


# export remaining tip labels for common format + taxonomy matcher
#saveRDS(fin[!is.na(fin$gbif_id),], file="SB_tip_labels.rds")

#source("common_format_creator_SB.R")
# data/input_tip_labels.rds
#source(taxonomic_matcher blabla.R)





## Load taxonomy matcher created GBIF tip labels #################################################### 
# data/fin_GBIF.rds.rds
matches <- readRDS("data/fin_species_match_GBIF.rds")
matches$elevated_to_species_id[which(is.na(matches$elevated_to_species_id))] <- matches$accepted_plant_name_id[which(is.na(matches$elevated_to_species_id))]
#matches <- matches[,c(tax_level, "id")]

# merge the GBIF accepted name IDS into fin
fin <- merge(fin, matches[,c("tip", tax_level)],
              by.x="tips_mod", by.y="tip", all.x=TRUE)
rm(matches)






# COMBINE NCBI AND GBIF MATCHES ######################################################################
# combine both accepted plant id columns
col_merge <- grep(tax_level, names(fin))
fin[is.na(fin[,col_merge[1]]),col_merge[1]] <- fin[is.na(fin[,col_merge[1]]),col_merge[2]]
fin <- fin[,-col_merge[2]]
names(fin)[grep(tax_level, names(fin))] <- "accepted_id"

table(is.na(fin$accepted_id))/nrow(fin) # 85% matches total


table(is.na(fin$accepted_id[fin$source=="ncbi"]))/nrow(fin[fin$source=="ncbi",]) # 85%
table(is.na(fin$accepted_id[fin$source=="gbif"]))/nrow(fin[fin$source=="gbif",]) # 84%





##### NON ACCEPTED INVESTIGATION ###################################################################################
# found 250 non accepted or non species taxa on the final tree. investigate where they come from: which matching procedure and so on
# check status in WCSP
fin <- fin[-which(is.na(fin$accepted_id)),]
wcsp <- readRDS("data/WCSP.apg.rds")
wcsp <- wcsp %>% 
  select(-family) %>%
  rename(family = family.apg)

table(fin$accepted_id %in% wcsp$plant_name_id)
phylo.dat<- wcsp[wcsp$plant_name_id %in% fin$accepted_id,]
nrow(phylo.dat)
table(phylo.dat$taxon_status, phylo.dat$taxon_rank)
table(phylo.dat$taxon_rank)
pd.invalid <- phylo.dat[phylo.dat$taxon_status!="Accepted" | phylo.dat$taxon_rank!="Species",]
nrow(pd.invalid)
table(fin$source[fin$accepted_id %in% pd.invalid$plant_name_id])

gbif.matches <- readRDS("data/fin_species_match_GBIF.rds")
ncbi.matches <- readRDS("data/fin_species_match_NCBI.rds")

View(gbif.matches[which(gbif.matches$accepted_plant_name_id %in% pd.invalid$plant_name_id),])
gbif.sub <- gbif.matches[which(gbif.matches$accepted_plant_name_id %in% pd.invalid$plant_name_id),]
table(gbif.sub$taxon_rank)
# all species from gbif
wcsp[wcsp$plant_name_id=="1010729-az",]
wcsp$taxon_status[wcsp$plant_name_id %in% gbif.sub$accepted_plant_name_id]

# ok the gbif matching is wrong. there is synonyms in the accepted plant name IDs

ncbi.sub <- ncbi.matches[which(ncbi.matches$accepted_plant_name_id %in% pd.invalid$plant_name_id),]
table(ncbi.sub$taxon_rank)
# varieties and ssps are from NCBI matching side - elevated to species part mistake?

# how is the status of all outcomes from gbif matching?
table(wcsp$taxon_status[wcsp$plant_name_id %in% gbif.matches$accepted_plant_name_id])
# ok that`s not good....
# mistakes happen in strict match no author, strict match, and matching authors

######################################################################################################################





 ### REPLACE TIP LABELS and RESOlVE MULTIPLE MATCHES #################################################################

# put wcsp IDs on tree
phylo_org <- phylo
phylo$tip.label <- fin$accepted_id[match(phylo$tip.label, fin$tips)]
  
## duplicates
length(which(table(phylo$tip.label)>1)) # 20089 duplicated tip labels (differs between ALLMB and ALLOTB)
  
# load genbank only tree:
phylogb_a <- read.tree("trees/GBMB.tre") # Magallon
phylogb_b <- read.tree("trees/GBOTB.tre") # Opentree
wcsp <- readRDS("data/WCSP.apg.rds") # WCSP

resolve_multiple <- function(MATCHES, wcsp, phylo, phylo_org, phylogb){
    # matches = accepted WCSP ID tip labels, phylo=phylogeny with labels replaced, phylo_org=original tip labels
    # remove multiple linkages, i.e. when multiple tips in the tree are assigned to the same accepted name, using following criteria
    # 1) preferably keep tips that have molecular data behind them
    # 2) preferably keep tips that have the same genus name as the species they link to in WCSP - important for tips that have been added based on taxonomy
    # 3) randomly thereafter
    for(n in names(table(MATCHES)[table(MATCHES)>1])){
      counter <- 1
      tips <- which(MATCHES == n)
      erase <- rep(FALSE,length(tips))
      #first criterion: if any of the tips comes from GenBank, erase all that don't
      if(sum(phylo$tip.label[tips] %in% phylogb$tip.label)>0){
        erase <- !phylo$tip.label[tips] %in% phylogb$tip.label
        #break
      } else {#second criterion (only relevant if choosing among non-genbank tips): if at least one of the original tips has the same genus as in wcsp, erase the rest
        #retrieve WCSP genus
        gen <- as.vector(wcsp[wcsp$plant_name_id==n,"genus"])
        #retrieve tip genera
        GEN <- vector("character", length(tips))
        for(tip in 1:length(tips)){
            
          GEN[tip] <- strsplit(phylo_org$tip.label[tips],"_")[[tip]][1]
        }
        #if at least one of the tips has the same genus as in wcsp, erase the rest
        if(gen %in% GEN){
          erase <- !GEN == gen
        }
      }
      if(sum(!erase)>1){#if there are still n>1 items to be kept (i.e. not erased), chose n-1 randomly for erasing
        erase[sample(which(erase==FALSE), size=(sum(!erase)-1))] <- TRUE
      } 
      
      #erase tips 
      MATCHES[tips[erase]] <- NA
    }
    counter <- counter + 1
    if(!counter%%100) print(counter)# cat(counter,"\r")
    return(MATCHES)
  }
# function as from Wolf, adjusted to work on the phylogeny with replaced tip labels

# 20 MINUTES ###
res_multi <- resolve_multiple(phylo$tip.label, wcsp, phylo, phylo_org, phylogb_a) 
table(is.na(res_multi)) # 84 479 NAs introduced
  
  
  
# DROPPING THE UNUSED TIP LABELS ######################################################################
  
res_multi_noNA <- na.omit(res_multi)
clean_tree <- keep.tip(phylo, as.character(res_multi_noNA))
write.tree(clean_tree, "trees/allmb_matched_no_multi.tre")

# SAVE DATA ############################################################################################

#save(clean_tree, file = "matched_tree.RData")



# # reformat matches for export
# MATCHES_a_a <- data.frame(tip = phylo$tip.label, conservative = MATCHES_a_a)
# for(i in 1:ncol(MATCHES_a_a)) MATCHES_a_a[,i] <- as.vector(MATCHES_a_a[,i])
# MATCHES_b_a <- data.frame(tip = phylo_b$tip.label, conservative = MATCHES_b_a)
# for(i in 1:ncol(MATCHES_b_a)) MATCHES_b_a[,i] <- as.vector(MATCHES_b_a[,i])
# 
# save(MATCHES_a_a, MATCHES_b_a, phylo_a, phylo_b, wcsp, file = "MATCHES.RData")

#save.image("match_data.RData")


  