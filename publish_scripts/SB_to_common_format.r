############################
# Taxon matching taxonomy to WCSP #
############################
rm(list=ls())
library(tidyverse)

# chose dataset: BIEN or NCBI 
DB.name <- "SB"



######## comments Melanie ######################################################
# line 8: Moved dataset choice to the top of the script
# lines 39+40 are mandatory to get genus and species into the tax_comb-->taxon_rank column - needs to stay
# multimatch handling: clarfiy matching procedures & adapt the matching workflow chart
# renamed all done[1-9] to "done" to update in stead of producing more variables
# changed order: put no infra rank matching before no family matching
# I removed all entries from the WCSP that do not point to an accepted taxon ID: avoids matching resulting in NA later on!
# added taxon rank for genus + species to the NCBI common format, was missing

### TESTED WITH NCBI + BIEN ###

##### READ WCSP DATA ###########################################################

#wc_all <- readRDS("./data/wcp_dec_19.rds")
wc_all <- readRDS("data/wcp_dec_19.rds")




#### RENAMING WCSP #############################################################

# add new column "tax_comb" to unify taxon ranks and infraspecific ranks column
wc_all$tax_comb <- wc_all$infraspecific_rank # transfer content
wc_all$tax_comb <- gsub("\\.|,", "", wc_all$tax_comb) # clean up, remove dots and commas
empty_ranks <- which(wc_all$tax_comb=="")
wc_all$tax_comb[empty_ranks] <- as.character(wc_all$taxon_rank[empty_ranks])
# MEL: the last 2 lines need to stay; they add the species and genus into the tax_comb column

# remove capitals
wc_all$tax_comb <- gsub("Species", "species", wc_all$tax_comb)
wc_all$tax_comb <- gsub("Genus", "genus", wc_all$tax_comb)

# replace empty cells with NA
wc_all[wc_all==""]<-NA 

# fix hybrid notation
wc_all$genus_hybrid <- as.character(wc_all$genus_hybrid)
wc_all$species_hybrid <- as.character(wc_all$species_hybrid)
wc_all$genus_hybrid[!is.na(wc_all$genus_hybrid)] <- "x"
wc_all$species_hybrid[!is.na(wc_all$species_hybrid)] <- "x"

# subset the dataset with key columns and remove all NA accepted IDs
wc_all_sub <- wc_all %>% 
  select(tax_comb, taxon_status, family, genus, genus_hybrid, species, 
         species_hybrid, infraspecies, taxon_authors, accepted_plant_name_id) %>% 
  unique() %>% 
  filter(!is.na(accepted_plant_name_id))

# rename columns
names(wc_all_sub)[names(wc_all_sub) %in% c("tax_comb", "infraspecies", "taxon_authors")] <-
  c("taxon_rank", "infra_name", "author")

wc_all_sub[,1:ncol(wc_all_sub)] <- 
  lapply(wc_all_sub[,1:ncol(wc_all_sub)], as.character)
# saveRDS(wc_all_sub, file="./results/wc_all_sub.rds")
#wc_all_sub <- readRDS("./results/wc_all_sub.rds")




#### define the dataset to reconcile with WCSP ##############################

#### NCBI ####

if(DB.name=="NCBI"){
  
  dataset <- read.csv(paste0("./data/",DB.name, ".csv"), header=TRUE, stringsAsFactors = FALSE)
  
  #cleaning
  dataset$infraspecific_rank <- gsub("\\.|,", "", dataset$infraspecific_rank)
  dataset[dataset==""]<-NA
  dataset$genus_hybrid[!is.na(dataset$genus_hybrid)] <- "x"
  dataset$species_hybrid[!is.na(dataset$species_hybrid)] <- "x"
  
  #reorder, maybe not neccesary
  dataset <- dataset %>% 
    select(infraspecific_rank, order, family, genus, genus_hybrid, species, 
           species_hybrid, infraspecies, taxon_authority, taxon_rank, ncbi_id) %>%
    unique()
  names(dataset)[names(dataset) %in% c("infraspecific_rank", "infraspecies", "taxon_authority", "taxon_rank", "ncbi_id")] <- 
    c("taxon_rank", "infra_name", "author", "Trank.ncbi.full", "id")
  dataset$author <- gsub("_", " ", dataset$author)
  
  dataset[,1:ncol(dataset)] <- lapply(dataset[,1:ncol(dataset)], as.character)
  
  dataset$taxon_rank[which(is.na(dataset$taxon_rank))] <- dataset$Trank.ncbi.full[which(is.na(dataset$taxon_rank))]
}

#### BIEN ####

if(DB.name=="BIEN"){
  
  bien_input <- readRDS("data/bien_input_vectorized3.rds")
  
  # RENAMING BIEN
  
  # straighten formats
  bien_input$taxon_rank[which(bien_input$taxon_rank=="[unranked]")] <- NA
  bien_input$taxon_rank <- gsub(".", "", bien_input$taxon_rank, fixed=TRUE)
  
  # rename
  bien_input$taxon_rank <- gsub("cv", "convar", bien_input$taxon_rank, fixed=TRUE)
  bien_input$taxon_rank <- gsub("fo", "f", bien_input$taxon_rank, fixed=TRUE)
  bien_input$taxon_rank <- gsub("ssp", "subsp", bien_input$taxon_rank, fixed=TRUE)
  bien_input$taxon_rank <- gsub("unranked", NA, bien_input$taxon_rank, fixed=TRUE)
  bien_input$taxon_rank[which(is.na(bien_input$taxon_rank))] <- "undefined"
  
  bien_input[,1:3] <- lapply(bien_input[,1:3], as.character)
  
  # exclude Bryophyta
  bry <- read.csv("data/bryophyta.csv")
  bryos <- as.character(bry[,1])
  bryos <- gsub(" ", "", bryos)
  bien_input <- bien_input[!bien_input$family %in% bryos,]
  
  # exclude taxa without extracted genus
  bien_input <- bien_input[-which(is.na(bien_input$genus)),]
  
  dataset <- bien_input
  rm(bien_input)
}  




#### TAXONOMY MATCHING - Shared Functionality ####################################

##### STRICT MATCHING, no author #################################################

res <- left_join(dataset, wc_all_sub, all.x=TRUE,
                 by=c("taxon_rank", "family", "genus", "genus_hybrid",
                      "species", "species_hybrid", "infra_name"))

res.tmp <- res %>% 
  select(-author.x, -author.y) %>% 
  unique() %>% 
  mutate(match_type=ifelse(is.na(accepted_plant_name_id), NA, "strict match"))

no_match <- res.tmp %>% 
  filter(is.na(accepted_plant_name_id))
match <- res.tmp %>% 
  filter(!is.na(accepted_plant_name_id))

# get IDs that have multiple matches
mm_strict_id <- unique(match$id[duplicated(match$id)])

done <- match %>% 
  filter(!(id %in% mm_strict_id)) %>% 
  select(id, accepted_plant_name_id, match_type)
mm <- res %>% 
  filter(id %in% mm_strict_id) %>% 
  arrange(id)

# remove taxa from no_match if they occur in done or mm
no_match <- no_match %>% 
  filter(!(no_match$id %in% done$id))
no_match <- no_match %>% 
  filter(!(no_match$id %in% mm$id))

# last overall check
try(if(all(unique(c(mm$id, no_match$id, done$id)) %in% dataset$id)){
  print("all IDs represented!")
}else{
  stop("NOT all IDs represented!")
})




###### MULTIMATCHES ##############################################################

# simplify author names: remove all spaces, dots and numbers
mm$author.x <- gsub(" |\\.|[0-9]*|,", "", mm$author.x) # author names from dataset
mm$author.y <- gsub(" |\\.|[0-9]*|,", "", mm$author.y) # author names from WCSP

# assign matching authors status 
mm$match_type <- ifelse(mm$author.x == mm$author.y, "matching authors", "no matching authors")

### I realized ""matching authors"" as not as efficent to resue all the names, unless comnined with "taxon_status=="Accepted"
### since you can tell the name is corretly matched, even though mm$author.x != mm$author.y (e.g., id==1099790, 1288263)
### or run code below
# Talked to Wolf about this, and it was decided that we will go for the conservative matching
# and count cases that cannot be matched by author as unresolvable. We also added another 
# check: if multimatches point to the same checklist ID, they also count as solved.

#### define a function to check multiple matches ####
multi_match_checker <- function(mm_ids, subdata){
  valid_matches <- c("matching authors", "strict match no family", "no infra (same same)")
  one_match <- NULL
  more_match <- NULL
  zero_match <- NULL
  resolved_mm <- NULL
  unresolved_mm <- NULL
  
  for (id in unique(mm_ids)){
    temp <- subdata[subdata$id==id,]
    # all entries point to the same accepted ID
    if(length(unique(temp$accepted_plant_name_id)) == 1){
      one_match <- c(one_match, id)
      resolved_mm <- rbind(resolved_mm, temp)
      # entries point to different accepted IDs
    }else{
      # one entry with matching authors
      if(sum(temp$match_type %in% valid_matches, na.rm=TRUE)==1){    
        one_match <- c(one_match, id)
        resolved_mm <- rbind(resolved_mm, temp[which(temp$match_type %in% valid_matches),])
      }
      # > 1 entries with matching authors
      if(sum(temp$match_type %in% valid_matches, na.rm=TRUE) > 1){
        more_match <- c(more_match, id)
        unresolved_mm <- rbind(unresolved_mm, temp)
      }
      # no matches because missing authors
      if(all(is.na(temp$match_type))){
        zero_match <- c(zero_match, id)
        unresolved_mm <- rbind(unresolved_mm, temp)
      }
    }
  }
  results <- list(
    "one_match" = one_match,
    "more_match" = more_match,
    "zero_match" = zero_match,
    "resolved_mm" = resolved_mm,
    "unresolved_mm" = unresolved_mm
  )
  return(results)
}

mm_results <- multi_match_checker(mm_strict_id, mm)

wcp_conflicts1 <- mm_results$unresolved_mm %>% 
  filter(id %in% mm_results$more_match)


#output some summary reports
# 
# check results
# length(unique(mm$id))
# length(mm_results$one_match) + length(mm_results$more_match) + length(mm_results$zero_match)
# there is species missing: those with no matching authors that are also not NA
# View(mm[!mm$id %in% c(mm_results$one_match, mm_results$more_match, mm_results$zero_match),])
#
# if(!dir.exists("results")){
#   dir.create("results")
#   write.csv(wcp_conflicts1, "./results/wcp_conflicts1.csv", row.names = FALSE, quote=FALSE)
#   
#   pdf(paste0("./results/wcp_unresolved_conflicts_", DB.name, ".pdf"))
#   par(mfrow=c(2,1), mar = c(7, 4, 2, 1))
#   barplot(rev(tail(sort(table(mm_results$unresolved_mm$family)), 30)),
#           las=2, ylab="Taxa", cex.axis=0.65, cex.names = 0.8)
#   title(main = "Famies with the most unresolved names")
#   
#   barplot(rev(tail(sort(table(wcp_conflicts$family)), 40)),
#           las=2, ylab="Taxa", cex.names = 0.8)
#   title(main = "Families with the most conflict names")
#   
#   dev.off()
#   
# }


# update done data set with resolved entries
done <- rbind(done, mm_results$resolved_mm %>% 
                select(id, accepted_plant_name_id, match_type))





###### NO STRICT MATCHES #########################################################

## no infraspecific match: if infra == species name → match on species level

no_match <- res %>%
  filter(res$id %in% no_match$id) %>%
  select(-author.y, -taxon_status, -accepted_plant_name_id) %>%
  unique() %>%
  rename(author = author.x)

## exclude taxa that have been matched to avoid double assignment
nms <- no_match %>% 
  filter(infra_name == species & !(id %in% done$id))

res3 <- left_join(nms, wc_all_sub, all.x=TRUE, 
                  by=c("family", "genus", "genus_hybrid",
                       "species", "species_hybrid", "author")) %>% 
  mutate(match_type=ifelse(is.na(accepted_plant_name_id), NA, "no infra (same same)"))

## get species ID with multiple matches
mm_strict_id <- unique(res3$id[duplicated(res3$id)])
res3_sub <- res3 %>% 
  filter(!(id %in% mm_strict_id) & !is.na(match_type))

# attach single matches to the resolved dataframe
done <- rbind(done, res3_sub %>% 
                select(id, accepted_plant_name_id, match_type))

# generate multimatches dataset
res3_mm <- res3 %>% 
  filter(id %in% mm_strict_id)

mm_results3 <- multi_match_checker(mm_strict_id, res3_mm)

# attach the resolved multimaches to done
done <- rbind(done, mm_results3$resolved_mm %>% 
                select(id, accepted_plant_name_id, match_type))




## no family, but authors

## exclude taxa that have been matched already to avoid double assignment
no_match <- no_match %>% 
  filter(!(id %in% done$id))

res2 <- left_join(no_match, wc_all_sub, all.x=TRUE, 
                  by=c("taxon_rank", "genus", "genus_hybrid",
                       "species", "species_hybrid", "infra_name", "author")) %>% 
  mutate(match_type=ifelse(is.na(accepted_plant_name_id), NA, "strict match no family"))

## get species ID with multiple matches
mm_strict_id <- unique(res2$id[duplicated(res2$id)])

## attach single matches to the resolved dataframe
res2_sub <- res2 %>% 
  filter(!(res2$id %in% mm_strict_id) & !is.na(match_type))
done <- rbind(done, res2_sub %>% 
                select(id, accepted_plant_name_id, match_type))

# generate multimatches dataset, altough differtent families
res2_mm <- res2 %>% 
  filter(id %in% mm_strict_id)

if(length(mm_strict_id)>0){
  mm_results2 <- multi_match_checker(mm_strict_id, res2_mm)
  
  # attach the resolved multimaches to done
  done <- rbind(done, mm_results2$resolved_mm %>% 
                  select(id, accepted_plant_name_id, match_type))
}









####### final check and accepted ID assignment ##############################

# no double assignments?
table(tapply(done$accepted_plant_name_id, done$id, function(x)length(unique(x))), useNA="ifany")

# remove NAs + no matching authors from converging multimatch merge
dupli_match_types <- done$id[which(duplicated(done$id))]
done <- done %>% 
  filter(!is.na(match_type) | !is.na(id)) %>%
  unique()
done <- done[-which(done$id %in% dupli_match_types & done$match_type=="no matching authors"),]

# check for double assignments
if(any(duplicated(done$id))){
  dupl.id <- done$id[duplicated(done$id)]
  print("double assignment happend - check double_assignments for cases")
  double_assignments <- done %>% filter(id %in% dupl.id)
}else print("No double assignments. All is well!")

#write.csv(done, paste0("./results/wcsp_",DB.name," _matching_list_done.csv", sep=""), row.names = FALSE, quote=FALSE)

# attach done dataframe to input dataset
fin <- left_join(dataset, done, by="id", all.x=TRUE)

# combine unresolved into one file
unsolved <- bind_rows(mm_results$unresolved_mm, mm_results2$unresolved_mm, mm_results3$unresolved_mm)
saveRDS(unsolved, file=paste0("./results/unsolved_", DB.name , ".rds"))

# assign multimatch 
if(any(fin$id %in% unsolved$id)){
  fin$match_type[fin$id %in% unsolved$id] <- "multimatch"
}

# save output
saveRDS(fin, file=paste0("./results/fin_", DB.name , ".rds"))

# table(fin$match_type, useNA = "ifany")
# table(fin$match_type, useNA = "ifany")/nrow(fin)
# table(fin$family[fin$match_type %in% c("multimatch")], useNA = "ifany")
# table(fin$taxon_rank[is.na(fin$match_type)])









#### SPECIES LEVEL MATCHING #########################################################
# matching logic: get accepted ID for all subspecies entries 
# --> get genus + species name from those accepted IDs names
# --> create another dataframe from wcsp only containing species level entries
# --> match the subspecies-level dataframe with the species-level dataframe based on genus and species name
# --> attach the species-level accepted ID to the final dataframe

## subset wcsp to species occurring in BIEN 
ids <- unique(fin$accepted_plant_name_id)
wc_sub <- wc_all[wc_all$plant_name_id %in% ids,]

## subset to taxon rank < species
wc_subspecies <- wc_sub[!wc_sub$tax_comb %in% c("species", "genus"),]

# get species-level entries from wcsp. (do they have to be accepted?)
wc_species <- wc_all[wc_all$tax_comb %in% c("species", "genus") & wc_all$taxon_status=="Accepted",]

# match subsepcies with species based on which have the same genus + species entry
## plant name id from subspecies will be used to match to the fin dataframe using accepted ID
species_match <- left_join(wc_subspecies[,c("plant_name_id", "species", "genus")],
                           wc_species[,c("accepted_plant_name_id", "species", "genus")],
                           by=c("species", "genus"),
                           all.x=TRUE)

nrow(species_match)-nrow(wc_subspecies)
# there is some multimatches

species_match$elevated_to_species_id <- species_match$accepted_plant_name_id

## get species ID with multiple matches
multimatch_id <- unique(species_match$plant_name_id[which(duplicated(species_match$plant_name_id))])
species_mms <- species_match[species_match$plant_name_id %in% multimatch_id,]

species_match_sub <- species_match[!species_match$plant_name_id %in% multimatch_id,]
species_match_sub <- species_match_sub[!is.na(species_match_sub$accepted_plant_name_id),]

#add to fin
fin_species_match <- merge(fin, species_match_sub[,c("plant_name_id", "elevated_to_species_id")],
                           by.x="accepted_plant_name_id",
                           by.y="plant_name_id", all.x=TRUE)

# save output
saveRDS(fin_species_match, file=paste0("./results/fin_species_match_", DB.name, ".rds"))

# View(fin_species_match[,c("scrubbed_taxon_name_no_author", "match_type", "accepted_plant_name_id", "elevated_to_species_id")])

# Note: excample "Equisetum variegatum var. alaskanum" has no elevated to species level ID, because the accepted plant name ID is at species level already
# If the original taxon name input was species level, and it still gets an elevated to species level ID, that means the species name points to an accepted subspecies, 
# which again was matched with its parent sepcies name  (see for example "Centaurea asperula")
