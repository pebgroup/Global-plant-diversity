############################
# Taxon matching taxonomy to WCVP #
############################
library(tidyr)
library(data.table)


##### READ + ADJUST WCVP DATA #################################################

wc_all <- readRDS(paste0(data_folder_path, wcvp_input_filename))

# add new column "tax_comb" to unify taxon ranks and infraspecific ranks column
wc_all$tax_comb <- wc_all$infraspecific_rank # transfer content
wc_all$tax_comb <- gsub("\\.|,", "", wc_all$tax_comb) # clean up, remove dots and commas
empty_ranks <- which(is.na(wc_all$tax_comb)) # get empty infraspec ranks
wc_all$tax_comb[empty_ranks] <- as.character(wc_all$taxon_rank[empty_ranks])

# remove capitals
wc_all$tax_comb <- gsub("Species", "species", wc_all$tax_comb)
wc_all$tax_comb <- gsub("Genus", "genus", wc_all$tax_comb)

# # replace empty cells with NA
# wc_all[wc_all==""]<-NA 
# 
# fix hybrid notation
#wc_all$genus_hybrid <- as.character(wc_all$genus_hybrid)
#wc_all$species_hybrid <- as.character(wc_all$species_hybrid)
wc_all$genus_hybrid[!is.na(wc_all$genus_hybrid)] <- "x"
wc_all$species_hybrid[!is.na(wc_all$species_hybrid)] <- "x"

# subset the dataset with key columns and remove all NA accepted IDs
wc_all_sub <- wc_all %>% 
  dplyr::select(tax_comb, taxon_status, family.apg, genus, genus_hybrid, species, 
         species_hybrid, infraspecies, taxon_authors, accepted_plant_name_id) %>% 
  unique() %>% 
  filter(!is.na(accepted_plant_name_id))

# rename columns
names(wc_all_sub)[names(wc_all_sub) %in% c("tax_comb", "infraspecies", "taxon_authors")] <-
  c("taxon_rank", "infra_name", "author")


#### NCBI ####

if(DB.name=="NCBI"){
  
  # dataset <- read.csv(paste0("./data/",DB.name, ".csv"), header=TRUE, stringsAsFactors = FALSE)
  dataset <- readRDS(paste0(data_folder_path, ncbi_input_filename))
  #cleaning
  dataset$infraspecific_rank <- gsub("\\.|,", "", dataset$infraspecific_rank)
  dataset[dataset==""]<-NA
  dataset$genus_hybrid[!is.na(dataset$genus_hybrid)] <- "x"
  dataset$species_hybrid[!is.na(dataset$species_hybrid)] <- "x"
  
  #reorder, maybe not necessary
  dataset <- dataset %>% 
    dplyr::select(infraspecific_rank, order, family.apg, genus, genus_hybrid, species, 
           species_hybrid, infraspecies, taxon_authority, taxon_rank, ncbi_id) %>%
    unique()
  names(dataset)[names(dataset) %in% c("infraspecific_rank", "infraspecies", "taxon_authority", "taxon_rank", "ncbi_id")] <- 
    c("taxon_rank", "infra_name", "author", "Trank.ncbi.full", "id")
  dataset$author <- gsub("_", " ", dataset$author)
  
  dataset[,1:ncol(dataset)] <- lapply(dataset[,1:ncol(dataset)], as.character)
  
  dataset$taxon_rank[which(is.na(dataset$taxon_rank))] <- dataset$Trank.ncbi.full[which(is.na(dataset$taxon_rank))]
  }


#### GBIF ####

if(DB.name=="GBIF"){
  
  input <- readRDS(paste0(data_folder_path, gbif_input_filename))
  
  
  # RENAMING
  # remove old family field
  input <- input %>% 
    dplyr::select(-family)
  
  # rename
  input$taxon_rank <- gsub("variety", "var", input$taxon_rank, fixed=TRUE)
  input$taxon_rank[which(is.na(input$taxon_rank))] <- "undefined"
  
  # set data classes right: all factors need to be changed to characters
  classes <- as.vector(unlist(lapply(input, class)))
  input[,which(classes=="factor")] <- lapply(input[,which(classes=="factor")], as.character)
  
  # exclude taxa without extracted genus
  if(any(is.na(input$genus))){
    input <- input[-which(is.na(input$genus)),]      
  }
  
  dataset <- input
  rm(input)
}  




#### TAXONOMY MATCHING - Shared Functionality #################################

##### STRICT MATCHING, no author ##############################################

res <- left_join(dataset, wc_all_sub, all.x=TRUE,
              by=c("taxon_rank", "family.apg", "genus", "genus_hybrid",
                   "species", "species_hybrid", "infra_name"))

res.tmp <- res %>% 
  dplyr::select(-author.x, -author.y) %>% 
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
  dplyr::select(id, accepted_plant_name_id, match_type)
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



###### MULTIMATCHES #########################################################

# simplify author names: remove all spaces, dots and numbers
mm$author.x <- gsub(" |\\.|[0-9]*|,", "", mm$author.x) # author names from dataset
mm$author.y <- gsub(" |\\.|[0-9]*|,", "", mm$author.y) # author names from WCVP

# assign matching authors status 
mm$match_type <- ifelse(mm$author.x == mm$author.y, "matching authors", "no matching authors")

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

Sys.time()
mm_results <- multi_match_checker(mm_strict_id, mm) # should not take longer than 2 minutes
Sys.time()

wcp_conflicts1 <- mm_results$unresolved_mm %>% 
  dplyr::filter(id %in% mm_results$more_match)

# update done data set with resolved entries
done <- rbind(done, mm_results$resolved_mm %>% 
                dplyr::select(id, accepted_plant_name_id, match_type))





###### NO STRICT MATCHES #####################################################

## no infraspecific match: if infra == species name → match on species level

no_match <- res %>%
  dplyr::filter(res$id %in% no_match$id) %>%
  dplyr::select(-author.y, -taxon_status, -accepted_plant_name_id) %>%
  unique() %>%
  dplyr::rename(author = author.x)

## exclude taxa that have been matched to avoid double assignment
nms <- no_match %>% 
  filter(infra_name == species & !(id %in% done$id))

res3 <- left_join(nms, wc_all_sub, all.x=TRUE, 
              by=c("family.apg", "genus", "genus_hybrid",
                   "species", "species_hybrid", "author")) %>% 
  mutate(match_type=ifelse(is.na(accepted_plant_name_id), NA, "no infra (same same)"))

## get species ID with multiple matches
mm_strict_id <- unique(res3$id[duplicated(res3$id)])
res3_sub <- res3 %>% 
  filter(!(id %in% mm_strict_id) & !is.na(match_type))

# attach single matches to the resolved dataframe
done <- rbind(done, res3_sub %>% 
                dplyr::select(id, accepted_plant_name_id, match_type))

# generate multimatches dataset
res3_mm <- res3 %>% 
  filter(id %in% mm_strict_id)

mm_results3 <- multi_match_checker(mm_strict_id, res3_mm)

# attach the resolved multimaches to done
if(!all(lapply(mm_results3, length)==0)){
  done <- rbind(done, mm_results3$resolved_mm %>% 
                  dplyr::select(id, accepted_plant_name_id, match_type))
}




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
                dplyr::select(id, accepted_plant_name_id, match_type))

# generate multimatches dataset, altough differtent families
res2_mm <- res2 %>% 
  filter(id %in% mm_strict_id)

if(length(mm_strict_id)>0){
  mm_results2 <- multi_match_checker(mm_strict_id, res2_mm)
  
  # attach the resolved multimaches to done
  done <- rbind(done, mm_results2$resolved_mm %>% 
                  dplyr::select(id, accepted_plant_name_id, match_type))
}else{
  mm_results2 <- NULL
}


####### final check and accepted ID assignment ##############################

# no double assignments? should all be 1
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


# attach done data frame to input data set
fin <- left_join(dataset, done, by="id", all.x=TRUE)

# combine unresolved into one file
unsolved <- bind_rows(mm_results$unresolved_mm, mm_results2$unresolved_mm, mm_results3$unresolved_mm)

# assign multimatch 
if(any(fin$id %in% unsolved$id)){
  fin$match_type[fin$id %in% unsolved$id] <- "multimatch"
}

wcsp.acc.name <- wc_all %>% 
  dplyr::select(taxon_name, taxon_status, accepted_plant_name_id) %>% 
  filter(taxon_status=="Accepted")

wcsp.w.acc.name <- left_join(fin, wcsp.acc.name, by="accepted_plant_name_id")



#### SPECIES LEVEL MATCHING ##################################################
# matching logic: get accepted ID for all subspecies entries 
# --> get genus + species name from those accepted IDs names
# --> create another dataframe from wcsp only containing species level entries
# --> match the subspecies-level dataframe with the species-level dataframe based on genus and species name
# --> attach the species-level accepted ID to the final dataframe

## subset wcsp to species present in input data
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
# there are XXX multimatches

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
saveRDS(fin_species_match, file=paste0(results_folder_path, fin_file_species))


# Note that some species will not have species level ID assigned, that is
# because the accepted plant name ID is at species level already (example:
# "Equisetum variegatum var. alaskanum") If the original taxon name input was
# species level, and it still gets an elevated to species level ID assigned,
# that means the species name points to an accepted subspecies, which again was
# matched with its parent species name  (see for example "Centaurea asperula")
