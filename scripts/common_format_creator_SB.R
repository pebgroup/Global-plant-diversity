# builds a common format based on the Smith & Brown ALLMB tree tip labels from GBIF to feed into the taxonomy matcher
rm(list=ls())
library(tidyverse)


# Load Smith & Brown tip labels not covered by NCBI 
sb <- readRDS("SB_tip_labels.rds")



# GBIF ##################################################################################
# download runs via server:
# gbif_server.rds

gbif <- readRDS("data/gbif_all.rds")

# transform list to dataframe
## removing lines that do not include all required columns
cols <- c("key", "scientificName", "authorship",
          "taxonomicStatus", "rank", "family", "genus",
          "species")
func1 <- function(x){all(cols %in% names(x))}
remo <- which(unlist(lapply(gbif, func1))==FALSE)
non.null.res <- gbif[-remo]

taxonID <- sapply(non.null.res, "[[", "key")
scientificName <- sapply(non.null.res, "[[", "scientificName")
authorship <- sapply(non.null.res, "[[", "authorship")
taxonomicStatus <- sapply(non.null.res, "[[", "taxonomicStatus")
rank <- sapply(non.null.res, "[[", "rank")
#kingdom <- sapply(res, "[[", "kingdom")
family <- sapply(non.null.res, "[[", "family")
genus <- sapply(non.null.res, "[[", "genus")
species <- sapply(non.null.res, "[[", "species")

gbif <- data.frame(taxonID, 
                   scientificName,
                   authorship,
                   taxonomicStatus,
                   rank,
                   family,
                   genus,
                   species)

all(gbif$taxonID %in% sb$gbif_id)

sb <- merge(sb, gbif, by.x="gbif_id", by.y="taxonID", all.x=TRUE)
table(is.na(sb$family))/nrow(sb)




# set up common format ###################################################################
split_length <- unlist(lapply(strsplit(as.character(sb$tips_mod), split = " "), length))
input <- data.frame(tip=sb$tips_mod,
                    taxon_name = sb$scientificName,
                    family = sb$family,
                    author= sb$authorship, 
                    split_length=split_length,
                    genus_hybrid = NA,
                    species_hybrid = NA,
                    species = sub(".* ", "", sb$species),
                    genus = sb$genus,
                    taxon_rank = tolower(sb$rank),
                    infra_name = NA,
                    comment = NA,
                    usable = NA)


## Cleaning
# remove sp. species
input <- input[!grepl("sp\\.", input$tip),]

# remove aggregates
input <- input[!grepl("aggr$", input$tip),]

# order the dataframe by split length
input <- input[order(input$split_length),]
input$id <- c(1:nrow(input))

# create split list
split_list <- strsplit(as.character(input$tip), split = " ")
names(split_list) <- input$id








# Vectorize + loop mix ###############################################################
  
## split length == 1
ind <- which(input$split_length==1)
input$usable[input$split_length==1] <- "no"

## split length == 2
ind <- which(input$split_length==2)

# # get authors
# input$author[ind] <- gsub("^[A-Z][a-z| ]+", "", input[ind, "taxon_name"])

# # get the rest
# input$genus[ind] <- sapply(split_list[ind], "[[", 1)
# input$species[ind] <- sapply(split_list[ind], "[[", 2)
# input$taxon_rank[ind] <- "species"

ind.hyb <- which(grepl("^×", input$genus))
input$genus[ind.hyb] <- gsub("^×", "", sapply(split_list[ind.hyb], "[[", 1))
input$genus_hybrid[ind.hyb] <- "x"
input$usable[ind.hyb] <- "no"


## split length == 3
ind <- which(input$split_length==3)
for(i in 1:length(ind)){
  if(split_list[[ind[i]]][2]=="sp."){
#    input$genus[ind[i]] <- split_list[[ind[i]]][1]
#    input$species[ind[i]] <- split_list[[ind[i]]][2]
#    input$taxon_rank[ind[i]] <- "species"
    input$usable[ind[i]] <- "no"
    
  }
  if(split_list[[ind[i]]][2]=="x"){ # does not occur but you never know
    if(grepl("[A-Z]",split_list[[ind[i]]][3])){
      input$genus_hybrid[ind[i]] <- "x"
      input$usable[ind[i]] <- "no"
    }else{
      input$species_hybrid[ind[i]] <- "x"
#      input$genus[ind[i]] <- split_list[[ind[i]]][[1]]
#      input$species[ind[i]] <- split_list[[ind[i]]][[3]]
#      input$taxon_rank[ind[i]] <- "species"
    }
  }
  if(split_list[[ind[i]]][1]!="x" & split_list[[ind[i]]][2]!="x" & !grepl("[A-Z]",split_list[[ind[i]]][3])){ 
#    input$genus[ind[i]] <- split_list[[ind[i]]][[1]]
#    input$species[ind[i]] <- split_list[[ind[i]]][[2]]
#    input$taxon_rank[ind[i]] <- NA
    input$infra_name[ind[i]] <- split_list[[ind[i]]][[3]]
    }
  if(split_list[[ind[i]]][2]!="x" & grepl("[A-Z]",split_list[[ind[i]]][3])){
#    input$genus[ind[i]] <- split_list[[ind[i]]][[1]]
#    input$species[ind[i]] <- split_list[[ind[i]]][[2]]
#    input$taxon_rank[ind[i]] <- "species"
#    input$author[ind[i]] <- split_list[[ind[i]]][[3]]
  }
  if(grepl("sp\\.", split_list[[ind[i]]][[2]])){
    input$usable[ind[i]] <- "no"
  }
}

# # get authors
# input$author[ind] <- gsub("^[A-Z][a-z| |-]+", "", input[ind, "taxon_name"])

# # manually fix la Croix
# a <- which(grepl("^Croix.*", input$author))
# input$author[a] <- paste0("la ", input$author[a])
# input$infra_name[a] <- NA
# input$taxon_rank[a] <- "species"

# manually remove 2 letter third parts that are not Author names - this would be fixing tip labels, do we want that?


# remove all unusable taxa
input <- input[-which(input$usable=="no"),]

  
saveRDS(input, "data/input_tip_labels.rds")

