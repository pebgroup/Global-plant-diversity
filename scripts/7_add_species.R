# Taxonomic adding of missing species, creating polytomies
# this code cannot be paralleled as the tree changes with each iteration :(

wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str())) 
library(phytools)
library(dplyr)
library(tidyr)
source("0_functions.R")

phylo <- read.tree("../processed_data/allmb_matched_clean.tre")
wcp <- readRDS("../processed_data/apg_wcp_jul_21_clean.rds")
ferns <- as.vector(read.csv("../data/fern_list.txt")[,1])


wcp <- wcp %>% 
  dplyr::select(-family) %>%
  dplyr::rename("family" = "family.apg")

# remove ferns
wcp <- wcp[!wcp$family %in% ferns,]
wcp <- wcp[!wcp$family %in% "Isoetaceae",] # remove fern family not caught because of special character


# get list of all accepted species
# goodspp <- wcp[wcp$genus_hybrid == "" &
#                   wcp$species_hybrid == "" &
#                   wcp$infraspecific_rank == "" &
#                   wcp$species != "" &
#                   wcp$taxon_status == "Accepted",]
goodspp <- wcp[is.na(wcp$genus_hybrid) &
              is.na(wcp$species_hybrid) &
              is.na(wcp$infraspecific_rank) &
              !is.na(wcp$species) &
              wcp$taxon_status == "Accepted",]
saveRDS(goodspp, "../processed_data/goodspp.rds")

# get list of good species that are not in matching
toadd <- goodspp[!goodspp$plant_name_id %in% phylo$tip.label,] #nrow=63,677
toadd <- as.data.frame(toadd)
goodspp <- as.data.frame(goodspp)

# create index for orders
{unresolved_families <- list()
  unresolved_families[["Cynomoriaceae"]] <- c("Altingiaceae", "Aphanopetalaceae", "Cercidiphyllaceae", "Crassulaceae", "Daphniphyllaceae", "Grossulariaceae", "Haloragaceae", "Hamamelidaceae", "Iteaceae", "Paeoniaceae", "Penthoraceae", "Peridiscaceae", "Saxifragaceae", "Tetracarpaeaceae")
  unresolved_families[["Circaeasteraceae"]] <- c("Berberidaceae", "Eupteleaceae", "Lardizabalaceae", "Menispermaceae", "Papaveraceae", "Ranunculaceae")
  unresolved_families[["Mitrastemonaceae"]] <- c("Actinidiaceae", "Balsaminaceae", "Cyrillaceae", "Clethraceae", "Diapensiaceae", "Ebenaceae", "Ericaceae", "Fouquieriaceae", "Lecythidaceae", "Marcgraviaceae", "Pentaphylacaceae", "Polemoniaceae", "Primulaceae", "Roridulaceae", "Sapotaceae", "Sarraceniaceae", "Sladeniaceae", "Styracaceae", "Symplocaceae", "Tetrameristaceae", "Theaceae")
  unresolved_families[["Cytinaceae"]] <- c("Bixaceae", "Cistaceae", "Dipterocarpaceae", "Malvaceae", "Muntingiaceae", "Neuradaceae", "Sarcolaenaceae", "Sphaerosepalaceae", "Thymelaeaceae")
  unresolved_families[["Physenaceae"]] <- c("Achatocarpaceae", "Aizoaceae", "Amaranthaceae", "Anacampserotaceae", "Ancistrocladaceae", "Asteropeiaceae", "Barbeuiaceae", "Basellaceae", "Cactaceae", "Caryophyllaceae", "Didiereaceae", "Dioncophyllaceae", "Droseraceae", "Drosophyllaceae", "Frankeniaceae", "Gisekiaceae", "Halophytaceae", "Kewaceae", "Limeaceae", "Lophiocarpaceae", "Macarthuriaceae", "Microteaceae", "Molluginaceae", "Montiaceae", "Nepenthaceae", "Nyctaginaceae", "Phytolaccaceae", "Plumbaginaceae", "Polygonaceae", "Portulacaceae", "Rhabdodendraceae", "Petiveriaceae", "Sarcobataceae", "Simmondsiaceae", "Stegnospermataceae", "Talinaceae", "Tamaricaceae")
  unresolved_families[["Tetracarpaeaceae"]] <- c("Altingiaceae", "Aphanopetalaceae", "Cercidiphyllaceae", "Crassulaceae", "Cynomoriaceae", "Daphniphyllaceae", "Grossulariaceae", "Haloragaceae", "Hamamelidaceae", "Iteaceae", "Paeoniaceae", "Penthoraceae", "Peridiscaceae", "Saxifragaceae")
}

# get the families of the species that remain to be added
families <- as.vector(unique(toadd$family))

# get a list of genera and families that are already represented in the tree
genera <- unique(as.vector(goodspp[goodspp$plant_name_id %in% as.vector(phylo$tip.label),"genus"]))
fams <- unique(as.vector(goodspp[goodspp$plant_name_id %in% as.vector(phylo$tip.label),"family"]))

# families not in the tree yet: fams=tree families, families=to be added families
lefties <- families[which(!families %in% fams)]
lefties[which(lefties %in% unresolved_families)]
nrow(toadd[toadd$family %in% lefties,]) # 9 entries with families that are not in the tree






# reduce toadd for testing purposes
#toadd <- toadd[sample(c(1:nrow(toadd)), 200, replace=FALSE),]

# initialize vector to record which method to use to find MRCA (genus = 1, family = 2, or order = 3)
method <- vector("numeric", nrow(toadd))

# species that can be added to genus
method[toadd$genus %in% genera] <- 1


# species that can be added to family
method[method==0][as.vector(toadd[method==0,"family"]) %in% fams] <- 2

# species that can be added to order
method[method==0][!as.vector(toadd[method==0,"family"]) %in% fams] <- 3





MATCHES <- phylo$tip.label
taxa_not_added <- c()
nodeheight.phylo <- max(nodeHeights(phylo))

for(i in 1:nrow(toadd)){
  
  # define new taxon name
  newname <- toadd[i,"plant_name_id"]
    #gsub(pattern = " ", replacement = "_", toadd[i,"plant_name_id"])
  # if name already exists in tree, append "_WLE" to make it unique
  if(newname %in% phylo$tip.label) newname <- paste(newname, "_WLE", sep="")
  
  # gather the tips representing the higher taxon of interest
  if(method[i] == 1){ # from genus
    tips = phylo$tip.label[as.vector(phylo$tip.label) %in% as.vector(goodspp[goodspp$genus == as.vector(toadd[i,"genus"]),"plant_name_id"])]
  }
  if(method[i] == 2){ # from family
    tips = phylo$tip.label[as.vector(phylo$tip.label) %in% as.vector(goodspp[goodspp$family == as.vector(toadd[i,"family"]),"plant_name_id"])]
  }
  if(method[i] == 3){ # from order
    tips = phylo$tip.label[as.vector(phylo$tip.label) %in% as.vector(goodspp[goodspp$family %in% unresolved_families[[as.vector(toadd[i,"family"])]],"plant_name_id"])]
  }
  
  if(length(tips) < 1){
    print("Attempting to add to a group that is not represented in the tree!")
    taxa_not_added <- c(taxa_not_added, as.character(newname))
    next
  }
  
  if(length(tips) == 1){ # if adding to a terminal branch
    nn <- which(phylo$tip.label == tips)
    #len <- phylo$edge.length[phylo$tip.label==tips] # add same length as terminal branch  
    phylo <- bind.tip2(phylo, as.character(newname), where = nn,
                         position = 0.5 * phylo$edge.length[which(phylo$edge[, 2] == nn)],
                         edge.length = 0.5 * phylo$edge.length[which(phylo$edge[, 2] == nn)])
  } else { # if adding to an internal node  
    nn <- findMRCA(phylo, tips = tips)
    len <- nodeheight.phylo - nodeheight(phylo, nn)
    phylo <- bind.tip2(phylo, as.character(newname), where = nn, edge.length = len)
  }
  
  MATCHES <- c(MATCHES, newname)
  
  if(i %% ceiling(nrow(toadd)/100) == 0) print(paste(i/ceiling(nrow(toadd)/100), "% complete", Sys.time()))
}


# write.tree(phylo, "../processed_data/allmb_matched_added_species_Sep21.tre")




# test results from adding species ############################################
#
# check ultrametric?
is.ultrametric(tree_mod)
is.ultrametric(phylo)
library(castor)
all_distances = get_all_distances_to_root(tree_mod)
tip_distances = all_distances[1:length(tree_mod$tip.label)]
range(tip_distances)
all_distances.sb = get_all_distances_to_root(phylo)
tip_distances.sb = all_distances.sb[1:length(phylo$tip.label)]
range(tip_distances.sb)
# same as in the original tree, not changed by adding taxa

# make ultrametric (use tip extension to avoid NxN matrix storage issues)
## from the vignett: This function provides a quick-and-dirty way to make a tree ultrametric, or to correct small numericalinaccuracies in supposed-to-be ultrametric trees.
tree_mod <- extend_tree_to_height(tree_mod, new_height=NULL)
ultra_tree <- tree_mod[[1]]
is.ultrametric(ultra_tree)
write.tree(ultra_tree, "../processed_data/allmb_matched_added_species_Sep21_clean.tre")

