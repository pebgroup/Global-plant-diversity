
### 15 HOURS ON SERVER ####################################################

rm(list=ls())
bind.tip2 <- function (tree, tip.label, edge.length = NULL, where = NULL, 
                       position = 0, interactive = FALSE, ...) 
{
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  use.edge.length <- if (is.null(tree$edge.length)) 
    FALSE
  else TRUE
  if (use.edge.length == FALSE) 
    tree <- compute.brlen(tree)
  if (interactive == TRUE) {
    plotTree(tree, ...)
    cat(paste("Click where you would like to bind the tip \"", 
              tip.label, "\"\n", sep = ""))
    flush.console()
    obj <- get.treepos(message = FALSE)
    where <- obj$where
    position <- obj$pos
  }
  else if (is.null(where)) 
    where <- Ntip(tree) + 1
  if (where <= Ntip(tree) && position == 0) {
    pp <- 1e-12
    if (tree$edge.length[which(tree$edge[, 2] == where)] <= 
        1e-12) {
      tree$edge.length[which(tree$edge[, 2] == where)] <- 2e-12
      ff <- TRUE
    }
    else ff <- FALSE
  }
  else pp <- position
  if (is.null(edge.length) && is.ultrametric(tree)) {
    H <- nodeHeights(tree)
    if (where == (Ntip(tree) + 1)) 
      edge.length <- max(H)
    else edge.length <- max(H) - H[tree$edge[, 2] == where, 
                                   2] + position
  }
  tip <- list(edge = matrix(c(2, 1), 1, 2), tip.label = tip.label, 
              edge.length = edge.length, Nnode = 1)
  class(tip) <- "phylo"
  obj <- bind.tree(tree, tip, where = where, position = pp)
  if (where <= Ntip(tree) && position == 0) {
    nn <- obj$edge[which(obj$edge[, 2] == which(obj$tip.label == 
                                                  tip$tip.label)), 1]
    obj$edge.length[which(obj$edge[, 2] == nn)] <- obj$edge.length[which(obj$edge[, 
                                                                                  2] == nn)] + 1e-12
    obj$edge.length[which(obj$edge[, 2] == which(obj$tip.label == 
                                                   tip$tip.label))] <- 0
    obj$edge.length[which(obj$edge[, 2] == which(obj$tip.label == 
                                                   tree$tip.label[where]))] <- 0
  }
  #root.time <- if (!is.null(obj$root.time)) 
  #  obj$root.time
  #else NULL
  #obj <- untangle(obj, "read.tree")
  #if (!is.null(root.time)) 
  #  obj$root.time <- root.time
  #if (interactive) 
  #  plotTree(obj, ...)
  if (!use.edge.length) 
    obj$edge.length <- NULL
  obj
}
      

library(ape)
library(phytools)
library(tidyverse)

phylo <- read.tree("trees/allmb_matched_no_multi.tre")
wcsp <- readRDS("data/apg_wcp_jun_20_clean.rds")
wcsp <- wcsp %>% select(-family) %>%
  rename("family" = "family.apg")

# remove ferns + moss
ferns <- as.vector(read.csv("fern_list.txt")[,1])
wcsp <- wcsp[!wcsp$family %in% ferns,]
wcsp <- wcsp[!wcsp$family %in% "Isoetaceae",] # special character thing


# get list of all accepted species
goodspp <- wcsp[wcsp$genus_hybrid == "" &
                  wcsp$species_hybrid == "" &
                  wcsp$infraspecific_rank == "" &
                  wcsp$species != "" &
                  wcsp$taxon_status == "Accepted",]

# get list of good species that are not in matching
toadd <- goodspp[!goodspp$plant_name_id %in% phylo$tip.label,]

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
nrow(toadd[toadd$family %in% lefties,])






# reduce toadd for testing purposes
#toadd <- toadd[sample(c(1:nrow(toadd)), 500, replace=FALSE),]

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

saveRDS(taxa_not_added, "taxa_not_added_2.rds")
write.tree(phylo, "allmb_matched_added_species_Nov20.tre")
