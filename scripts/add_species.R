rm(list=ls())
source("scripts/functions.R")

library(ape)
library(phytools)
#library(tidyverse)
library(dplyr)

phylo <- read.tree("trees/allmb_matched_no_multi.tre")
wcsp <- readRDS("data/WCSP.apg.rds")
wcsp <- wcsp %>% 
  select(-family) %>%
  rename(family = family.apg)

# remove ferns + moss
ferns <- as.vector(read.csv("data/fern_list.txt")[,1])
wcsp <- wcsp[!wcsp$family %in% ferns,]
wcsp <- wcsp[!wcsp$family %in% "Isoetaceae",] # special character thing


# get list of all accepted species
goodspp <- wcsp[wcsp$genus_hybrid == "" &
                  wcsp$species_hybrid == "" &
                  wcsp$infraspecific_rank == "" &
                  wcsp$species != "" &
                  wcsp$taxon_status == "Accepted",]
goodspp2 <- wcsp[wcsp$taxon_status == "Accepted" & wcsp$taxon_rank!="Genus",]
nrow(goodspp2)

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
fams_in_tree <- unique(as.vector(goodspp[goodspp$plant_name_id %in% as.vector(phylo$tip.label),"family"]))

# families not in the tree yet: fams=tree families, families=to be added families
lefties <- families[which(!families %in% fams_in_tree)]
# the weirdos... monospecific families e.g.
lefties[which(lefties %in% unresolved_families)]
nrow(goodspp[goodspp$family %in% lefties,])
# 31 species cannot be added





# reduce toadd for testing purposes
nrow(toadd)
toadd <- toadd[sample(c(1:nrow(toadd)), 500, replace=FALSE),]

# initialize vector to record which method to use to find MRCA (genus = 1, family = 2, or order = 3)
method <- vector("numeric", nrow(toadd))

# species that can be added to genus
method[toadd$genus %in% genera] <- 1


# species that can be added to family
method[method==0][as.vector(toadd[method==0,"family"]) %in% fams_in_tree] <- 2

# species that can be added to order
method[method==0][!as.vector(toadd[method==0,"family"]) %in% fams_in_tree] <- 3





MATCHES <- phylo$tip.label
taxa_not_added <- c()
nodeheight.phylo <- max(nodeHeights(phylo))
tip_number <- c()
nh <- c()

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
    tip_number <- c(tip_number, length(tips)) # diagnostics
    next()
  }
  
  if(length(tips) == 1){ # if adding to a terminal branch (terminal branches happen when there is a single species in a taxon e.g.)
    tip_number <- c(tip_number, length(tips)) # diagnostics
    nn <- which(phylo$tip.label == tips)
    
    # fix edge.length removal
    # len <- phylo$edge.length[which(phylo$edge[, 2] == nn)] *0.5
    # len <- phylo$edge.length[which(phylo$tip.label==tips)] # add same length as terminal branch  
    # phylo2 <- bind.tip2(phylo, newname, where = nn,
    #                      position = 0.5 * phylo$edge.length[which(phylo$edge[, 2] == nn)],
    #                      edge.length = len)
          # position is setting where the new tip would be attached on the branch. conservatively this should be 0.5 and would cause polytomies if not specified, but not if there is only one tip in the group!
    phylo <- bind.tip2(phylo, newname, where = nn,
                        position = 0.5 * phylo$edge.length[which(phylo$edge[, 2] == nn)],
                        edge.length = 0.5 * phylo$edge.length[which(phylo$edge[, 2] == nn)])
  } else { # if adding to an internal node  
    tip_number <- c(tip_number, length(tips)) # diagnostics
    nn <- findMRCA(phylo, tips = tips)
    #tree_sub <- extract.clade(phylo, nn)
    
    #TEST plot old node ################################################
    # plotTree(tree_sub, ftype="off")
    # tip.names <- tips # store for repeated usage
    # text<-tip.names
    # tips<-sapply(text,function(x,y) which(y==x),y=tree_sub$tip.label)
    # linklabels(text,tips,link.type="curved")
    # nodelabels()
    ####################################################################3
    
    # fix edge.length removal
#    len <- max(makeL(tree_sub)) # root of subtree to any tip
#    len <- distRoot(tree_sub, tips[1]) # root of subtree to one tip
    len <- nodeheight.phylo - nodeheight(phylo, nn)
    phylo <- bind.tip2(phylo, as.character(newname), where = nn, edge.length = len)
  }
  
  #TEST plot new clade ################################################
#   tree_sub_new <- extract.clade(phylo, nn+1) # nodenumber changes by +1 when tip gets added
#   plotTree(tree_sub_new, ftype="off")
#   text<-c(tip.names, as.character(newname))
# #  a<-c(rep("black",length(tip.names)), "red")
#   tips<-sapply(text,function(x,y) which(y==x),y=tree_sub_new$tip.label)
#   linklabels(text,tips,link.type="curved") #, col=a
#   nodelabels()
  ####################################################################3  
  
  MATCHES <- c(MATCHES, newname)
  
  if(i %% ceiling(nrow(toadd)/100) == 0) print(paste(i/ceiling(nrow(toadd)/100), "% complete", Sys.time()))
  nh <- c(nh, max(nodeHeights(phylo)))
}

#write.tree(phylo, "allmb_added_species.tre")

nh
table(phylo$edge.length<0)
table(method)
phylo$tip.label[which(phylo$edge.length<0)]
table(tip_number==1) # how many terminal branches

nodeheight.phylo
max(nodeHeights(phylo))






#### CHECK SERVER RESULTS ##############################################################################

fin_tree <- read.tree("trees/allmb_matched_added_species.tre")
fin_tree2 <- read.tree("trees/allmb_matched_added_species_2.tre")

er <- readRDS("data/taxa_not_added.rds")
er2 <- readRDS("data/taxa_not_added_2.rds")

setdiff(er, er2)
all(er2 %in% er)

phylo <- read.tree("trees/allmb_matched_no_multi.tre")
# wcsp <- readRDS("data/WCSP.apg.rds")
# wcsp <- wcsp %>% select(-family) %>%
#   rename("family" = "family.apg")

length(fin_tree$tip.label)-length(phylo$tip.label)
length(fin_tree2$tip.label)-length(phylo$tip.label)
any(duplicated(fin_tree$tip.label))
any(duplicated(fin_tree2$tip.label))

#View(wcsp[wcsp$plant_name_id %in% er,]) # Osmundaceae = fern
#View(wcsp[wcsp$plant_name_id %in% er2,]) # Osmundaceae = fern

# Check branch lengths
table(fin_tree$edge.length<0)
hist(log(fin_tree$edge.length))
table(fin_tree2$edge.length<0)

#fin_tree$edge[which(fin_tree$edge.length<0)]

# 34 species not added

plot(fin_tree, "fan", show.tip.label = FALSE)
plot(fin_tree2, "fan", show.tip.label = FALSE)

library(castor)
dist.phylo <- get_all_distances_to_root(phylo, as_edge_count=FALSE)
dist <- get_all_distances_to_root(fin_tree, as_edge_count=FALSE)
dist2 <- get_all_distances_to_root(fin_tree2, as_edge_count=FALSE)
hist(dist[1:Ntip(fin_tree)])
hist(dist.phylo[1:Ntip(phylo)])
hist(dist2[1:Ntip(fin_tree2)])

table(dist.phylo[1:Ntip(phylo)]>325 & dist.phylo[1:Ntip(phylo)]<326)
table(dist[1:Ntip(fin_tree)]>325 & dist[1:Ntip(fin_tree)]<326)
table(dist2[1:Ntip(fin_tree2)]>325 & dist2[1:Ntip(fin_tree2)]<326)
table(dist2>325 & dist2<326)

hist(dist2[which(fin_tree2$tip.label %in% toadd$plant_name_id)])
hist(dist2[1:Ntip(fin_tree2)])

fin_tree2$tip.label[which(dist2[1:Ntip(fin_tree2)] > 50000)]


fin_tree <- untangle(fin_tree)











# Check wich status the tip labels have in the WCSP
phylo.dat <- wcsp[which(as.character(wcsp$plant_name_id) %in% fin_tree2$tip.label),] 
nrow(phylo.dat)
table(phylo.dat$taxon_status)

## before adding species
phylo.dat.org <- wcsp[which(as.character(wcsp$plant_name_id) %in% phylo$tip.label),] 
nrow(phylo.dat.org)
table(phylo.dat.org$taxon_status)

## There is some non- accepted species in the Snith&Brown matched version.
table(phylo.dat.org$taxon_rank, phylo.dat.org$taxon_status)







