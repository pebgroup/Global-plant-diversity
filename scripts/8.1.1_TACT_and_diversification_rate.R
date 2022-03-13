
# Diversification Rate ----------------------------------------------------


wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str())) 
library(phytools)
library(data.table)
library(phytools)
library(castor)


# PREP ##################################################################
# *** GENBANK TREE -------------------------------------------------------

phylo <- read.tree("../data/GBMB.tre") # genbank tree to identify tips with molecular data
matches <- readRDS("../data/fin_species_match_NCBI_wcvp21.rds") # load WCVP NCBI matching
ott <- readRDS("../data/ott.rds") # OTL taxonomy to identify tip labels
wcvp <- readRDS("../processed_data/apg_wcp_jul_21_clean.rds")


# is binary?
is.binary(phylo)
table(tabulate(phylo$edge[,1]))
nr <- tabulate(phylo$edge[,1])[which(tabulate(phylo$edge[,1])>2)]
(poly <- which(tabulate(phylo$edge[,1])>2))
edge.nr <- which(phylo$edge[,1] %in% poly) # identify edges
edge.col <- rep("black", length=nrow(phylo$edge))
edge.col[edge.nr] <- "red" # color those poly edges
edge.wid <- rep(1, length=nrow(phylo$edge))
edge.wid[edge.nr] <- 1
plot(phylo, typ="fan", show.tip.label = FALSE, edge.color = edge.col, 
     edge.width = edge.wid)
labelnodes("3", node=poly[1], interactive=F, shape="rect")
labelnodes("3", node=poly[2], interactive=F, shape="rect")
labelnodes("47", node=poly[3], interactive=F, shape="rect")
labelnodes("7", node=poly[4], interactive=F, shape="rect")
labelnodes("3", node=poly[5], interactive=F, shape="rect")
labelnodes("3", node=poly[6], interactive=F, shape="rect")
labelnodes("5", node=poly[7], interactive=F, shape="rect")
labelnodes("3", node=poly[8], interactive=F, shape="rect")
labelnodes("4", node=poly[9], interactive=F, shape="rect")

forty7 <- getDescendants(phylo, poly[3])
phylo$tip.label[forty7]
five <- getDescendants(phylo, poly[7])
phylo$tip.label[five]

set.seed(1234)
phylo <- multi2di(phylo)


fin <- data.frame(tips = phylo$tip.label)
fin$tips_mod <- gsub("_", " ", fin$tips) # replace underscore with space

ott <- ott[ott$name %in% fin$tips_mod,]
ott_sub <- ott[-which(!ott$uniqname==""),] # remove non-plant taxa
rm(ott)

# mark tip labels according to source
ott_ncbi <- ott_sub[grepl("ncbi", ott_sub$sourceinfo),]
fin$source <- NA
fin$source[which(fin$tips_mod %in% ott_ncbi$name)] <- "ncbi"
fin$source[which(is.na(fin$source))] <- "gbif"






## *** get NCBI id from OTT ###################################################
one <- gsub("ncbi:", "", ott_ncbi$sourceinfo)
ott_ncbi$ncbi_id <- gsub(",gbif:[0-9].*", "", one)
ott_ncbi$ncbi_id <- gsub(",.*", "", ott_ncbi$ncbi_id)

# MATCH NCBI TIP LABELS 
# combine elevated to species column with the rest - what happens if species is not defined?
matches$elevated_to_species_id[which(is.na(matches$elevated_to_species_id))] <- matches$accepted_plant_name_id[which(is.na(matches$elevated_to_species_id))]

# do all duplicated IDs have the same genus?? (crit 2) --> NO
lunique <- function(x){length(unique(x))}
test <- tapply(matches$genus, matches$elevated_to_species_id, lunique)
table(test)

matches <- matches[,c("elevated_to_species_id", "id")]
ott_ncbi <- merge(ott_ncbi, matches, 
                  by.x="ncbi_id", by.y="id", all.x=TRUE)
table(is.na(ott_ncbi$elevated_to_species_id)) / nrow(ott_ncbi) # 92% matches

## merge into fin
fin <- merge(fin, ott_ncbi[,c("name", "ncbi_id", "elevated_to_species_id")], 
             by.x="tips_mod", by.y="name", all.x=TRUE)
names(fin)[grep("elevated_to_species_id", names(fin))] <- "accepted_id"


# actual replacement
phylo_org <- phylo
phylo$tip.label <- fin$accepted_id[match(phylo$tip.label, fin$tips)]

## duplicates (we have duplicates because we are working on species level,
## therefore subspecies etc have been assigned all the same ID)
length(which(table(phylo$tip.label)>1)) # duplicated tip labels: 3,518



# in contrast to the ALLMB tree, we know that all tips have genetic info behind
# them, so we modify our multi resolver function for this
resolve_multiple_gbmb <- function(MATCHES, wcp, phylo, phylo_org){
  # matches = accepted WCSP ID tip labels, phylo=phylogeny with labels replaced,
  # phylo_org=original tip labels remove multiple linkages, i.e. when multiple
  # tips in the tree are assigned to the same accepted name, using following
  # criteria 1) preferably keep tips that have molecular data behind them 2)
  # preferably keep tips that have the same genus name as the species they link
  # to in WCSP - important for tips that have been added based on taxonomy 3)
  # randomly thereafter
  for(n in names(table(MATCHES)[table(MATCHES)>1])){
    counter <- 1
    tips <- which(MATCHES == n)
    erase <- rep(FALSE,length(tips))
    # #first criterion: if any of the tips comes from GenBank, erase all that don't
    # if(sum(phylo$tip.label[tips] %in% phylogb$tip.label)>0){
    #   erase <- !phylo$tip.label[tips] %in% phylogb$tip.label
    #   #break
    #} else {
    #second criterion: if at least one of the original tips has the same genus
    #as in wcvp, erase others
    #retrieve WCVP genus
    gen <- as.vector(wcp[wcp$plant_name_id==n,"genus"])
    #retrieve tip genera
    GEN <- vector("character", length(tips))
    for(tip in 1:length(tips)){
      
      GEN[tip] <- strsplit(phylo_org$tip.label[tips],"_")[[tip]][1]
    }
    #if at least one of the tips has the same genus as in wcsp, erase the rest
    if(gen %in% GEN){
      erase <- !GEN == gen
    }
    #    }
    if(sum(!erase)>1){#if there are still n>1 items to be kept (i.e. not erased), chose n-1 randomly for erasing
      erase[sample(which(erase==FALSE), size=(sum(!erase)-1))] <- TRUE
    } 
    
    #erase tips 
    MATCHES[tips[erase]] <- NA
  }
  counter <- counter + 1
  if(!counter%%1) print(counter)# cat(counter,"\r")
  return(MATCHES)
}

# 2 MINUTES ###
Sys.time()
(res_multi <- resolve_multiple_gbmb(phylo$tip.label, wcvp, phylo, phylo_org)) 
table(is.na(res_multi)) # NAs introduced, 7285 multis
Sys.time()

# drop unused tips
res_multi_noNA <- na.omit(res_multi)
clean_tree <- keep.tip(phylo, as.character(res_multi_noNA))
write.tree(clean_tree, "../processed_data/gbmb_matched.tre")


# *** Fix some non monophyletic clades ----------------------------------------

tree <- read.tree("../processed_data/gbmb_matched.tre")
goodsp <- fread("../processed_data/goodsp.csv")
goodsp.sub <- goodsp[goodsp$accepted_plant_name_id%in%tree$tip.label,]
any(tree$edge.length==0)
tree$edge.length[tree$edge.length==0] <- 1.6e-05
out <- castor::extend_tree_to_height(tree)
out$max_extension
t <- out$tree
is.ultrametric(t, tol=1e-20) # T
is.binary(t)

# check monophylies
t0 <-t
goodsp <- fread("../processed_data/goodsp.csv")
goodsp.t0 <- goodsp[goodsp$accepted_plant_name_id %in% t0$tip.label,]

orders <- sort(unique(goodsp.t0$order))
families <- sort(unique(goodsp.t0$family))
genera <- sort(unique(goodsp.t0$genus))
res <- data.frame(unique(goodsp.t0[,c("order", "family", "genus")], ord.mono=NA, fam.mono=NA, gen.mono=NA))
res$ord.mono <- NA
res$fam.mono <- NA
res$gen.mono <- NA
for(i in 1:length(orders)){
  temp <- goodsp.t0[goodsp.t0$order %in% orders[i],]
  res[which(res$order==orders[i]), 4] <- is_monophyletic(t0, temp$accepted_plant_name_id)
  if(!i%%1)cat(i,"\r")
}
for(i in 1:length(families)){
  temp <- goodsp.t0[goodsp.t0$family %in% families[i],]
  res[which(res$family==families[i]), 5] <- is_monophyletic(t0, temp$accepted_plant_name_id)
  if(!i%%1)cat(i,"\r")
}
for(i in 1:length(genera)){
  temp <- goodsp.t0[goodsp.t0$genus %in% genera[i],]
  res[which(res$genus==genera[i]), 6] <- is_monophyletic(t0, temp$accepted_plant_name_id)
  if(!i%%1)cat(i,"\r")
}

nm <- res
nm.orders <- unique(nm$order[!nm$ord.mono])

res <- data.frame(nm.orders, all.in=NA, others.in=NA)
too.much <- list()
missing <- list()
for(i in 1:length(nm.orders)){
  wcvp.taxa <- goodsp.sub$accepted_plant_name_id[goodsp.sub$order%in%nm.orders[i]]
  #tkp <- keep.tip(t, goodsp.sub$accepted_plant_name_id[goodsp.sub$order%in%nm.orders[i]])
  #ord <- to.be.tacted$accepted_plant_name_id[to.be.tacted$order %in% nm.orders[i]]
  # get the tree order clade
  node.name <- t$node.label[grep(nm.orders[i], t$node.label)]
  t.ord.clade <- extract.clade(t, node.name)
  
  # all tips in the actual order node?
  res[i, "all.in"] <- all(wcvp.taxa%in%t.ord.clade$tip.label)
  if(!all(wcvp.taxa%in%t.ord.clade$tip.label)){
    missing[[nm.orders[i]]] <- wcvp.taxa[!wcvp.taxa%in%t.ord.clade$tip.label]
  }
  
  # extra species in the order node?
  res[i, "others.in"] <- any(!t.ord.clade$tip.label%in%wcvp.taxa)
  if(any(!t.ord.clade$tip.label%in%wcvp.taxa)){
    too.much[[nm.orders[i]]] <- t.ord.clade$tip.label[!t.ord.clade$tip.label%in%wcvp.taxa]
  }
  
  if(!i%%1)cat(i,"\r")
}
res
all(unlist(too.much) %in% unlist(missing))
rogue <- unique(c(unlist(too.much),unlist(missing)))

# fix, check + save
tp <- drop.tip(t, rogue)
goodsp.tp <- goodsp[goodsp$accepted_plant_name_id%in%tp$tip.label,]

for(i in 1:length(unique(goodsp.tp$order))){
  temp <- goodsp.tp[goodsp.tp$order %in% unique(goodsp.tp$order)[i],]
  print(is_monophyletic(tp, temp$accepted_plant_name_id))
}
# all TRUE

write.tree(tp, "../processed_data/gbmb_matched_no_misplaced.tre")




# *** Write Goodsp as clean csv -----------------------------------------------


goodsp <- readRDS("../processed_data/goodspp.rds")
# remove those where accepted != plant id
goodsp <- goodsp[-which(goodsp$accepted_plant_name_id != goodsp$plant_name_id),]
names(goodsp)
mod <- goodsp[,c("order", "family", "genus", "accepted_plant_name_id")]
any(duplicated(mod$accepted_plant_name_id))

# all tips in the taxonomy file?
table(clean_tree$tip.label %in% mod$accepted_plant_name_id)
wcvp[wcvp$plant_name_id %in% clean_tree$tip.label[which(!clean_tree$tip.label %in% mod$accepted_plant_name_id)],]
# marked as hybrids, remove them from the tree
clean_tree <- drop.tip(clean_tree, clean_tree$tip.label[which(!clean_tree$tip.label %in% mod$accepted_plant_name_id)])



fwrite(mod, sep=",", file="../processed_data/goodsp.csv")


