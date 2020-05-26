# arguments for manual parallelization on server
#args = commandArgs(trailingOnly=TRUE)
#x = as.numeric(args[1])

rm(list=ls())
library(ape)
library(parallel)
source("scripts/functions.R")

# retrieve link tables, phylogenies, and wcsp data (from add_species.R)
phylo <- read.tree("trees/allmb_matched_added_species_2.tre") 
wcsp <- readRDS("data/WCSP.apg.rds")  

# retrieve community matrix (from process_geography.R)
load("data/comm.RData") # object name = comm

# stats ######
# wcsp_acc <- wcsp[wcsp$taxon_status=="Accepted",]
# table(wcsp_acc$plant_name_id %in% phylo$tip.label)/nrow(wcsp_acc)
# table(wcsp_acc$plant_name_id %in% phylo$tip.label)
# 
# table(colnames(comm) %in% phylo$tip.label)

####################################################
# 1. Computing tree-level variables (RD and EDGES) #
####################################################

 keep <- sample(1:length(phylo$tip.label), 10000)
 drop <- 1:length(phylo$tip.label)
 drop <- drop[-keep]
 testtree <- drop.tip(phylo,phylo$tip.label[drop])
 rm(drop, keep)

trees <- c("testtree")

# this takes >1h

for(i in 1:length(trees)){
  print(Sys.time())
  assign(paste0("RD.", trees[i]), unlist(root.distance(get(trees[i]), mc.cores = 4)))
  #saveRDS(get(paste0("RD.", trees[i])), paste0("RD.", trees[i], ".rds"))
  print(paste(paste0("RD.", trees[i]), "calculated.", Sys.time()))

  # assign(paste("EDGES.", trees[i], sep=""), mclapply(1:Ntip(get(trees[i])), get_edges, phylo=get(trees[i]), mc.cores=4))
  # #saveRDS(get(paste("EDGES.", trees[i], sep="")), paste("EDGES.", trees[i], ".rds", sep=""))
  # print(paste(paste("EDGES.", trees[i], sep=""), "calculated.", Sys.time()))
}
rm(i)

# 10 000 tip label: 20 seconds for RD, 
    
####################
# 2. Calculate MRD #
####################

# analyses <- cbind(
#   c("phylo_a_a_a", "phylo_b_a_a"),
#   c("MATCHES_a_a_a", "MATCHES_b_a_a")
# )
# rownames(analyses) <- c("a_a_a", "b_a_a")
#   
# stub function for parallelization
mrd <- function(i, phylo, RD){
  # take tip labels that are in comm matrix with species occurrences, take the mean of their root distances
  MRD <- mean(RD[which(phylo$tip.label %in% phylo$tip.label[phylo$tip.label %in% colnames(comm)[comm[i,] == 1]])])
  return(MRD)
}


# calculate observed MRD
MRD <- unlist(mclapply(1:nrow(comm), mrd, phylo=testtree, RD=RD.testtree, mc.cores = 4))
# in the testtree there is NaN , that is normal behaviour


# calculate randomized MRDs
for(i in 1:99){
  phylo.rnd <- get(trees)
  phylo.rnd$tip.label <- sample(phylo.rnd$tip.label)
  MRD <- cbind(MRD, unlist(mclapply(1:nrow(comm), mrd, phylo=testtree, RD=RD.testtree, mc.cores = 4)))
  print(paste("replicate", i, "finished"))
}

colnames(MRD) <- c("obs", paste("rnd.", 1:99, sep=""))
rownames(MRD) <- rownames(comm)

# this happends on the server now: resolve_polyotmies.R


# 
# 
# 
# for(n in rownames(analyses)){
# 
#   # choose relevant matching
#   MATCHES <- get(analyses[n,2])
#   phylo   <- get(analyses[n,1])
#   RD      <- get(paste("RD.", analyses[n,1], sep=""))
#   print(paste("Starting calculation of analysis", n))
# 
#   # calculate observed MRD
#   MRD <- unlist(mclapply(1:nrow(comm), mrd, MATCHES=MATCHES, phylo=phylo, RD=RD, mc.cores = 28))
#   print("observed MRD calculated")
# 
#   # calculate randomized MRDs
#   MATCHES.rnd <- MATCHES
#   for(i in 1:99){
#     MATCHES.rnd[!is.na(MATCHES.rnd[,2]),2] <- sample(MATCHES.rnd[!is.na(MATCHES.rnd[,2]),2])
#     MRD <- cbind(MRD, unlist(mclapply(1:nrow(comm), mrd, MATCHES=MATCHES.rnd, phylo=phylo, RD=RD, mc.cores = 28)))
#     print(paste("replicate", i, "finished"))
#   }
# 
#   colnames(MRD) <- c("obs", paste("rnd.", 1:99, sep=""))
#   rownames(MRD) <- rownames(comm)
# 
#   assign(paste("MRD.", n, sep=""), as.data.frame(MRD))
# }
# rm(mrd, MRD, phylo, MATCHES, MATCHES.rnd, RD, i, n)

#######################
# 2. Compute phylosim #
#######################

# for(n in rownames(analyses)){
#   #n = rownames(analyses)[x]
# 
#   MATCHES <- get(analyses[n,2])
#   phylo   <- get(analyses[n,1])
#   EDGES   <- get(paste("EDGES.", analyses[n,1], sep=""))
#   print(paste("Starting calculation of analysis", n))
#   
#   # reduce community matrix to species included in the relevant matching
#   comm_red <- comm[,colnames(comm) %in% MATCHES[,2]]
#   
#   # rename community matrix with phylogeny tip labels
#   comm.obs <- comm_red
#   MATCHES.tmp <- MATCHES[!is.na(MATCHES[,2]),]
#   idx <- as.vector(MATCHES.tmp$tip)
#   names(idx) <- as.vector(MATCHES.tmp[,2])
#   rm(MATCHES.tmp)
#   colnames(comm.obs) <- idx[colnames(comm.obs)]
#   rm(idx)
#   
#   ps_no <- list()
#   ps_bl <- list()
#   ps <- phylosim(phylo, comm.obs, EDGES = EDGES)
#   ps_no[[1]] <- ps[[1]]
#   ps_bl[[1]] <- ps[[2]]
#   rm(ps)
#   print("observed ps calculated")
#  
#   # calculate randomised phylosim
#   
#   for(i in 1:99){
#     
#     # shuffle matching  
#     MATCHES.rnd <- MATCHES
#     MATCHES.rnd[!is.na(MATCHES.rnd[,2]),2] <- sample(MATCHES.rnd[!is.na(MATCHES.rnd[,2]),2])
#     MATCHES.tmp <- MATCHES.rnd[!is.na(MATCHES.rnd[,2]),]
#     idx <- as.vector(MATCHES.tmp$tip)
#     names(idx) <- as.vector(MATCHES.tmp[,2])
#     rm(MATCHES.tmp)
#     # create new (randomized) community matrix
#     comm.rnd <- comm_red
#     colnames(comm.rnd) <- idx[colnames(comm.rnd)]
#     rm(idx)
#     
#     ps <- phylosim(phylo, comm.rnd, EDGES = EDGES)
#     ps_no[[i+1]] <- ps[[1]]
#     ps_bl[[i+1]] <- ps[[2]]
#     print(paste("replicate", i, "finished"))
#     rm(comm.rnd, ps)
#   }
#   
#   rm(i, comm_red, comm.obs, MATCHES, phylo, EDGES)
#   assign(paste("ps_no.", n, sep=""), ps_no)
#   assign(paste("ps_bl.", n, sep=""), ps_bl)
#   #saveRDS(get(paste("ps.", n, sep="")), paste("ps.", n, ".rds", sep=""))
#   rm(ps_no, ps_bl)
# }
# rm(n)
# 
# # for(n in rownames(analyses)){
# #   assign(paste("ps.", n, sep=""), readRDS(paste("ps.", n, ".rds", sep="")))
# # }
# # rm(n)

save.image("phylostruct.RData")

save(comm, MRD.a_a_a, MRD.b_a_a, ps_no.a_a_a, ps_no.b_a_a, ps_bl.a_a_a, ps_bl.b_a_a, file="MRD_ps.RData")
