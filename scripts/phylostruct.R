# arguments for manual parallelization on server
#args = commandArgs(trailingOnly=TRUE)
#x = as.numeric(args[1])

rm(list=ls())
library(ape)
library(parallel)
library(castor)
library(phytools)
library(scales)
source("scripts/functions.R")

# retrieve link tables, phylogenies, and wcsp data (from add_species.R)
phylo <- read.tree("trees/allmb_matched_added_species_2.tre") 
#wcsp <- readRDS("data/WCSP_clean.apg.rds")  

# retrieve community matrix (from process_geography.R)
comm <- readRDS("data/comm.rds") # object name = comm

# stats ######
# wcsp_acc <- wcsp[wcsp$taxon_status=="Accepted",]
# table(wcsp_acc$plant_name_id %in% phylo$tip.label)/nrow(wcsp_acc)
# table(wcsp_acc$plant_name_id %in% phylo$tip.label)

#table(colnames(comm) %in% phylo$tip.label)

####################################################
# 1. Computing RD #
####################################################


# #### test tree ####
# keep <- sample(1:length(phylo$tip.label), 100)
# drop <- 1:length(phylo$tip.label)
# drop <- drop[-keep]
# testtree <- drop.tip(phylo,phylo$tip.label[drop])
# #rm(drop, keep)
# 
# plot(testtree, align.tip.label=TRUE, cex=0.6, use.edge.length = FALSE)
# 
# 
# # Get number of edges for each tip and node - remove node counts afterwards
# ## dist counts the edges, which is equivalent to the number of nodes if you include the root.
# dist <- get_all_distances_to_root(testtree, as_edge_count=TRUE)
# plot(testtree, align.tip.label = TRUE)
# nodelabels(frame="none",adj=c(1.1,-0.4))
# 
# # dist[1:Ntip(testtree)] shows tip label edge counts
# 
# # combine root distance and tip label
# dist[1:Ntip(testtree)]
# res <- data.frame(tip.label = testtree$tip.label, rd = dist[1:Ntip(testtree)])
# plot(testtree, align.tip.label = TRUE, tip.color = res$rd)
# RD.testtree <- res
################################################################################




#source("scripts/resolve_polytomies.R")

res <- readRDS("data/polytomie_RD_results.rds")
dist.all <- get_all_distances_to_root(phylo, as_edge_count=TRUE)
res$org.rd <- dist.all[1:Ntip(phylo)]

# Correlation of root distance ranges with the former polytomie root distance
cor.test(res$range, res$org.rd, method="s")
cor.test(res$mean.rd, res$org.rd, method="s")
cor.test(res$sd.rd, res$org.rd, method="s")

plot(res$range, res$org.rd, col = alpha("black",0.01))
hist(res$range)
plot(res$range, res$mean.rd, col = alpha("black",0.01),
     xlab="Range RD per species",
     ylab="Mean RD")



####################
# 2. Calculate MRD #
####################

RD <- res$mean.rd
# stub function for parallelization
mrd <- function(i, phylo, RD){
  # take tip labels that are in comm matrix with species occurrences, take the mean of their root distances
  # subsets the matrix to actual occurrences (=1) and then gets relevant colnames
  MRD <- mean(RD[which(phylo$tip.label %in% phylo$tip.label[phylo$tip.label %in% colnames(comm)[comm[i,] == 1]])])
  return(MRD)
}

# calculate observed MRD

# test
MRD <- unlist(mclapply(1:nrow(comm), mrd, phylo=testtree, RD=RD.testtree$rd, mc.cores = 4))
# row 8 produces NaN
RD.testtree$rd[which(testtree$tip.label %in% testtree$tip.label[testtree$tip.label %in% colnames(comm)[comm[8,] == 1]])]

## for real
Sys.time()
MRD <- unlist(mclapply(1:nrow(comm), mrd, phylo=phylo, RD=RD, mc.cores = 4))
Sys.time()

# calculate randomized MRDs 

for(i in 1:99){
  phylo.rnd <- phylo
  phylo.rnd$tip.label <- sample(phylo.rnd$tip.label)
  MRD <- cbind(MRD, unlist(mclapply(1:nrow(comm), mrd, phylo=phylo.rnd, RD=RD, mc.cores = 4)))
  print(paste("replicate", i, "finished"))
}

colnames(MRD) <- c("obs", paste("rnd.", 1:99, sep=""))
rownames(MRD) <- rownames(comm)

saveRDS(MRD, "data/mrd.rds")






#save.image("phylostruct.RData")
#save(comm, MRD.a_a_a, MRD.b_a_a, ps_no.a_a_a, ps_no.b_a_a, ps_bl.a_a_a, ps_bl.b_a_a, file="MRD_ps.RData")















# this takes >1h

# for(i in 1:length(trees)){
#   print(Sys.time())
#   assign(paste0("RD.", trees[i]), unlist(root.distance(get(trees[i]), mc.cores = 8)))
#   #saveRDS(get(paste0("RD.", trees[i])), paste0("RD.", trees[i], ".rds"))
#   print(paste(paste0("RD.", trees[i]), "calculated.", Sys.time()))
# 
#   # assign(paste("EDGES.", trees[i], sep=""), mclapply(1:Ntip(get(trees[i])), get_edges, phylo=get(trees[i]), mc.cores=8))
#   # #saveRDS(get(paste("EDGES.", trees[i], sep="")), paste("EDGES.", trees[i], ".rds", sep=""))
#   # print(paste(paste("EDGES.", trees[i], sep=""), "calculated.", Sys.time()))
# }
# rm(i)
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

