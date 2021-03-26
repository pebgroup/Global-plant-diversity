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
phylo <- read.tree("trees/allmb_matched_added_species_Nov20.tre") 

# retrieve community matrix (from process_geography.R)
comm <- readRDS("data/comm_Feb2021.rds") # object name = comm # use new one

# stats ######
# wcsp_acc <- wcsp[wcsp$taxon_status=="Accepted",]
# table(wcsp_acc$plant_name_id %in% phylo$tip.label)/nrow(wcsp_acc)
# table(wcsp_acc$plant_name_id %in% phylo$tip.label)

table(colnames(comm) %in% phylo$tip.label)
colnames(comm)[!colnames(comm) %in% phylo$tip.label]
# "385449-wcs" is not in the tree: Lorales - Gomortegaceae - Gormortega: monotypic family, not included in tree is according to protocol

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

res <- readRDS("data/polytomie_RD_results_Nov20.rds")
res.es <- readRDS("data/polytomie_ES_results.rds")
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
RD.sd <- res$sd.rd
# stub function for parallelization
mrd.sd <- function(i, phylo, RD.sd){
  # take tip labels that are in comm matrix with species occurrences, take the mean of their root distances
  # subsets the matrix to actual occurrences (=1) and then gets relevant colnames
  MRD.sd <- mean(RD.sd[which(phylo$tip.label %in% phylo$tip.label[phylo$tip.label %in% colnames(comm)[comm[i,] == 1]])])
  return(MRD.sd)
}

DR <- res.es$mean.dr
# stub function for parallelization
mdr <- function(i, phylo, DR){
  # take tip labels that are in comm matrix with species occurrences, take the mean of their root distances
  # subsets the matrix to actual occurrences (=1) and then gets relevant colnames
  MDR <- mean(DR[which(phylo$tip.label %in% phylo$tip.label[phylo$tip.label %in% colnames(comm)[comm[i,] == 1]])])
  return(MDR)
}
DR.sd <- res.es$sd.es
# stub function for parallelization
mdr.sd <- function(i, phylo, DR.sd){
  # take tip labels that are in comm matrix with species occurrences, take the mean of their root distances
  # subsets the matrix to actual occurrences (=1) and then gets relevant colnames
  MDR.sd <- mean(DR.sd[which(phylo$tip.label %in% phylo$tip.label[phylo$tip.label %in% colnames(comm)[comm[i,] == 1]])])
  return(MDR.sd)
}
# calculate observed MRD

# test
#MRD <- unlist(mclapply(1:nrow(comm), mrd, phylo=testtree, RD=RD.testtree$rd, mc.cores = 4))
# row 8 produces NaN
#RD.testtree$rd[which(testtree$tip.label %in% testtree$tip.label[testtree$tip.label %in% colnames(comm)[comm[8,] == 1]])]

## for real
Sys.time()
MRD <- unlist(mclapply(1:nrow(comm), mrd, phylo=phylo, RD=RD, mc.cores = 4))
MRD.sd <- unlist(mclapply(1:nrow(comm), mrd.sd, phylo=phylo, RD.sd=RD.sd, mc.cores = 4))
MDR <- unlist(mclapply(1:nrow(comm), mdr, phylo=phylo, DR=DR, mc.cores = 4))
MDR.sd <- unlist(mclapply(1:nrow(comm), mdr.sd, phylo=phylo, DR.sd=DR.sd, mc.cores = 4))
Sys.time()

saveRDS(MRD.sd, "results/MRD_standard_deviation.rds")
saveRDS(MDR.sd, "results/MDR_standard_deviation.rds")
# calculate randomized MRDs 

Sys.time()
for(i in 1:99){
  phylo.rnd <- phylo
  phylo.rnd$tip.label <- sample(phylo.rnd$tip.label)
  #MRD <- cbind(MRD, unlist(mclapply(1:nrow(comm), mrd, phylo=phylo.rnd, RD=RD, mc.cores = 4)))
  MDR <- cbind(MDR, unlist(mclapply(1:nrow(comm), mdr, phylo=phylo.rnd, DR=DR, mc.cores = 4)))
  #print()
  cat(paste("replicate", i, "finished"),"\r")
}
Sys.time()

# colnames(MRD) <- c("obs", paste("rnd.", 1:99, sep=""))
# rownames(MRD) <- rownames(comm)
colnames(MDR) <- c("obs", paste("rnd.", 1:99, sep=""))
rownames(MDR) <- rownames(comm)

#saveRDS(MRD, "data/mrd_Nov2020.rds")
saveRDS(MDR, "data/mdr_Mar2021.rds")


# compare to old mrd that included ferns:
mrd_old <- readRDS("data/mrd_Nov2020.rds")
mrd_new <- readRDS("data/mrd_feb2021.rds")
plot(mrd_new[,1], mrd_old[,1])
cor(na.omit(mrd_new[,1]), na.omit(mrd_old[,1]))

# compare MRD with DR metric
mrd_new <- readRDS("data/mrd_feb2021.rds")
mdr <- readRDS("data/mdr_Mar2021.rds")

par(mfrow=c(2,2))
plot(mrd_new[,1], mdr[,1], xlab="MRD", ylab="1/ES")
hist(mdr[,1], breaks=20, xlab="1/ES", main="")
hist(mrd_new[,1], breaks=20, xlab="MRD", main="")
plot(res.es$sd.es, res$sd.rd)
#plot(res.es$sd.es/res.es$mean.es, res$sd.rd/res$mean.rd)
dev.off()
#save.image("phylostruct.RData")
#save(comm, MRD.a_a_a, MRD.b_a_a, ps_no.a_a_a, ps_no.b_a_a, psf_bl.a_a_a, ps_bl.b_a_a, file="MRD_ps.RData")















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

# my interpret #####
# subset community matrix to species in phylogeny
comm_sub <- comm[,which(colnames(comm) %in% phylo$tip.label)]
all(colnames(comm_sub) %in% phylo$tip.label)
tokeep <- phylo$tip.label[which(phylo$tip.label %in% colnames(comm_sub))]
phylo <- keep.tip(phylo, tokeep)

ps_no <- list()
ps_bl <- list()
system.time(
ps <- phylosim(phylo, comm_sub, mc.cores=4) # 
)

if(is.null(EDGES)){
  # get the edges connecting each tip to the root of phylo
  EDGES <- mclapply(1:length(phylo$tip.label), get_edges, phylo=phylo, mc.cores=4)
}

# create list of vectors with edges connecting all tips to the root for each community
comm_edges <- list()
for(i in 1:nrow(comm)){
  #extract all edges for species in community i
  spp <- colnames(comm)[comm[i,]==1]
  unique(unlist(EDGES[phylo$tip.label %in% spp])) -> comm_edges[[i]]
}
rm(spp,i)

# create distance matrix
ps <- matrix(nrow=nrow(comm), ncol=nrow(comm))
rownames(ps) <- colnames(ps) <- rownames(comm)
ps[diag(ps)] <- 0

# create distance matrix for branch lengths
ps.bl <- ps

# populate distance matrix using formula by Holt et al
for(i in 1:(nrow(ps)-1)){
  for(j in (i+1):ncol(ps)){
    shared <- intersect(comm_edges[[i]], comm_edges[[j]])
    # calculate without branch lengths
    a <- length(shared)
    b <- length(comm_edges[[i]])# - a
    c <- length(comm_edges[[j]])# - a
    ps[i,j] <- 1 - a/min(b,c)# + a)
    # calculate with branch lengths
    A <- sum(phylo$edge.length[shared])
    B <- sum(phylo$edge.length[comm_edges[[i]]])# - A
    C <- sum(phylo$edge.length[comm_edges[[j]]])# - A
    ps.bl[i,j] <- 1 - A/min(B,C)# + A)
  }
}


######### end check


ps_no[[1]] <- ps[[1]]
ps_bl[[1]] <- ps[[2]]
rm(ps)
print("observed ps calculated")

# calculate randomised phylosim
for(i in 1:99){
  # create new (randomized) community matrix by shuffling colnames
  comm.rnd <- comm_sub
  idx <- colnames(comm.rnd)
  colnames(comm.rnd) <- sample(idx)
  
  ps <- phylosim(phylo, comm.rnd, mc.cores = 4)
  # Returns list of two matrices:
  # [[1]] phylogenetic Simpson dissimilarity without branch lengths
  # [[2]] phylogenetic Simpson dissimilarity with branch lengths
  ps_no[[i+1]] <- ps[[1]]
  ps_bl[[i+1]] <- ps[[2]]
  print(paste("replicate", i, "finished"))
  rm(comm.rnd, ps)
}

rm(i, comm_red, comm.obs, MATCHES, phylo, EDGES)
assign(paste("ps_no.", n, sep=""), ps_no)
assign(paste("ps_bl.", n, sep=""), ps_bl)
#saveRDS(get(paste("ps.", n, sep="")), paste("ps.", n, ".rds", sep=""))
rm(ps_no, ps_bl)
}
rm(n)

# for(n in rownames(analyses)){
#   assign(paste("ps.", n, sep=""), readRDS(paste("ps.", n, ".rds", sep="")))
# }
# rm(n)
# my interpret END #####







for(n in rownames(analyses)){
  #n = rownames(analyses)[x]

  MATCHES <- get(analyses[n,2])
  phylo   <- get(analyses[n,1])
  EDGES   <- get(paste("EDGES.", analyses[n,1], sep=""))
  print(paste("Starting calculation of analysis", n))

  # reduce community matrix to species included in the relevant matching
  comm_red <- comm[,colnames(comm) %in% MATCHES[,2]]

  # rename community matrix with phylogeny tip labels
  comm.obs <- comm_red
  MATCHES.tmp <- MATCHES[!is.na(MATCHES[,2]),]
  idx <- as.vector(MATCHES.tmp$tip)
  names(idx) <- as.vector(MATCHES.tmp[,2])
  rm(MATCHES.tmp)
  colnames(comm.obs) <- idx[colnames(comm.obs)]
  rm(idx)

  ps_no <- list()
  ps_bl <- list()
  ps <- phylosim(phylo, comm.obs, EDGES = EDGES)
  ps_no[[1]] <- ps[[1]]
  ps_bl[[1]] <- ps[[2]]
  rm(ps)
  print("observed ps calculated")

  # calculate randomised phylosim

  for(i in 1:99){

    # shuffle matching
    MATCHES.rnd <- MATCHES
    MATCHES.rnd[!is.na(MATCHES.rnd[,2]),2] <- sample(MATCHES.rnd[!is.na(MATCHES.rnd[,2]),2])
    MATCHES.tmp <- MATCHES.rnd[!is.na(MATCHES.rnd[,2]),]
    idx <- as.vector(MATCHES.tmp$tip)
    names(idx) <- as.vector(MATCHES.tmp[,2])
    rm(MATCHES.tmp)
    # create new (randomized) community matrix
    comm.rnd <- comm_red
    colnames(comm.rnd) <- idx[colnames(comm.rnd)]
    rm(idx)

    ps <- phylosim(phylo, comm.rnd, EDGES = EDGES)
    ps_no[[i+1]] <- ps[[1]]
    ps_bl[[i+1]] <- ps[[2]]
    print(paste("replicate", i, "finished"))
    rm(comm.rnd, ps)
  }

  rm(i, comm_red, comm.obs, MATCHES, phylo, EDGES)
  assign(paste("ps_no.", n, sep=""), ps_no)
  assign(paste("ps_bl.", n, sep=""), ps_bl)
  #saveRDS(get(paste("ps.", n, sep="")), paste("ps.", n, ".rds", sep=""))
  rm(ps_no, ps_bl)
}
rm(n)

# for(n in rownames(analyses)){
#   assign(paste("ps.", n, sep=""), readRDS(paste("ps.", n, ".rds", sep="")))
# }
# rm(n)

