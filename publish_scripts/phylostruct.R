####################################################
# 1. Computing RD #
####################################################

dist.all <- get_all_distances_to_root(phylo, as_edge_count=TRUE)
res$org.rd <- dist.all[1:Ntip(phylo)]

# Correlation of root distance ranges with the former polytomie root distance
# cor.test(res$range, res$org.rd, method="s")
# cor.test(res$mean.rd, res$org.rd, method="s")
# cor.test(res$sd.rd, res$org.rd, method="s")
# 
# plot(res$range, res$org.rd, col = alpha("black",0.01))
# hist(res$range)
# plot(res$range, res$mean.rd, col = alpha("black",0.01),
#      xlab="Range RD per species",
#      ylab="Mean RD")
# 


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

# calculate observed MRD

Sys.time()
MRD <- unlist(mclapply(1:nrow(comm), mrd, phylo=phylo, RD=RD, mc.cores = n_cores))
MRD.sd <- unlist(mclapply(1:nrow(comm), mrd.sd, phylo=phylo, RD.sd=RD.sd, mc.cores = n_cores))
Sys.time()

mrd.res <- data.frame(level3=rownames(comm), mrd=MRD, mrd_sd=MRD.sd)


# 
# # calculate randomized MRDs 
# 
# Sys.time()
# for(i in 1:99){
#   phylo.rnd <- phylo
#   phylo.rnd$tip.label <- sample(phylo.rnd$tip.label)
#   MRD <- cbind(MRD, unlist(mclapply(1:nrow(comm), mrd, phylo=phylo.rnd, RD=RD, mc.cores = n_cores)))
#   cat(paste("replicate", i, "finished"),"\r")
# }
# Sys.time()
# 
# colnames(MRD) <- c("obs", paste("rnd.", 1:99, sep=""))
# rownames(MRD) <- rownames(comm)


