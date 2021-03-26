# Resolve polytomies and get species-level lineage diversification rate (DR) for every species i as the inverse of its mean equal splits measure (see https://www.nature.com/articles/nature11631#Sec2)

# From supplemnt (https://static-content.springer.com/esm/art%3A10.1038%2Fnature11631/MediaObjects/41586_2012_BFnature11631_MOESM105_ESM.pdf): The ES measure for a focal tip on a rooted bifurcating tree is the sum of the edge lengths from the species i to the root, with each consecutive edge discounted by a factor of ½. 
# For a more general equation for non-bifurcating trees, see Redding et al. (2008). --> On a fully bifurcating tree, if we denote the edges between a PE and the root by λ0 through to λr, the ES measure of a species (Eq. (3)) can be reduced to...
# The inverse of this measure can be seen as a measure of the splitting rate of the path to a tip: species in rapidly-diversifying clades will have short edge lengths shared among many species and low ES values, while isolated species on a tree have no evidence of recent diversification and large ES values.  We term the 1/ES metric for a species its species-level lineage diversification rate, or DR

# https://barnabasdaru.com/2020/12/13/fast-evolutionary-distinctiveness-in-r/ "For very large trees, such calculation is a memory-intensive operation and a bottle-neck for these algorithms. Because of these challenges, we developed a new method in our phyloregion package that speeds up the process significantly to produce results in seconds!"


## Load data
library(ape)
library(castor)
library(phytools)
library(parallel)
library(picante)
library(phyloregion)


# n=20 takes 20 minutes. 4 parallel cores take 70seconds for 4 runs --> 30 minutes for 100 runs
# running on server with 16 cores 1000 reps: 70min


phylo <- read.tree("trees/allmb_matched_added_species_Nov20.tre") 

# keep <- sample(1:length(phylo$tip.label), 20000)
# drop <- 1:length(phylo$tip.label)
# drop <- drop[-keep]
# testtree <- drop.tip(phylo,phylo$tip.label[drop])
# 
# system.time({
#   #test <- evol.distinct(testtree, type="equal.splits")
#   set.seed(1234)
#   temp <- multi2di(testtree)
#   es <- evol.distinct(temp, type="equal.splits")
# })
# 
# system.time({
#   #test <- evol.distinct(testtree, type="equal.splits")
#   set.seed(1234)
#   temp <- multi2di(testtree)
#   es2 <- evol_distinct(temp, type="equal.splits")
# })
# 
# es2 <- as.data.frame(es2)
# plot(es$w, es2$es2)
# 
# system.time({
#   set.seed(1234)
#   temp <- multi2di(phylo)
#   es2 <- evol_distinct(phylo, type="equal.splits")
# })
# 43.400s for one complete multi2di run + equal splits

## GO
# increase numbers when ready to go    
reps <- 1000
res <- matrix(ncol=length(phylo$tip.label), nrow=reps)
system.time({
# parallelize it
es <- function(i, phylo){
  # resolve polytomie and calculate root distance
  temp <- multi2di(phylo)
  #es <- evol.distinct(temp, type="equal.splits")
  es2 <- evol_distinct(temp, type="equal.splits")
  ES <- es2
  #dist.all2 <- get_all_distances_to_root(temp, as_edge_count=TRUE)
  #RD <- dist.all2[1:Ntip(temp)]
  return(ES)
} 

ES <- mclapply(1:reps, es, phylo=phylo, mc.cores = 16)
})
# 173.5s for 4reps on 4 cores, 732s for 20reps on 4 cores

output <- matrix(unlist(ES), ncol = length(ES[[1]]), byrow = TRUE)


# # attach original distances
# dist.all <- get_all_distances_to_root(phylo, as_edge_count=TRUE)
# res <- rbind(output, dist.all[1:Ntip(phylo)])
res <- data.frame(species = names(ES[[2]]),
                  mean.es = colMeans(output), 
                  sd.es = apply(output, 2, sd), 
                  range = as.numeric(diff(apply(output, 2, range))),
                  mean.dr = colMeans(1/output)
                  )

#saveRDS(res, "polytomie_ES_results.rds")

# hist(res$mean.dr) #[res$mean.es<5] # some really high values
# hist(res$sd.es)
# 
# 
res.mrd <- readRDS("data/polytomie_RD_results_Nov20.rds")
res <- readRDS("data/polytomie_ES_results.rds")
hist(res.mrd$mean.rd)
hist(log(res$mean.dr))
cor(res.mrd$mean.rd, res$mean.dr, method="s") # rho=0.36, not linear at all
# ggplot(data=data.frame(rd=rd.test, dr=1/es.test), aes(rd,dr))+
#   geom_point(alpha=0.1)+
#   scale_y_log10()



# # scaleable test 
# keep <- sample(1:length(phylo$tip.label), length(phylo$tip.label))
# drop <- 1:length(phylo$tip.label)
# drop <- drop[-keep]
# temp <- drop.tip(phylo,phylo$tip.label[drop])
# temp <- multi2di(temp)
# #plot(temp)
# es.test <- evol_distinct(temp, type="equal.splits")
# rd.all <- get_all_distances_to_root(temp, as_edge_count=TRUE)
# rd.test <- rd.all[1:Ntip(temp)]
# rd.test # looks fine
# es.test # fine too
# library(ggplot2)
# cor.test(rd.test, 1/es.test, method="s")
# ggplot(data=data.frame(rd=rd.test, dr=1/es.test), aes(rd,dr))+
#   geom_point(alpha=0.1)+
#   scale_y_log10()

# 
# # polytomies error bigger than the diffs in root distance between countries?
# hist(res$mean.rd)
# hist(res$range)
# plot(res$mean.rd~res$range, main="polytomie summary mean~range")
# plot(res$mean.rd~res$sd.rd, main="polytomie summary mean~sd")
# 
# hist(res$sd.rd / res$mean.rd)
# mean(res$sd.rd / res$mean.rd)
# 
# 
# 
# analysis polytomie splits
# hist(apply(output, 2, sd))
# hist(apply(output, 2, range))
# hist(apply(output, 2, mean))
# table(apply(output, 2, sd, na.rm=TRUE)==0) # only 15 tip labels not affected by polytomies
# 
# # Correlation of root distance ranges with the former polytomie root distance
# cor.test(diff(apply(output, 2, range)),tail(res, 1), method="s")
# 
# # Correlation of root distance ranges with the mean new root distance
# cor.test(as.numeric(diff(apply(output, 2, range))), apply(output, 2, mean, na.rm=TRUE), method="s")
# 
# library(scales)
# plot(tail(res, 1), diff(apply(res, 2, range)), col = alpha("black",0.01))
# 
# hist(res[,which.max(diff(apply(res, 2, range)))])
