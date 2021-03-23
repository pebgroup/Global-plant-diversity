# Resolve polytomies

## Load data
library(ape)
library(castor)
library(phytools)
library(parallel)

# n=20 takes 20 minutes. 4 parallel cores take 70seconds for 4 runs --> 30 minutes for 100 runs
# running on server with 16 cores 1000 reps: 70min


phylo <- read.tree("trees/allmb_matched_added_species_Nov20.tre") 


## GO
# increase numbers when ready to go    
reps <- 2
res <- matrix(ncol=length(phylo$tip.label), nrow=reps)
Sys.time()
# for(i in 1:reps){
#   # phytools
#   temp <- multi2di(phylo)
#   dist.all2 <- get_all_distances_to_root(temp, as_edge_count=TRUE)
#   res[i,] <- dist.all2[1:Ntip(temp)]
#   if(!i%%1)cat(i,"\r")
# }
# parallelize it
rd <- function(i, phylo){
  # resolve polytomie and calculate root distance
  temp <- multi2di(phylo)
  dist.all2 <- get_all_distances_to_root(temp, as_edge_count=TRUE)
  RD <- dist.all2[1:Ntip(temp)]
  return(RD)
} 
RD <- mclapply(1:reps, rd, phylo=phylo, mc.cores = 4)
Sys.time()
output <- matrix(unlist(RD), ncol = length(RD[[1]]), byrow = TRUE)

# # attach original distances
# dist.all <- get_all_distances_to_root(phylo, as_edge_count=TRUE)
# res <- rbind(output, dist.all[1:Ntip(phylo)])
res <- data.frame(mean.rd = colMeans(output), 
                  sd.rd = apply(output, 2, sd), 
                  range = as.numeric(diff(apply(output, 2, range)))
                  )


saveRDS(res, "polytomie_RD_results.rds")
#res <- readRDS("data/polytomie_RD_results.rds")

# polytomies error bigger than the diffs in root distance between countries?
hist(res$mean.rd)
hist(res$range)
plot(res$mean.rd~res$range, main="polytomie summary mean~range")
plot(res$mean.rd~res$sd.rd, main="polytomie summary mean~sd")

hist(res$sd.rd / res$mean.rd)
mean(res$sd.rd / res$mean.rd)



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
