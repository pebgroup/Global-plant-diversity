# Resolve polytomies


wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str())) 
library(phytools)
library(parallel)
library(castor)
n_cores = 16
reps <- 1000


phylo <- read.tree("../processed_data/allmb_matched_added_species_Sep21_clean.tre")


## GO
# increase numbers when ready to go    
res <- matrix(ncol=length(phylo$tip.label), nrow=reps)
Sys.time()
# parallelize it
rd <- function(i, phylo){
  # resolve polytomie and calculate root distance
  temp <- multi2di(phylo)
  dist.all2 <- get_all_distances_to_root(temp, as_edge_count=TRUE)
  RD <- dist.all2[1:Ntip(temp)]
  return(RD)
} 
RD <- mclapply(1:reps, rd, phylo=phylo, mc.cores = n_cores)
Sys.time()
output <- matrix(unlist(RD), ncol = length(RD[[1]]), byrow = TRUE)

# # attach original distances
# dist.all <- get_all_distances_to_root(phylo, as_edge_count=TRUE)
# res <- rbind(output, dist.all[1:Ntip(phylo)])
res <- data.frame(mean.rd = colMeans(output), 
                  sd.rd = apply(output, 2, sd), 
                  range = as.numeric(diff(apply(output, 2, range)))
                  )

saveRDS(res, "../processed_data/polytomie_RD_results_Sep21_2.rds")
#saveRDS(res, "processed_data/polytomie_RD_results.rds")

