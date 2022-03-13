# Get diversification metrics

wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str())) 
library(ape)
library(castor)
library(parallel)
n_cores=4

res <- readRDS("../processed_data/polytomie_RD_results_Sep21_2.rds")
temp <- readRDS("../processed_data/DR_tip_rates.rds")
phylo <- read.tree("../processed_data/allmb_matched_added_species_Sep21_clean.tre") 
comm <- readRDS("../processed_data/comm_September2021.rds")

# compare phylogenetic with geographic data again
table(colnames(comm) %in% phylo$tip.label) # should be all represented



# Root distance

dist.all <- get_all_distances_to_root(phylo, as_edge_count=TRUE)
res$org.rd <- dist.all[1:Ntip(phylo)]



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
saveRDS(mrd.res, "../processed_data/mrd_Sep2021.rds")


# Calculate mean DR --------------------------------------------------------

row.names(temp) <- temp$plant_id
all(row.names(temp) %in% colnames(comm))
all(colnames(comm) %in% row.names(temp))
comm <- comm[,-which(! colnames(comm) %in% row.names(temp))]
# function for parallelization
av.per.tdwg <- function(i, data, comm.matrix){
  # data = df with metrics in columns and row names = species ID
  # returns a named numeric vector
  presences <- which(row.names(data) %in% colnames(comm.matrix)[comm.matrix[i,]==1])
  tmp <- apply(data[presences,], 2, mean, na.rm=TRUE)
  #  tmp <- mean(data[which(data[,colns] %in% data[presences,]),coln])
  return(tmp)
}


temp <- temp[,-grep("plant_id|ES\\.", names(temp))]
#rm(res, res.df, sorted_dfl)

n_cores=2
system.time(
  test <- mclapply(1:nrow(comm), av.per.tdwg, data=temp,
                   comm.matrix=comm, mc.cores=n_cores)
)
res <- data.frame(
  MDR = sapply(test, "[[", "MDR"),
  level3 = row.names(comm))

saveRDS(res, "../processed_data/DR_rates_Feb2022.rds")


