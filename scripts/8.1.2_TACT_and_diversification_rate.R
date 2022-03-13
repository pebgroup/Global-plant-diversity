

# PROCESS TACTed TREES #####################################################

wd <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
rm(list = setdiff(ls(), lsf.str())) 
library(ape)
library(phyloregion)
library(parallel)
library(picante)
library(castor)
library(data.table)


nt <- grep("*newick\\.tre", dir("../processed_data/tact/"))
trees <- dir("../processed_data/tact/", full.names = T)[nt]
length(trees)
myfun <- function(x){return(read_tree(file=x))}
alltrees=lapply(trees, myfun)

# polish tip labels
alltrees <- lapply(alltrees, function(x){
  x$tip.label <- gsub("'", "", x$tip.label)
  return(x)}
  )

# load presence absence matrix
#comm <- readRDS("../processed_data/comm_September2021.rds")
#comm <- comm[-which(!colnames(comm)%in%alltrees[[1]]$tip.label),]

# calculate equal splits
system.time(
  res <- lapply(alltrees, evol_distinct, type = "equal.splits")
)

rm(alltrees)


# Get mean per tip --------------------------------------------------------

# sort list for plant ID
res.df <- lapply(res, as.data.frame)
for(i in 1:length(res.df)){
  res.df[[i]][,2] <- row.names(res.df[[i]])
  colnames(res.df[[i]]) <- c("ES", "plant_id")
}
sorted_dfl <- lapply(res.df, function(df){
  df[order(df$plant_id),]
})
i=1;while(all(sorted_dfl[[i]]$plant_id==sorted_dfl[[i+1]]$plant_id) &
          i<=length(trees)-2){i <- i+1;if(i==length(trees)-2)print(
            "all data frames are in correct order (following plant_id)")}  

temp <- as.data.table(sorted_dfl, key="plant_id")
temp <- as.data.frame(temp)
temp <- temp[,grep("plant_id$|ES", names(temp))]
temp$MES <- rowMeans(temp[,grep("^ES", names(temp))])
temp$MDR <- apply(temp[,grep("^ES", names(temp))], 1, function(x)mean(x^(-1)))

saveRDS(temp, "../processed_data/DR_tip_rates.rds")



