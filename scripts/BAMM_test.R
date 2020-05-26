library(ape)
library(phangorn)
v <- read.tree("trees/allmb_matched_added_species.tre")
#tree <- v

# get a meaningful subsample of the tree: use palms
wcsp <- readRDS("data/WCSP.apg.rds")
palms <- wcsp$plant_name_id[wcsp$family=="Arecaceae" & wcsp$taxon_status=="Accepted"]
keeps <- v$tip.label[which(v$tip.label %in% palms)]
library(phytools)
n <- findMRCA(v, tips = keeps, type="node")



tree <- extract.clade(v, n)
plot(tree, show.tip.label = FALSE)

# make ultrametric (http://blog.phytools.org/2017/03/forceultrametric-method-for-ultrametric.html)
is.ultrametric(tree)
force.ultrametric<-function(tree,method=c("nnls","extend")){
  method<-method[1]
  if(method=="nnls") tree<-nnls.tree(cophenetic(tree),tree,
                                     rooted=TRUE,trace=0)
  else if(method=="extend"){
    h<-diag(vcv(tree))
    d<-max(h)-h
    ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
               y=tree$edge[,2])
    tree$edge.length[ii]<-tree$edge.length[ii]+d
  } else 
    cat("method not recognized: returning input tree\n\n")
  tree
}
tree <- force.ultrametric(tree)
is.ultrametric(tree)

# binary?
is.binary.tree(tree)
tree <- multi2di(tree)

## Now to check min branch length: must be greater than 0!
min(tree$edge.length)
  which(tree$edge.length==0)
sort(tree$edge.length)[1:10]

# quick fix:
tree$edge.length[which(tree$edge.length==0)] <- 0.0000001
is.binary.tree(tree)
min(tree$edge.length)


# tip labels need to be character values, also no dots etc
tree$tip.label <- sub("\\.", replacement = "", tree$tip.label)
write.tree(tree, "scripts/BAMM/BAMM_test.tre")


# get priors
library(BAMMtools)
setBAMMpriors(tree)

# write into control.txt
#run BAMM: bamm -c diversification_control.txt



# read results ##################################################3
 # 100000 generations, 6000 tip labels (~5 min  )
mcmcout <- read.csv("scripts/BAMM/mcmc_out.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

library(coda)
effectiveSize(postburn$N_shifts)
effectiveSize(postburn$logLik)
# both should be 200

edata <- getEventData(read.tree("scripts/BAMM/BAMM_test.tre"), eventdata = "scripts/BAMM/event_data.txt", burnin=0.1)

# speciation rates: edata$meanTipLambda
hist(edata$meanTipLambda)


plotRateThroughTime(edata, ratetype="speciation")

cols <- round(edata$meanTipLambda,2)
cols2 <- cols
cols2[cols==0.07] <- "blue"
cols2[cols==0.08] <- "green"
cols2[cols==0.09] <- "yellow"
cols2[cols==0.10] <- "orange"
cols2[cols==0.24] <- "red"
plot(tree, edge.color = cols2, "fan", show.tip.label = FALSE)

#You can also estimate individual tip-specific rates. For the whale example, this is actually   included as part of your bammdata object. If edata is the bammdata object for whales, the components edata$meanTipLambda and edata$meanTipMu are the relevant model-averaged mean rates of speciation and extinction at the tips of the tree.


