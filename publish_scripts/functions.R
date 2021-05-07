# function that gets all the matching wcsp entries for a given tip (i) in the phylogeny (phylo)
# requires separate lfists for wcsp species and infraspecific taxa 
# used in: match_data.R
get_matches <- function(i, phylo, wcsp_species, wcsp_infra){
  
  # initialize search for this tip label
  tip <- strsplit(phylo$tip.label[i], "_")[[1]] # split tip label by underscores
  matches <- c() # initialize vector to receive matching checklist ids
  
  # look for infraspecific taxa if enough information is available
  if(length(tip)>=4){
    #look for matches in infraspecific taxa (keeping all matches) and find the corresponding accepted species
    for(match in as.vector(wcsp_infra[wcsp_infra$genus == tip[1] & wcsp_infra$species == tip[2] 
                                      & wcsp_infra$infraspecific_rank == gsub("\\.", "", tip[3]) 
                                      & wcsp_infra$infraspecific_epithet == tip[4],"accepted_name_id"])){
      matches <- c(matches, as.vector(wcsp_species[wcsp_species$genus == as.vector(wcsp_infra[match,"genus"]) & wcsp_species$species == as.vector(wcsp_infra[match,"species"]) & wcsp_species$taxon_status_description == "Accepted", "checklist_id"]))
    }
  }
  
  # if there is no information on infraspecific taxa, or match at infraspecific level, move on to species level search
  if(length(matches)==0){
    # look for matches in species (keeping all matches) and find the corresponding accepted species in cases where a synonym of an ifrasp. taxon is hit
    # e.g. Macaranga schleinitziana (wcs-116515) is synonym of Macaranga involucrata var. involucrata (wcs-116515)
    for(match in as.vector(wcsp_species[wcsp_species$genus == tip[1] & wcsp_species$species == tip[2],"accepted_name_id"])){
      if(match %in% rownames(wcsp_species)){ # if it already is a species, keep match
        matches <- c(matches, match)
      } else { # if it is an infraspecific taxon, find corresponding species
        matches <- c(matches, as.vector(wcsp_species[wcsp_species$genus == as.vector(wcsp_infra[match,"genus"]) & wcsp_species$species == as.vector(wcsp_infra[match,"species"]) & wcsp_species$taxon_status_description == "Accepted","checklist_id"]))
      }
    }
  }
  
  return(unique(matches))
  
}

# function for collapsing a list into vector, removing any list entries of length != 1 (accounting for "")
# used in match_data.R
remove_multiples <- function(MATCHES){
  x <- rep(NA, length(MATCHES))
  for(i in 1:length(MATCHES)){
    if(length(MATCHES[[i]][MATCHES[[i]] != ""]) == 1) x[i] <- MATCHES[[i]][MATCHES[[i]] != ""]
  }
  return(x)
}

# Function for resolving cases in which multiple tips in the tree point to the same accepted wcsp name 
# Used in match_data.R
resolve_multiple <- function(MATCHES, wcsp, phylo, phylogb){
  # remove multiple linkages, i.e. when multiple tips in the tree are assigned to the same accepted name, using following criteria
  # 1) preferably keep tips that have molecular data behind them
  # 2) preferably keep tips that have the same genus name as the species they link to in WCSP - important for tips that have been added based on taxonomy
  # 3) randomly thereafter
  for(n in names(table(MATCHES)[table(MATCHES)>1])){
    tips <- which(MATCHES == n)
    erase <- rep(FALSE,length(tips))
    #first criterion: if any of the tips comes from GenBank, erase all that don't
    if(sum(phylo$tip.label[tips] %in% phylogb$tip.label)>0){
      erase <- !phylo$tip.label[tips] %in% phylogb$tip.label
      #break
    } else {#second criterion (only relevant if choosing among non-genbank tips): if at least one of the tips has the same genus as in wcsp, erase the rest
      #retrieve WCSP genus
      gen <- as.vector(wcsp[n,"genus"])
      #retrieve tip genera
      GEN <- vector("character", length(tips))
      for(tip in 1:length(tips)){
        GEN[tip] <- strsplit(phylo$tip.label[tips],"_")[[tip]][1]
      }
      #if at least on e of the tips has the same genus as in wcsp, erase the rest
      if(gen %in% GEN){
        erase <- !GEN == gen
      }
    }
    if(sum(!erase)>1){#if there are still n>1 items to be kept (i.e. not erased), chose n-1 randomly for erasing
      erase[sample(which(erase==FALSE), size=(sum(!erase)-1))] <- TRUE
    } 
    
    #erase tips 
    MATCHES[tips[erase]] <- NA
  }
  return(MATCHES)
}

# Function for finding all edges that connect a given tip (tip) to the root in a phylogeny (phylo)
# Returns the row numbers of the edges in phylo$edge
# Used in: phylosim.R
get_edges <- function(tip,phylo){
  
  #number of species in tree
  nspp <- Ntip(phylo)
  
  #number of internal branches
  n.int = nspp-2 
  
  #vector to receive edges
  edges <- rep(NA,n.int)
  
  #set starting node to tip
  node <- tip
  
  for(j in 1:n.int){
    edges[j] <- which(phylo$edge[,2]==node)
    node <- phylo$edge[phylo$edge[,2]==node,1]
    if(node == nspp+1) break
  }
  
  edges <- edges[!is.na(edges)]
  
  #returns vector of edges = row numbers of phylo$edge
  return(edges)
}

# Function for calculating phylogenetic Simpson dissimilarity among a set of communities (comm) based on a phylogeny (phylo)
# Specificiations for phylo and comm follow picante standards. 
# Requires that all species in comm (column names) are in phylo$tip.label
# Requires get_edges()
# Returns list of two matrices:
# [[1]] phylogenetic Simpson dissimilarity without branch lengths
# [[2]] phylogenetic Simpson dissimilarity with branch lengths
# Used in: phylosim.R
phylosim <- function(phylo, comm, EDGES = NULL, mc.cores = 8){
  
  if(is.null(EDGES)){
    # get the edges connecting each tip to the root of phylo
    EDGES <- mclapply(1:length(phylo$tip.label), get_edges, phylo=phylo, mc.cores=mc.cores)
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
  rm(a,b,c,i,j)
  
  return(list(ps, ps.bl))
}

# A modified version of bind.tip from phytools, in which the "untangling" of the tree is silenced
# This saves time when mulitple tips need to be bound to the tree. However, this means the tree neeeds to be 
# "untangled" manually after adding the last tip using untangle()
# Used in add_species.R
# Requires phytools
bind.tip2 <- function (tree, tip.label, edge.length = NULL, where = NULL, 
                       position = 0, interactive = FALSE, ...) 
{
  if (!inherits(tree, "phylo")) 
    stop("tree should be an object of class \"phylo\".")
  use.edge.length <- if (is.null(tree$edge.length)) 
    FALSE
  else TRUE
  if (use.edge.length == FALSE) 
    tree <- compute.brlen(tree)
  if (interactive == TRUE) {
    plotTree(tree, ...)
    cat(paste("Click where you would like to bind the tip \"", 
              tip.label, "\"\n", sep = ""))
    flush.console()
    obj <- get.treepos(message = FALSE)
    where <- obj$where
    position <- obj$pos
  }
  else if (is.null(where)) 
    where <- Ntip(tree) + 1
  if (where <= Ntip(tree) && position == 0) {
    pp <- 1e-12
    if (tree$edge.length[which(tree$edge[, 2] == where)] <= 
        1e-12) {
      tree$edge.length[which(tree$edge[, 2] == where)] <- 2e-12
      ff <- TRUE
    }
    else ff <- FALSE
  }
  else pp <- position
  if (is.null(edge.length) && is.ultrametric(tree)) {
    H <- nodeHeights(tree)
    if (where == (Ntip(tree) + 1)) 
      edge.length <- max(H)
    else edge.length <- max(H) - H[tree$edge[, 2] == where, 
                                   2] + position
  }
  tip <- list(edge = matrix(c(2, 1), 1, 2), tip.label = tip.label, 
              edge.length = edge.length, Nnode = 1)
  class(tip) <- "phylo"
  obj <- bind.tree(tree, tip, where = where, position = pp)
  if (where <= Ntip(tree) && position == 0) {
    nn <- obj$edge[which(obj$edge[, 2] == which(obj$tip.label == 
                                                  tip$tip.label)), 1]
    obj$edge.length[which(obj$edge[, 2] == nn)] <- obj$edge.length[which(obj$edge[, 
                                                                                  2] == nn)] + 1e-12
    obj$edge.length[which(obj$edge[, 2] == which(obj$tip.label == 
                                                   tip$tip.label))] <- 0
    obj$edge.length[which(obj$edge[, 2] == which(obj$tip.label == 
                                                   tree$tip.label[where]))] <- 0
  }
  #root.time <- if (!is.null(obj$root.time)) 
  #  obj$root.time
  #else NULL
  #obj <- untangle(obj, "read.tree")
  #if (!is.null(root.time)) 
  #  obj$root.time <- root.time
  #if (interactive) 
  #  plotTree(obj, ...)
  if (!use.edge.length) 
    obj$edge.length <- NULL
  obj
}

# Function that calculates root distance for a given tip i in a tree tree, accounting for polytomies using a vector poly that contains for 
# each node how many additional nodes need to be added due to its polytomy status. 
# Not to be used on its own - used in root.distance
# Requires ape
countnodes <- function(i, tree, polys){
  # initiate number of nodes
  nnodes = 0 #vector("numeric",nspp-1) 
  # set current node to terminal node of interest
  node = i
  for(j in 1:(Ntip(tree)-1)){
    node = tree$edge[tree$edge[,2] == node,1] #find new node
    if(node == (Ntip(tree)+1)) break # break loop once it hits the root node
    nnodes = nnodes + 1 + polys[as.character(node)] #add current node to nodes, plys any nodes implicated by the polytomy status
  }
  return(unname(nnodes))
}

# function that estimates the root distance (RD) for all tips in a phylo object
# polytomies are accounted for assuming a birth-death process and using an empirically estimated function (-2/3 + 2*log(N))
# requires ape
# used in phylostruct.R
root.distance <- function(tree, mc.cores=8){
  # score polytomies in tree and estimate how many nodes need to be added for each one
  polys <- table(tree$edge[,1]) #calculates for each edge how many descendents it has
  polys <- -1.677857 + log(polys) * 2.00427 #transforms this into expected RD
  polys[polys <0] = 0
  # calculate root distances for all tips in parallel
  RDs <- mclapply(1:Ntip(tree), countnodes, tree=tree, polys=polys, mc.cores = mc.cores)
  return(RDs)
}

# Function that calculates non-phylogenetic Simpson dissimilarity based on a community matrix
betasim <- function(comm){
  bs_spp <- matrix(ncol = nrow(comm), nrow = nrow(comm))
  rownames(bs_spp) <- colnames(bs_spp) <- rownames(comm)
  
  for(i in 1:(nrow(bs_spp)-1)){
    for(j in (i+1):ncol(bs_spp)){
      a <- sum(comm[i,] == 1 & comm[j,] == 1)
      b <- sum(comm[i,])
      c <- sum(comm[j,])
      bs_spp[i,j] <- 1 - a/min(b,c)
    }
    print(i)
  }
  
  return(bs_spp)
}

# Function that calculates the ratio between the sum of between cluster dissimilarity and the sum of within cluster dissimilarity
# requires cluster vector and dissimilarity matrix (NOT distance object!)
get_ratio <- function(clust,diss){
  inside = 0
  outside = 0
  for(i in 1:nrow(diss)){
    inside = inside + sum(diss[i,clust == clust[i]], na.rm=T)
    outside = outside + sum(diss[i,clust != clust[i]], na.rm=T)
  }
  return(outside/(inside+outside))
}

# Function to match cluster numbers in nested clusterings, where N(new) = N(reference)+1
match_nested <- function(reference,new){
  a <- c()
  l <- c()
  for(i in 1:max(new)){
    a <- c(a,reference[new == i][1])
    l <- c(l,length(reference[new == i]))
  }
  if(l[a == as.numeric(names(table(a))[table(a)==2])][1] == l[a == as.numeric(names(table(a))[table(a)==2])][2]){ #tie
    a[a == as.numeric(names(table(a))[table(a)==2])][1] <- max(new) #arbitrarily broken
  }else{
    a[a == as.numeric(names(table(a))[table(a)==2]) & l == min(l[a == as.numeric(names(table(a))[table(a)==2])])] <- max(new)
  }
  names(a) <- 1:max(new)
  new <- unname(a[as.character(new)])
  return(new)
}
