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
