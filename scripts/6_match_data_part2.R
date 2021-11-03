## Load taxonomy matcher created GBIF tip labels #################################################### 

matches$elevated_to_species_id[which(is.na(matches$elevated_to_species_id))] <- matches$accepted_plant_name_id[which(is.na(matches$elevated_to_species_id))]
matches <- matches[,c(tax_level, "tip")]

# merge the GBIF accepted name IDS into fin
fin <- merge(fin, matches[,c("tip", tax_level)],
             by.x="tips_mod", by.y="tip", all.x=TRUE) # should be matched by tip labels
rm(matches)






# COMBINE NCBI AND GBIF MATCHES ######################################################################
# combine both accepted plant id columns
col_merge <- grep(tax_level, names(fin))
fin[is.na(fin[,col_merge[1]]),col_merge[1]] <- fin[is.na(fin[,col_merge[1]]),col_merge[2]]
fin <- fin[,-col_merge[2]]
names(fin)[grep(tax_level, names(fin))] <- "accepted_id"

table(is.na(fin$accepted_id))/nrow(fin) # 85.2% matches


table(is.na(fin$accepted_id[fin$source=="ncbi"]))/nrow(fin[fin$source=="ncbi",]) # 85.5%
table(is.na(fin$accepted_id[fin$source=="gbif"]))/nrow(fin[fin$source=="gbif",]) # 85.1%

# should be all accepted, check:
wcp <- readRDS("processed_data/apg_wcp_jul_21_clean.rds")
gbif_labels <- fin[fin$source=="gbif",]
ncbi_labels <- fin[fin$source=="ncbi",]

table(wcp$taxon_status[wcp$plant_name_id %in% gbif_labels$accepted_id])
table(wcp$taxon_status[wcp$plant_name_id %in% ncbi_labels$accepted_id])

# # manual correct Citrus pennivesiculata (724235-az --> 724129-az)
# fin$accepted_id[fin$accepted_id=="724235-az"] <- "724129-az"

wcp_tip_sub <- wcp[which(wcp$plant_name_id %in% fin$accepted_id),]
table(wcp_tip_sub$taxon_status)






if(all(wcp_tip_sub$taxon_status=="Accepted")){print("All taxa have status accepted")
}else{stop("There are non-accepted taxa in the tip labels!")}  


#View(fin[which(!fin$accepted_id %in% wcp$accepted_plant_name_id),])





### REPLACE TIP LABELS and RESOlVE MULTIPLE MATCHES #################################################################

# put wcsp IDs on tree
phylo_org <- phylo
phylo$tip.label <- fin$accepted_id[match(phylo$tip.label, fin$tips)]

## duplicates
length(which(table(phylo$tip.label)>1)) # duplicated tip labels: 21,963

resolve_multiple <- function(MATCHES, wcp, phylo, phylo_org, phylogb){
  # matches = accepted WCSP ID tip labels, phylo=phylogeny with labels replaced, phylo_org=original tip labels
  # remove multiple linkages, i.e. when multiple tips in the tree are assigned to the same accepted name, using following criteria
  # 1) preferably keep tips that have molecular data behind them
  # 2) preferably keep tips that have the same genus name as the species they link to in WCSP - important for tips that have been added based on taxonomy
  # 3) randomly thereafter
  for(n in names(table(MATCHES)[table(MATCHES)>1])){
    counter <- 1
    tips <- which(MATCHES == n)
    erase <- rep(FALSE,length(tips))
    #first criterion: if any of the tips comes from GenBank, erase all that don't
    if(sum(phylo$tip.label[tips] %in% phylogb$tip.label)>0){
      erase <- !phylo$tip.label[tips] %in% phylogb$tip.label
      #break
    } else {#second criterion (only relevant if choosing among non-genbank tips): if at least one of the original tips has the same genus as in wcsp, erase the rest
      #retrieve WCSP genus
      gen <- as.vector(wcp[wcp$plant_name_id==n,"genus"])
      #retrieve tip genera
      GEN <- vector("character", length(tips))
      for(tip in 1:length(tips)){
        
        GEN[tip] <- strsplit(phylo_org$tip.label[tips],"_")[[tip]][1]
      }
      #if at least one of the tips has the same genus as in wcsp, erase the rest
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
  counter <- counter + 1
  if(!counter%%1) print(counter)# cat(counter,"\r")
  return(MATCHES)
}

# 10 MINUTES ###
Sys.time()
(res_multi <- resolve_multiple(phylo$tip.label, wcp, phylo, phylo_org, phylogb_a)) 
table(is.na(res_multi)) # NAs introduced
Sys.time()


# DROPPING THE UNUSED TIP LABELS ######################################################################

res_multi_noNA <- na.omit(res_multi)
clean_tree <- keep.tip(phylo, as.character(res_multi_noNA))


