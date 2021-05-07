# # read gbif server
# gbif <- read.csv("Taxon.tsv", sep="\t")
# gbif <- gbif[gbif$kingdom=="Plantae",]
# gbif_sub <- gbif[,c("taxonID", "scientificName", "scientificNameAuthorship", 
#                    "taxonomicStatus", "taxonRank", "kingdom", "family", "genus", "specificEpithet")]
# saveRDS(gbif_sub, "gbif_sub.rds")

    # read GBIF complete backbone taxonomy and transform into .rds file
    library(rgbif)
    #dataset key backbone taxonomy plants: d7dddbf4-2cf0-4f39-9b2a-bb099caae36c
    #test <- name_usage(datasetKey = "d7dddbf4-2cf0-4f39-9b2a-bb099caae36c", limit = 10000)$data
    #View(test)
    
    # DEFINE HERE TAXA YOU WANT TO DOWNLOAD (via gbif ID)
    sb <- read.csv("../Teaching/jeppe/taxonKeys_jeppe.csv")
    ids <- sb$taxonID
    #sb <- readRDS("SB_tip_labels.rds")
    #ids <- sb$gbif_id
    res <- list()
    for(i in 1:length(ids)){
      tryCatch({q
      res[[i]] <- name_usage(ids[i], return="data")
      if(!i%%10)cat(i,"\r")
      # include failsafe
      if(i %in% c(50000, 100000, 150000, 200000)){saveRDS(res, "gbif.rds")}
      if(i==length(ids)){saveRDS(res, "gbif_jeppe.rds")}
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }