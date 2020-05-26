# read GBIF complete backbone taxonomy and transform into .rds file
library(rgbif)
#dataset key backbone taxonomy plants: d7dddbf4-2cf0-4f39-9b2a-bb099caae36c
#test <- name_usage(datasetKey = "d7dddbf4-2cf0-4f39-9b2a-bb099caae36c", limit = 10000)$data
#View(test)


sb <- readRDS("SB_tip_labels.rds")
ids <- sb$gbif_id
res <- list()
for(i in 1:1000){
  res[[i]] <- name_usage(ids[i], return="data")
  if(!i%%10)cat(i,"\r")
  # include failsafe
  if(i %in% c(50000, 100000, 150000, 200000)){safeRDS(res, "gbif.rds")}
  if(i==length(ids)){safeRDS(res, "gbif_all.rds")}
}

#manual switches
non.null.res <- res[-59]

# extract stuff to dataframe
taxonID <- sapply(non.null.res, "[[", "key")
scientificName <- sapply(non.null.res, "[[", "scientificName")
authorship <- sapply(non.null.res, "[[", "authorship")
taxonomicStatus <- sapply(non.null.res, "[[", "taxonomicStatus")
rank <- sapply(non.null.res, "[[", "rank")
#kingdom <- sapply(res, "[[", "kingdom")
family <- sapply(non.null.res, "[[", "family")
genus <- sapply(non.null.res, "[[", "genus")
species <- sapply(non.null.res, "[[", "species")

gbif <- data.frame(taxonID, 
                   scientificName,
                   authorship,
                   taxonomicStatus,
                   rank,
                   family,
                   genus,
                   species)

#res[[59]] has no family - how to deal with that
#--> remove manuall?

View(gbif)




# # replace nulls with NAs
# non.null.res <- lapply(res, lapply, function(x)ifelse(is.null(x), NA, x))






library(readr)
test <- read_tsv("data/Taxon.tsv", col_names = TRUE, 
         col_types = cols(
            taxonRemarks = col_skip(),
            originalNameUsageID = col_skip(), 
            datasetID = col_skip(),
            parentNameUsageID = col_skip(),
            canonicalName = col_skip(),
            nameAccordingTo = col_skip(),
            namePublishedIn = col_skip(),
            genericName = col_skip(),
            nomenclaturalStatus = col_skip(),
            phylum = col_skip(),
            class = col_skip(),
            order = col_skip(),
            namePublishedIn = col_character()
            ), 
         n_max=161023
)
saveRDS(test, "test.rds")
rm(test)
test2 <- read_tsv("data/Taxon.tsv",
                 skip=161024, 
                 n_max=10000,
                 head()
)
         # 
         # col_types = NULL,
         # locale = default_locale(), na = c("", "NA"), quoted_na = TRUE,
         # quote = "\"", comment = "", trim_ws = TRUE, skip = 0,
         # n_max = Inf, guess_max = min(1000, n_max),
         # progress = show_progress(), skip_empty_rows = TRUE)


skip <- seq(0, 6500000, 250000)
nline <- seq(250000, 6586626, 250000)
nline[length(nline)] <- 6586626
sub_cols <- c(1, 6, 7, 10, 12, 15, 18, 22, 23)
for(i in 1:length(nline)){
  gbif <- read.csv("data/Taxon.tsv", sep="\t", skip=skip[i], nrows=500000, header=!i>1)
  gbif <- gbif[gbif[,18]=="Plantae",]
  if(i==1){
    gbif_sub <- gbif[,c("taxonID", "scientificName", "scientificNameAuthorship",
                        "taxonomicStatus", "taxonRank", "kingdom", "family", "genus", "specificEpithet")]
    saveRDS(gbif_sub, paste0("gbif_sub", i, ".rds"))
    cols <- names(gbif_sub)
    sub_cols <- which(names(gbif) %in% cols)
    rm(list=c("gbif", "gbif_sub"))
  }else{
    gbif_sub <- gbif[,sub_cols]
    saveRDS(gbif_sub, paste0("gbif_sub", i, ".rds"))
    rm(list=c("gbif", "gbif_sub"))
  }

  print(i)
}


# # access OTT API
# curl -X POST https://api.opentreeoflife.org/v3/taxonomy/taxon_info \
# -H 'content-type:application/json' -d '{"ott_id":515698}'
# 
# 
# library(rotl)
# req <- taxonomy_taxon_info(c(614208, 3574))
# tax_name(req)
# 
# deuterostomes <- tnrs_match_names(names=c("echinodermata", "xenacoelomorpha",
#                                           "chordata", "hemichordata"))
#   
#
# gbif <- read.csv("data/Taxon.tsv", sep="\t", nrow=500000)
# gbif <- gbif[gbif$kingdom=="Plantae",]
# gbif_sub <- gbif[,c("taxonID", "scientificName", "scientificNameAuthorship", 
#                     "taxonomicStatus", "taxonRank", "kingdom", "family", "genus", "specificEpithet")]
# saveRDS(gbif_sub, "gbif_sub1.rds")
# rm(list=c("gbif", "gbif_sub"))
# 
# gbif <- read.csv("data/Taxon.tsv", sep="\t", skip=500000, nrow=1000000)
# gbif <- gbif[gbif$kingdom=="Plantae",]
# gbif_sub <- gbif[,c("taxonID", "scientificName", "scientificNameAuthorship", 
#                     "taxonomicStatus", "taxonRank", "kingdom", "family", "genus", "specificEpithet")]
# saveRDS(gbif_sub, "gbif_sub1.rds")
# rm(list=c("gbif", "gbif_sub"))
# 
# gbif2 <- read.csv("data/Taxon.tsv", sep="\t", skip=500000, nrow=1000000)
# gbif3 <- read.csv("Taxon.tsv", sep="\t", skip=2000000, nrow=2000000)
# gbif4 <- read.csv("Taxon.tsv", sep="\t", skip=3000000, nrow=2000000)
# gbif5 <- read.csv("Taxon.tsv", sep="\t", skip=4000000, nrow=2000000)
# gbif6 <- read.csv("Taxon.tsv", sep="\t", skip=5000000, nrow=2000000)
# gbif7 <- read.csv("Taxon.tsv", sep="\t", skip=6000000, nrow=2000000)
# 
# gbif_sub <- gbif[gbif$kingdom=="Plantae",]
# gbif_sub <- gbif_sub[,c("taxonID", "scientificName", "scientificNameAuthorship", 
#                         "taxonomicStatus", "taxonRank", "kingdom", "family", "genus", "specificEpithet")]
# 
# saveRDS(gbif_sub, "gbif_plantae.rds")
