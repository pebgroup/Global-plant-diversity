rm(list=ls())
library("tidyverse")
library("taxonlookup")
library("stringr")

f.apg <- read.csv("./data/apgweb_parsed.csv", stringsAsFactors = F)
WCSP <- readRDS("./data/wcp_dec_19.rds")

# remove clade names without valid names ("near_XXX") or a fossil name ("fossil")
f.apg <- f.apg %>% 
  filter(!(str_detect(Clade, "near_") | str_detect(Clade, "fossil"))) %>% 
  select(Syn_Fam, Acc_Fam, Clade)

# manual edits 
f.apg[grep("Adoxaceae", f.apg$Syn_Fam),]$Acc_Fam <- "Adoxaceae"
f.apg[grep("Viburnaceae", f.apg$Syn_Fam),]$Acc_Fam <- "Adoxaceae"
f.apg <- rbind.data.frame(f.apg, c("Rhipogonaceae", "Ripogonaceae", "Liliales"),
                          c("Ripogonaceae", "Ripogonaceae", "Liliales")) %>% 
  arrange(Syn_Fam) %>% 
  unique()

names(f.apg) <- c("family", "family.apg", "order")


  

#check for no match and fix

setdiff(unique(WCSP$family), unique(f.apg$family))
# [1] "Aspleniaceae"   "Byxaceae"       "Gigaspermaceae" "Incertae_sedis" "Isoetaceae"    
# [6] "Oligomeris"     "Osmundaceae"    "Polypodiaceae"  "Schoberia"      "v" 

########################################################################
# caution serval invalid strings in fern clade names (check before do)#
########################################################################
# #ferns families
# "Aspleniaceae"
# "Osmundaceae" 
# "Polypodiaceae"
# "Isoetaceae"
# "v"=="Ophioglossum"="1142939-az"
# "Ophioglossaceae" 
# "Schizaeaceae"
# #mosse remove
# "Gigaspermaceae"
#Incertae_sedis
# > unique(WCSP[grep("Incertae_sedis", WCSP$family),]$genus)
# [1] "Angeja"      "Anonymos"    "Cipum"       "Euphrona"    "Ivonia"      "Pouslowia"  
# [7] "Theodoricea" "Thuraria"    "Urceola"

# remove mosses - uncomment this line below if want to keep ferns
moss.inc <- c("Gigaspermaceae", "Incertae_sedis")
WCSP1 <- WCSP %>% 
  filter(!(family %in% moss.inc)) %>% 
  filter(accepted_plant_name_id!="1142939-az") # what is special about this fern species? family=Schizaeaceae

#comment line below if want to keep ferns and uncomment the lines above
# ferns.moss.inc <- c("Aspleniaceae", "Osmundaceae", "Polypodiaceae", "Isoetaceae", "Ophioglossaceae", "Schizaeaceae", "Gigaspermaceae", "Incertae_sedis")
# WCSP1 <- WCSP %>% filter(!(family %in% ferns.moss.inc)) %>% filter(accepted_plant_name_id!="1142939-az")

setdiff(unique(WCSP1$family), unique(f.apg$family))
#rename
# "Oligomeris" "Resedaceae", "1012613-az," Genus name listed as family
# "Schoberia" "Amaranthaceae", "1036695-az", Genus name listes as family
# "Byxaceae" is actually "Buxaceae"
# Isoetaceae are listed as Isoëtaceae (special character) in f.apg
# Osmundaceae is a fern family, not listed in APG


# WCSP1$family <- mgsub::mgsub(WCSP1$family, c("Byxaceae", "Oligomeris", "Schoberia"), c("Buxaceae", "Resedaceae", "Amaranthaceae"))
# correct famlies names
WCSP1$family <- gsub("Byxaceae", "Buxaceae", WCSP1$family)
WCSP1$family <- gsub("Oligomeris", "Resedaceae", WCSP1$family)
WCSP1$family <- gsub("Schoberia", "Amaranthaceae", WCSP1$family)

f.apg$family <- gsub("Isoëtaceae", "Isoetaceae", f.apg$family)
f.apg$family.apg <- gsub("Isoëtaceae", "Isoetaceae", f.apg$family.apg)
# WCSP2 <- WCSP1 %>% filter(accepted_plant_name_id!="1142939-az") # ? 

setdiff(unique(WCSP1$family), unique(f.apg$family))


WCSP.apg <- left_join(WCSP1, f.apg)

saveRDS(WCSP.apg, "./data/WCSP.apg.rds")
